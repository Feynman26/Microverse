extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_extracellular_catalog_and_world_registration()
	_test_secondary_fields_begin_empty_and_diffuse_conservatively()
	_test_secondary_fields_do_not_grant_unscripted_uptake()
	_test_lysis_conserves_modeled_structural_material()
	_test_lysis_release_reaches_world_and_event_log()
	_test_same_state_replays_lysis_exactly()

	if failures == 0:
		print("PASS: %d M7 extracellular-substrate tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M7 extracellular-substrate tests failed" % [failures, tests_run])
		quit(1)

func _empty_environment_config():
	var config = SimConfigScript.new()
	config.world_width = 12
	config.world_height = 12
	config.initial_glucose = 0.0
	config.initial_oxygen = 0.0
	config.initial_nitrogen = 0.0
	config.initial_phosphorus = 0.0
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		config.secondary_extracellular_initial[metabolite_id] = 0.0
	return config

func _test_extracellular_catalog_and_world_registration() -> void:
	var config = SimConfigScript.new()
	var sim = SimulationEngineScript.new(config)
	var extracellular_ids: Array[String] = MetaboliteCatalogScript.extracellular_ids()
	var field_names: Dictionary = {}
	for metabolite_id in extracellular_ids:
		var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
		field_names[field_name] = true
		_assert_true(sim.world.has_field(field_name), "world registers extracellular field for %s" % metabolite_id)
		_assert_true(MetaboliteCatalogScript.extracellular_metabolite_for_field(field_name) == metabolite_id, "extracellular field mapping round-trips for %s" % metabolite_id)

	_assert_true(extracellular_ids.size() == 14, "M7-A exposes fourteen chemically explicit extracellular metabolites")
	_assert_true(field_names.size() == extracellular_ids.size(), "every extracellular metabolite has a unique chamber field")
	_assert_true(sim.world.fields.size() == extracellular_ids.size(), "world contains no unregistered ad-hoc chemical fields")
	_assert_true(MetaboliteCatalogScript.extracellular_field("G") == "glucose", "legacy glucose field identity remains stable")
	_assert_true(MetaboliteCatalogScript.extracellular_field("O2") == "oxygen", "legacy oxygen field identity remains stable")
	_assert_true(MetaboliteCatalogScript.extracellular_field("NH4") == "nitrogen", "legacy nitrogen field identity remains stable")
	_assert_true(MetaboliteCatalogScript.extracellular_field("P") == "phosphorus", "legacy phosphorus field identity remains stable")
	_assert_true(not MetaboliteCatalogScript.has_extracellular_field("ATP"), "energy currency is not silently turned into a diffusible nutrient")
	_assert_true(not MetaboliteCatalogScript.has_extracellular_field("BIO"), "structural biomass is decomposed on lysis instead of becoming a magic extracellular BIO field")

func _test_secondary_fields_begin_empty_and_diffuse_conservatively() -> void:
	var config = _empty_environment_config()
	var sim = SimulationEngineScript.new(config)
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
		_assert_close(sim.world.get_field(field_name).total_amount(), 0.0, 1e-12, "secondary field %s starts empty by default" % field_name)

	var field = sim.world.get_field(MetaboliteCatalogScript.extracellular_field("W1"))
	field.add_amount(6, 6, 10.0)
	var before: float = field.total_amount()
	sim.world.diffuse(config.tick_dt_min)
	var after: float = field.total_amount()
	_assert_close(after, before, 1e-10, "secondary extracellular diffusion conserves closed-chamber material")
	_assert_true(field.minimum_value() >= -1e-12, "secondary extracellular diffusion remains nonnegative")
	_assert_true(field.get_value(5, 6) > 0.0 and field.get_value(6, 5) > 0.0, "secondary metabolite spreads through ordinary local diffusion")

func _test_secondary_fields_do_not_grant_unscripted_uptake() -> void:
	var config = _empty_environment_config()
	config.mutation_enabled = false
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	cell.set_pool("W1", 0.0)
	var field_name: String = MetaboliteCatalogScript.extracellular_field("W1")
	sim.world.release(field_name, cell.position, 25.0)
	var extracellular_before: float = sim.world.get_field(field_name).total_amount()
	sim._allocate_membrane_transport(config.tick_dt_min)
	_assert_close(cell.pool("W1"), 0.0, 1e-12, "creating a W1 field does not gift the ancestor secondary-metabolite uptake")
	_assert_close(sim.world.get_field(field_name).total_amount(), extracellular_before, 1e-12, "secondary material is untouched until matching transport machinery exists")

func _test_lysis_conserves_modeled_structural_material() -> void:
	var config = _empty_environment_config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	cell.set_pool("C3", 1.25)
	cell.set_pool("C2", 0.75)
	cell.set_pool("W1", 2.0)
	cell.set_pool("W2", 1.5)
	cell.set_pool("AA", 0.8)
	cell.set_pool("LIP", 0.6)
	cell.set_pool("NUC", 0.4)

	var chemical_before: Dictionary = MetabolicSolverScript.structural_totals(cell.metabolites)
	var expression_before: Dictionary = ExpressionSystemScript.structural_storage_totals(cell.expression_state, config)
	var expected: Dictionary = {
		"C": float(chemical_before["C"]) + float(expression_before["C"]),
		"N": float(chemical_before["N"]) + float(expression_before["N"]),
		"P": float(chemical_before["P"]) + float(expression_before["P"])
	}
	var release: Dictionary = cell.releasable_pools(config)
	var observed: Dictionary = _structural_totals_from_release(release)
	_assert_close(float(observed["C"]), float(expected["C"]), 1e-9, "lysis conserves modeled structural carbon including BIO and protein/RNA material")
	_assert_close(float(observed["N"]), float(expected["N"]), 1e-9, "lysis conserves modeled structural nitrogen including BIO and protein/RNA material")
	_assert_close(float(observed["P"]), float(expected["P"]), 1e-9, "lysis conserves modeled structural phosphorus including BIO and RNA material")

	var bio: float = cell.pool("BIO")
	var aa_field: String = MetaboliteCatalogScript.extracellular_field("AA")
	var lip_field: String = MetaboliteCatalogScript.extracellular_field("LIP")
	var nuc_field: String = MetaboliteCatalogScript.extracellular_field("NUC")
	_assert_true(float(release[aa_field]) >= cell.pool("AA") + 2.0 * bio, "lysis hydrolyzes BIO into reusable amino-acid precursor")
	_assert_close(float(release[lip_field]), cell.pool("LIP") + bio, 1e-12, "lysis hydrolyzes BIO into exact lipid precursor stoichiometry")
	_assert_true(float(release[nuc_field]) >= cell.pool("NUC") + 2.0 * bio, "lysis hydrolyzes BIO and RNA into reusable nucleotide precursor")

func _structural_totals_from_release(release: Dictionary) -> Dictionary:
	var result: Dictionary = {"C": 0.0, "N": 0.0, "P": 0.0}
	var fields: Array = release.keys()
	fields.sort()
	for field_variant in fields:
		var field_name := String(field_variant)
		var metabolite_id: String = MetaboliteCatalogScript.extracellular_metabolite_for_field(field_name)
		assert(not metabolite_id.is_empty(), "Release contains only catalogued extracellular fields")
		var units: Dictionary = MetaboliteCatalogScript.structural_units(metabolite_id)
		var amount: float = maxf(0.0, float(release[field_name]))
		for element in ["C", "N", "P"]:
			result[element] = float(result[element]) + amount * float(units[element])
	return result

func _test_lysis_release_reaches_world_and_event_log() -> void:
	var config = _empty_environment_config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	cell.set_pool("W1", 3.0)
	cell.set_pool("X", 1.25)
	var expected_all: Dictionary = cell.releasable_pools(config)
	var expected_positive: Dictionary = _positive_only(expected_all)
	cell.alive = false
	cell.death_reason = "controlled_lysis"
	sim._process_deaths()

	_assert_true(sim.population_size() == 0, "lysed cell is removed from the living population")
	for field_name_variant in expected_positive.keys():
		var field_name := String(field_name_variant)
		_assert_close(sim.world.get_field(field_name).total_amount(), float(expected_positive[field_name]), 1e-9, "lysis transfers %s from cell state into the chamber" % field_name)
	var death_event: Dictionary = sim.event_log[sim.event_log.size() - 1]
	_assert_true(String(death_event["kind"]) == "death" and String(death_event["reason"]) == "controlled_lysis", "lysis remains an ordinary inspectable death event")
	_assert_true(death_event.has("released_pools"), "death event records released molecular inventory")
	_assert_true(death_event["released_pools"] == expected_positive, "death event release ledger exactly matches material added to extracellular fields")

func _positive_only(values: Dictionary) -> Dictionary:
	var result: Dictionary = {}
	var keys: Array = values.keys()
	keys.sort()
	for key_variant in keys:
		var key := String(key_variant)
		var amount: float = maxf(0.0, float(values[key]))
		if amount > 0.0:
			result[key] = amount
	return result

func _test_same_state_replays_lysis_exactly() -> void:
	var config_a = _empty_environment_config()
	var config_b = _empty_environment_config()
	config_a.seed = 712991
	config_b.seed = 712991
	var first = SimulationEngineScript.new(config_a)
	var second = SimulationEngineScript.new(config_b)
	var first_cell = first.seed_ancestor(Vector2(5.0, 7.0))
	var second_cell = second.seed_ancestor(Vector2(5.0, 7.0))
	for metabolite_id in ["C2", "C3", "W1", "W2", "X"]:
		first_cell.set_pool(metabolite_id, 0.5 + float(metabolite_id.length()) * 0.1)
		second_cell.set_pool(metabolite_id, 0.5 + float(metabolite_id.length()) * 0.1)
	first_cell.alive = false
	second_cell.alive = false
	first_cell.death_reason = "replay_fixture"
	second_cell.death_reason = "replay_fixture"
	first._process_deaths()
	second._process_deaths()
	_assert_true(first.event_log == second.event_log, "same state and seed reproduce exact lysis event history")
	_assert_close(first.checksum(), second.checksum(), 1e-12, "same state and seed reproduce exact post-lysis world checksum")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
