extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_capacity_and_equilibrium_scale_with_biomass()
	_test_currency_capacity_scales_with_biomass_without_free_energy()
	_test_sensing_is_concentration_based()
	_test_division_preserves_machinery_concentration()
	_test_lysis_recycles_modeled_nutrients()
	if failures == 0:
		print("PASS: %d M10 cell-scale tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 cell-scale tests failed" % [failures, tests_run])
		quit(1)

func _test_capacity_and_equilibrium_scale_with_biomass() -> void:
	var config = SimConfigScript.new()
	var unit = _new_cell(1, 1.0, config)
	var double = _new_cell(2, 2.0, config)

	_assert_close(double.proteome_capacity(config), 2.0 * unit.proteome_capacity(config), 1e-10, "doubling biomass doubles physical proteome capacity")
	_assert_close(double.total_protein(), 2.0 * unit.total_protein(), 1e-10, "equilibrium protein inventory scales with biomass instead of remaining fixed per cell")
	_assert_close(double.total_mrna(), 2.0 * unit.total_mrna(), 1e-10, "equilibrium mRNA inventory scales with biomass")
	_assert_close(double.total_protein() / double.volume, unit.total_protein() / unit.volume, 1e-10, "equal genotype at equal state keeps protein concentration invariant across cell size")

	var genome = GenomeScript.create_ancestor()
	var state_unit: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 1.0)
	var state_double: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 2.0)
	var pools_unit: Dictionary = _abundant_expression_pools(1.0, config)
	var pools_double: Dictionary = _abundant_expression_pools(2.0, config)
	var summary_unit: Dictionary = ExpressionSystemScript.step(state_unit, genome, pools_unit, config.tick_dt_min, DeterministicRngScript.new(701), config, 1.0)
	var summary_double: Dictionary = ExpressionSystemScript.step(state_double, genome, pools_double, config.tick_dt_min, DeterministicRngScript.new(701), config, 2.0)
	_assert_close(float(summary_double["translation_capacity"]), 2.0 * float(summary_unit["translation_capacity"]), 1e-12, "shared translation throughput scales linearly with physical cell size")

func _test_currency_capacity_scales_with_biomass_without_free_energy() -> void:
	var config = SimConfigScript.new()
	var cell = _new_cell(11, 1.0, config)
	var initial_atp: float = cell.pool("ATP")
	var initial_nadh: float = cell.pool("NADH")
	var adenylate_density: float = cell.total_adenylate() / cell.volume
	var redox_density: float = cell.total_redox_currency() / cell.volume

	# Controlled biomass growth isolates carrier-capacity scaling from metabolic
	# flux. Carrier scaffolds are implicit material in this model; only their
	# charge/redox state is explicit. Therefore growth may add ADP/NAD capacity,
	# but must not manufacture ATP or NADH energy equivalents.
	cell.set_pool("BIO", 2.0 * float(config.biomass_units_per_volume))
	cell._sync_volume_from_biomass(config)
	_assert_close(cell.volume, 2.0, 1e-12, "controlled structural growth doubles cell volume")
	_assert_close(cell.total_adenylate() / cell.volume, adenylate_density, 1e-12, "adenylate carrier capacity scales with biomass instead of being founder-limited")
	_assert_close(cell.total_redox_currency() / cell.volume, redox_density, 1e-12, "redox carrier capacity scales with biomass instead of being founder-limited")
	_assert_close(cell.pool("ATP"), initial_atp, 1e-12, "biomass growth does not gift high-energy ATP")
	_assert_close(cell.pool("NADH"), initial_nadh, 1e-12, "biomass growth does not gift reducing power")
	_assert_true(cell.pool("ADP") > float(config.initial_adp_per_volume), "new adenylate capacity appears only in discharged ADP form")
	_assert_true(cell.pool("NAD") > float(config.initial_nad_per_volume), "new redox capacity appears only in oxidized NAD form")

func _test_sensing_is_concentration_based() -> void:
	var config = SimConfigScript.new()
	var genome = GenomeScript.create_ancestor()
	var state_unit: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 1.0)
	var state_double: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 2.0)
	var pools_unit: Dictionary = _abundant_expression_pools(1.0, config)
	var pools_double: Dictionary = _abundant_expression_pools(2.0, config)
	# Give both cells the same ligand concentrations by doubling every molecular
	# amount in the larger cell. Regulation must therefore see the same chemistry.
	var summary_unit: Dictionary = ExpressionSystemScript.step(state_unit, genome, pools_unit, config.tick_dt_min, DeterministicRngScript.new(991), config, 1.0)
	var summary_double: Dictionary = ExpressionSystemScript.step(state_double, genome, pools_double, config.tick_dt_min, DeterministicRngScript.new(991), config, 2.0)
	for locus_variant in summary_unit["regulation"].keys():
		var locus_id: int = int(locus_variant)
		_assert_close(
			float(summary_double["regulation"][locus_id]),
			float(summary_unit["regulation"][locus_id]),
			1e-12,
			"regulatory occupancy for locus %d depends on concentration rather than absolute cell size" % locus_id
		)

func _test_division_preserves_machinery_concentration() -> void:
	var config = SimConfigScript.new()
	config.partition_jitter = 0.0
	config.expression_partition_noise_scale = 0.0
	var mother = _new_cell(1, config.division_volume, config)
	var world = WorldStateScript.new(config.world_width, config.world_height, config.grid_cell_size_um)
	var rng = DeterministicRngScript.new(424242)
	var protein_before: float = mother.total_protein()
	var mrna_before: float = mother.total_mrna()
	var adenylate_before: float = mother.total_adenylate()
	var redox_before: float = mother.total_redox_currency()
	var protein_concentration_before: float = protein_before / mother.volume
	var mrna_concentration_before: float = mrna_before / mother.volume
	var adenylate_concentration_before: float = adenylate_before / mother.volume
	var redox_concentration_before: float = redox_before / mother.volume

	_assert_true(mother.ready_to_divide(config), "controlled biomass-scaled mother is ordinarily ready to divide")
	var daughters: Array = mother.create_daughters(2, 3, 10, rng, world, config)
	_assert_close(float(daughters[0].volume), 1.0, 1e-10, "symmetric division returns first daughter to reference cell volume")
	_assert_close(float(daughters[1].volume), 1.0, 1e-10, "symmetric division returns second daughter to reference cell volume")
	_assert_close(float(daughters[0].total_protein() + daughters[1].total_protein()), protein_before, 1e-10, "division conserves explicit protein inventory")
	_assert_close(float(daughters[0].total_mrna() + daughters[1].total_mrna()), mrna_before, 1e-10, "division conserves explicit mRNA inventory")
	_assert_close(float(daughters[0].total_adenylate() + daughters[1].total_adenylate()), adenylate_before, 1e-10, "division conserves adenylate carrier inventory")
	_assert_close(float(daughters[0].total_redox_currency() + daughters[1].total_redox_currency()), redox_before, 1e-10, "division conserves redox carrier inventory")
	for daughter in daughters:
		_assert_close(daughter.total_protein() / daughter.volume, protein_concentration_before, 1e-10, "daughter inherits maternal protein concentration without artificial post-division dilution")
		_assert_close(daughter.total_mrna() / daughter.volume, mrna_concentration_before, 1e-10, "daughter inherits maternal mRNA concentration without artificial post-division dilution")
		_assert_close(daughter.total_adenylate() / daughter.volume, adenylate_concentration_before, 1e-10, "daughter inherits biomass-scaled adenylate capacity")
		_assert_close(daughter.total_redox_currency() / daughter.volume, redox_concentration_before, 1e-10, "daughter inherits biomass-scaled redox capacity")
		_assert_true(daughter.total_protein() <= daughter.proteome_capacity(config) + 1e-8, "daughter inherited proteome remains inside its volume-scaled physical capacity")

func _test_lysis_recycles_modeled_nutrients() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(4.0, 4.0))
	cell.set_pool("G", 1.25)
	cell.set_pool("O2", 0.75)
	cell.set_pool("NH4", 2.50)
	cell.set_pool("P", 1.50)
	cell.set_pool("C2", 3.00)
	cell.set_pool("BIO", 2.0)
	cell._sync_volume_from_biomass(config)
	var expected: Dictionary = cell.releasable_pools(config)
	var before: Dictionary = {}
	for field_name in sim.world.field_order:
		before[field_name] = sim.world.get_field(field_name).total_amount()

	cell.alive = false
	cell.death_reason = "controlled_lysis"
	sim.step(1)
	_assert_true(sim.population_size() == 0, "controlled dead cell is removed from living population")
	for field_variant in expected.keys():
		var field_name: String = String(field_variant)
		var expected_delta: float = float(expected[field_name])
		var actual_delta: float = sim.world.get_field(field_name).total_amount() - float(before.get(field_name, 0.0))
		_assert_close(actual_delta, expected_delta, 1e-8, "lysis returns modeled nutrient pool %s to extracellular system" % field_name)
	var death_event: Dictionary = {}
	for event_variant in sim.event_log:
		if String(event_variant.get("kind", "")) == "death":
			death_event = event_variant
	_assert_true(not death_event.is_empty(), "lysis remains an explicit inspectable death event")
	_assert_true(String(death_event.get("reason", "")) == "controlled_lysis", "death ledger preserves physiological/experimental cause")
	_assert_true(death_event.has("released_pools"), "death ledger records exactly which nutrients returned to the chamber")

func _new_cell(cell_id: int, volume: float, config):
	var cell = CellStateScript.new(cell_id, -1, 0, 0, Vector2(4.0, 4.0), volume)
	cell.genome = GenomeScript.create_ancestor()
	cell.initialize_molecular_state(config)
	return cell

func _abundant_expression_pools(volume: float, config) -> Dictionary:
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(volume, config)
	# Scale all relevant molecular inventories with volume so the two controlled
	# cells begin at identical concentrations. The large values prevent resource
	# scarcity from obscuring the machinery-capacity test.
	pools["ATP"] = 1000.0 * volume
	pools["ADP"] = 1000.0 * volume
	pools["AA"] = 1000.0 * volume
	pools["NUC"] = 1000.0 * volume
	pools["O2"] = 2.0 * volume
	pools["X"] = 1.0 * volume
	return pools

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])