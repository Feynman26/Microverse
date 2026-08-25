extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const ChemicalFieldScript = preload("res://src/world/chemical_field.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_basal_reference_flux_is_molecularly_reproduced()
	_test_basal_transport_changes_only_when_realized_protein_changes()
	_test_stable_regime_keeps_exact_explicit_solver()
	_test_high_alpha_implicit_diffusion_is_mass_conserving_and_nonnegative()
	_test_high_alpha_diffusion_preserves_uniform_field()
	if failures == 0:
		print("PASS: %d M11 membrane/diffusion prerequisite tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M11 membrane/diffusion prerequisite tests failed" % [failures, tests_run])
		quit(1)

func _reference_cell_and_world() -> Dictionary:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	var world = WorldStateScript.new(8, 8, config.grid_cell_size_um)
	world.register_field("glucose", config.glucose_diffusion, config.initial_glucose)
	world.register_field("oxygen", config.oxygen_diffusion, config.initial_oxygen)
	world.register_field("nitrogen", config.nitrogen_diffusion, config.initial_nitrogen)
	world.register_field("phosphorus", config.phosphorus_diffusion, config.initial_phosphorus)
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2(4, 4), config.ancestor_volume)
	cell.genome = GenomeScript.create_ancestor()
	cell.initialize_molecular_state(config)
	return {"config": config, "world": world, "cell": cell}

func _test_basal_reference_flux_is_molecularly_reproduced() -> void:
	var fixture: Dictionary = _reference_cell_and_world()
	var config = fixture["config"]
	var world = fixture["world"]
	var cell = fixture["cell"]
	var dt: float = 0.1
	var cases: Array = [
		["G", "glucose", config.glucose_transport_vmax, config.glucose_transport_km],
		["O2", "oxygen", config.oxygen_transport_vmax, config.oxygen_transport_km],
		["NH4", "nitrogen", config.nitrogen_transport_vmax, config.nitrogen_transport_km],
		["P", "phosphorus", config.phosphorus_transport_vmax, config.phosphorus_transport_km]
	]
	for case in cases:
		var metabolite_id: String = String(case[0])
		var field_name: String = String(case[1])
		var vmax: float = float(case[2])
		var km: float = float(case[3])
		var external: float = float(world.sample(field_name, cell.position))
		var activity: float = MembraneTransportScript.basal_proteome_activity(cell.expression_state, metabolite_id, config)
		var molecular: float = MembraneTransportScript.basal_import_request(
			cell.pool(metabolite_id), cell.volume, external, activity, metabolite_id, vmax, km, dt, config
		)
		var legacy: float = minf(
			vmax * external / (km + external) * dt,
			float(config.intracellular_pool_capacity_per_volume) * cell.volume - cell.pool(metabolite_id)
		)
		_assert_close(molecular, legacy, 1e-12, "ancestral %s uptake reproduces historical reference through realized proteome" % metabolite_id)

func _test_basal_transport_changes_only_when_realized_protein_changes() -> void:
	var fixture: Dictionary = _reference_cell_and_world()
	var config = fixture["config"]
	var world = fixture["world"]
	var cell = fixture["cell"]
	var external: float = float(world.sample("glucose", cell.position))
	var before_activity: float = MembraneTransportScript.basal_proteome_activity(cell.expression_state, "G", config)
	var before: float = MembraneTransportScript.basal_import_request(
		cell.pool("G"), cell.volume, external, before_activity, "G",
		config.glucose_transport_vmax, config.glucose_transport_km, 0.1, config
	)

	# DNA changes first; the historical cohort remains physical and therefore
	# membrane transport must not change instantaneously.
	cell.genome.get_gene_by_locus(1).protein_signature ^= 1
	var dna_only_activity: float = MembraneTransportScript.basal_proteome_activity(cell.expression_state, "G", config)
	_assert_close(dna_only_activity, before_activity, 1e-15, "DNA mutation alone does not rewrite existing basal transporter proteins")

	# Replace only the realized locus-1 cohort with the one-bit mutant to model
	# post-turnover expression. Affinity, and therefore uptake, must now change.
	var cohorts: Dictionary = cell.expression_state[1]["protein"]
	var ancestral_amount: float = float(cohorts.get(0x1357, 0.0))
	cohorts.erase(0x1357)
	cohorts[0x1356] = ancestral_amount
	var after_activity: float = MembraneTransportScript.basal_proteome_activity(cell.expression_state, "G", config)
	var after: float = MembraneTransportScript.basal_import_request(
		cell.pool("G"), cell.volume, external, after_activity, "G",
		config.glucose_transport_vmax, config.glucose_transport_km, 0.1, config
	)
	_assert_true(after_activity < before_activity, "one-bit realized membrane mutation lowers target affinity")
	_assert_true(after < before, "changed realized membrane affinity changes basal uptake without a phenotype flag")

func _test_stable_regime_keeps_exact_explicit_solver() -> void:
	var field = ChemicalFieldScript.new(5, 5, 1.0, 1.0, 0.0)
	field.set_value(2, 2, 10.0)
	var expected := PackedFloat64Array()
	expected.resize(25)
	expected.fill(0.0)
	# alpha=0.1. Historical 5-point explicit update from a single center pulse.
	expected[2 + 2 * 5] = 6.0
	expected[1 + 2 * 5] = 1.0
	expected[3 + 2 * 5] = 1.0
	expected[2 + 1 * 5] = 1.0
	expected[2 + 3 * 5] = 1.0
	field.step_diffusion(0.1)
	_assert_true(field.last_diffusion_mode == "explicit", "stable historical alpha remains on exact explicit solver")
	for i in range(expected.size()):
		_assert_close(field.values[i], expected[i], 1e-15, "stable explicit lattice value %d remains exact" % i)

func _test_high_alpha_implicit_diffusion_is_mass_conserving_and_nonnegative() -> void:
	var field = ChemicalFieldScript.new(9, 9, 1.0, 40000.0, 0.0)
	field.set_value(4, 4, 100.0)
	var before: float = field.total_amount()
	field.step_diffusion(0.1) # alpha=4000: impossible for the old explicit step.
	var after: float = field.total_amount()
	_assert_true(field.last_diffusion_mode == "implicit_split", "high-alpha diffusion switches to stable implicit split")
	_assert_close(after, before, 1e-8, "high-alpha implicit diffusion conserves closed-field mass")
	_assert_true(field.minimum_value() >= -1e-12, "high-alpha implicit diffusion remains nonnegative")
	_assert_true(field.get_value(4, 4) < 100.0 and field.get_value(0, 0) > 0.0, "high-alpha solver spreads a localized pulse across the chamber")

func _test_high_alpha_diffusion_preserves_uniform_field() -> void:
	var field = ChemicalFieldScript.new(7, 6, 1.0, 100000.0, 3.25)
	var before: float = field.total_amount()
	field.step_diffusion(0.1)
	_assert_true(field.last_diffusion_mode == "implicit_split", "uniform high-alpha field uses implicit numerical path")
	_assert_close(field.total_amount(), before, 1e-8, "reflecting implicit diffusion preserves uniform-field mass")
	for value in field.values:
		_assert_close(float(value), 3.25, 1e-10, "uniform field is fixed point of reflecting diffusion")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
