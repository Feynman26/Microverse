extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_canonical_catalog_indices_preserve_historical_order()
	_test_dense_solver_matches_frozen_m10_substeps_exactly()
	_test_compiled_network_preserves_reaction_order_variants()
	_test_cell_workspace_is_reused_without_becoming_authoritative()
	_test_complete_simulation_matches_legacy_solver_every_tick()
	if failures == 0:
		print("PASS: %d P3 dense-metabolism tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d P3 dense-metabolism tests failed" % [failures, tests_run])
		quit(1)

func _test_canonical_catalog_indices_preserve_historical_order() -> void:
	var historical: Array[String] = []
	for key in MetaboliteCatalogScript.DEFINITIONS.keys():
		historical.append(String(key))
	historical.sort()
	_assert_true(MetaboliteCatalogScript.ids() == historical, "dense metabolite order exactly matches historical sorted IDs")
	for index in range(historical.size()):
		_assert_true(MetaboliteCatalogScript.index_of(historical[index]) == index, "metabolite index round-trips %s" % historical[index])
		_assert_true(MetaboliteCatalogScript.id_at(index) == historical[index], "metabolite ID round-trips index %d" % index)

func _test_dense_solver_matches_frozen_m10_substeps_exactly() -> void:
	var config = _config(true)
	var cell = _cell(config)
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var dense_pools: Dictionary = cell.metabolites.duplicate(true)
	var legacy_pools: Dictionary = cell.metabolites.duplicate(true)
	var workspace: Dictionary = {}
	var compiled = MetabolicSolverScript.compile_network(reactions)
	for _iteration in range(8):
		var dense_flux: Dictionary = MetabolicSolverScript.step(
			dense_pools, cell.genome, cell.expression_state, reactions, 0.37, cell.volume,
			config, compiled, workspace, true
		)
		var legacy_flux: Dictionary = MetabolicSolverScript.step_legacy_reference(
			legacy_pools, cell.genome, cell.expression_state, reactions, 0.37, cell.volume, config
		)
		_assert_true(dense_flux == legacy_flux, "dense fluxes exactly match frozen M10 across repeated substeps")
		_assert_true(dense_pools == legacy_pools, "dense pools exactly match frozen M10 across repeated substeps")

func _test_compiled_network_preserves_reaction_order_variants() -> void:
	var config = _config(true)
	var cell = _cell(config)
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	reactions.reverse()
	var variants: Array = [reactions, [reactions[0], reactions[4], reactions[8]]]
	var labels: Array[String] = ["reversed", "partial"]
	for variant_index in range(variants.size()):
		var variant: Array = variants[variant_index]
		var dense_pools: Dictionary = cell.metabolites.duplicate(true)
		var legacy_pools: Dictionary = cell.metabolites.duplicate(true)
		var dense_flux: Dictionary = MetabolicSolverScript.step(
			dense_pools, cell.genome, cell.expression_state, variant, 0.5, cell.volume,
			config, MetabolicSolverScript.compile_network(variant), {}, true
		)
		var legacy_flux: Dictionary = MetabolicSolverScript.step_legacy_reference(
			legacy_pools, cell.genome, cell.expression_state, variant, 0.5, cell.volume, config
		)
		_assert_true(dense_flux == legacy_flux, "compiled %s network retains exact legacy accumulation order" % labels[variant_index])
		_assert_true(dense_pools == legacy_pools, "compiled %s network retains exact legacy pool result" % labels[variant_index])

func _test_cell_workspace_is_reused_without_becoming_authoritative() -> void:
	var config = _config(true)
	var cell = _cell(config)
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var compiled = MetabolicSolverScript.compile_network(reactions)
	MetabolicSolverScript.step(
		cell.metabolites, cell.genome, cell.expression_state, reactions, 0.1, cell.volume,
		config, compiled, cell.metabolic_workspace, true
	)
	var keys_after_first: Array = cell.metabolic_workspace.keys()
	keys_after_first.sort()
	var checksum_before_workspace_reuse: float = cell.checksum()
	MetabolicSolverScript.step(
		cell.metabolites, cell.genome, cell.expression_state, reactions, 0.0, cell.volume,
		config, compiled, cell.metabolic_workspace, true
	)
	var keys_after_second: Array = cell.metabolic_workspace.keys()
	keys_after_second.sort()
	_assert_true(keys_after_first == keys_after_second and keys_after_first.size() == 8, "cell reuses the complete dense buffer set")
	_assert_close(cell.checksum(), checksum_before_workspace_reuse, 0.0, "derived workspace reuse does not enter authoritative checksum")

func _test_complete_simulation_matches_legacy_solver_every_tick() -> void:
	var dense_config = _config(true)
	var legacy_config = _config(false)
	var dense = SimulationEngineScript.new(dense_config)
	var legacy = SimulationEngineScript.new(legacy_config)
	dense.seed_ancestor()
	legacy.seed_ancestor()
	for tick in range(80):
		dense.step(1)
		legacy.step(1)
		_assert_close(dense.checksum(), legacy.checksum(), 0.0, "full dense trajectory matches legacy at tick %d" % (tick + 1))
		_assert_true(dense.event_log == legacy.event_log, "full dense event history matches legacy at tick %d" % (tick + 1))

func _cell(config):
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2(5.0, 5.0), 1.0)
	cell.genome = GenomeScript.create_ancestor()
	cell.initialize_molecular_state(config)
	cell.set_pool("G", 3.0)
	cell.set_pool("O2", 2.0)
	cell.set_pool("NH4", 2.0)
	cell.set_pool("P", 2.0)
	cell.set_pool("C3", 1.0)
	cell.set_pool("W1", 0.5)
	cell.set_pool("W2", 0.5)
	return cell

func _config(dense: bool):
	var config = SimConfigScript.new()
	config.seed = 930301
	config.world_width = 12
	config.world_height = 12
	config.max_cells = 16
	config.mutation_enabled = false
	config.metabolic_use_dense_solver = dense
	config.validate()
	return config

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
