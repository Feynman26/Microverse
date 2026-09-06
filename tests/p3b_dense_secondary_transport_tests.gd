extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_dense_activity_vector_matches_scalar_reference_exactly()
	_test_dense_energy_scaling_matches_dictionary_reference_exactly()
	_test_dense_allocator_matches_active_legacy_trajectory_exactly()
	if failures == 0:
		print("PASS: %d P3-B dense secondary-transport tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d P3-B dense secondary-transport tests failed" % [failures, tests_run])
		quit(1)

func _test_dense_activity_vector_matches_scalar_reference_exactly() -> void:
	var config = _config(true)
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(5.0, 5.0))
	_install_exact_transporter(cell, "W1")
	_install_exact_transporter(cell, "AA")
	_install_exact_transporter(cell, "ROS")
	var dense: PackedFloat64Array = MembraneTransportScript.proteome_activities(
		cell.expression_state,
		config.SECONDARY_EXTRACELLULAR_IDS,
		config
	)
	_assert_true(dense.size() == config.SECONDARY_EXTRACELLULAR_IDS.size(), "dense recognition returns one activity per canonical metabolite")
	for metabolite_index in range(config.SECONDARY_EXTRACELLULAR_IDS.size()):
		var metabolite_id: String = config.SECONDARY_EXTRACELLULAR_IDS[metabolite_index]
		var scalar: float = MembraneTransportScript.proteome_activity(cell.expression_state, metabolite_id, config)
		_assert_close(dense[metabolite_index], scalar, 0.0, "dense %s activity preserves exact scalar accumulation" % metabolite_id)

func _test_dense_energy_scaling_matches_dictionary_reference_exactly() -> void:
	var config = _config(true)
	var dense := PackedFloat64Array([0.4, -0.2, 0.0, 0.1])
	var legacy: Dictionary = {"a": 0.4, "b": -0.2, "c": 0.0, "d": 0.1}
	for available_atp in [0.0, 0.001, 0.1, 10.0]:
		_assert_close(
			MembraneTransportScript.energy_scale_dense(dense, available_atp, config),
			MembraneTransportScript.energy_scale(legacy, available_atp, config),
			0.0,
			"dense ATP scaling preserves exact canonical accumulation at ATP=%s" % available_atp
		)

func _test_dense_allocator_matches_active_legacy_trajectory_exactly() -> void:
	var dense = _active_simulation(true)
	var legacy = _active_simulation(false)
	_assert_close(dense.checksum(), legacy.checksum(), 0.0, "paired active transport scenarios begin identically")
	for tick in range(40):
		dense.step(1)
		legacy.step(1)
		if tick == 0:
			_assert_true(
				float(dense.last_secondary_transport_summary["total_moved"]) > 0.0,
				"paired trajectory exercises active secondary exchange"
			)
		_assert_true(
			dense.last_secondary_transport_summary == legacy.last_secondary_transport_summary,
			"dense transport ledger matches frozen M7 allocator at tick %d" % (tick + 1)
		)
		_assert_close(dense.checksum(), legacy.checksum(), 0.0, "dense full state matches frozen M7 allocator at tick %d" % (tick + 1))
		_assert_true(dense.event_log == legacy.event_log, "dense event history matches frozen M7 allocator at tick %d" % (tick + 1))

func _active_simulation(use_dense_allocator: bool):
	var config = _config(use_dense_allocator)
	var sim = SimulationEngineScript.new(config)
	var positions: Array[Vector2] = [
		Vector2(5.05, 5.05),
		Vector2(5.15, 5.15),
		Vector2(7.0, 6.0),
		Vector2(8.0, 8.0)
	]
	for index in range(positions.size()):
		var cell = sim.seed_ancestor(positions[index])
		_install_exact_transporter(cell, "W1")
		_install_exact_transporter(cell, "W2")
		_install_exact_transporter(cell, "ROS")
		cell.set_pool("W1", 0.15 * float(index))
		cell.set_pool("W2", 1.0 + 0.2 * float(index))
		cell.set_pool("ROS", 0.05 * float(index))
	for metabolite_id in ["W1", "W2", "ROS"]:
		var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
		sim.world.release(field_name, Vector2(5.0, 5.0), 2.0)
		sim.world.release(field_name, Vector2(8.0, 8.0), 0.4)
	return sim

func _config(use_dense_allocator: bool):
	var config = SimConfigScript.new()
	config.seed = 930302
	config.world_width = 12
	config.world_height = 12
	config.max_cells = 32
	config.mutation_enabled = false
	config.metabolic_use_dense_solver = true
	config.secondary_transport_use_dense_allocator = use_dense_allocator
	config.validate()
	return config

func _install_exact_transporter(cell, metabolite_id: String) -> void:
	var target: int = MembraneTransportScript.target_signature(metabolite_id)
	var best_locus: int = -1
	var best_distance: int = 99
	for gene in cell.genome.genes:
		var distance: int = MembraneTransportScript.hamming_distance(int(gene.protein_signature), target)
		if distance < best_distance:
			best_distance = distance
			best_locus = int(gene.locus_id)
	var selected_gene = cell.genome.get_gene_by_locus(best_locus)
	selected_gene.protein_signature = target
	var locus_state: Dictionary = cell.expression_state[best_locus]
	var abundance: float = 0.0
	for cohort_amount in locus_state["protein"].values():
		abundance += float(cohort_amount)
	locus_state["protein"] = {target: abundance}

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
