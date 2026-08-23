extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const CellMechanicsScript = preload("res://src/physics/cell_mechanics.gd")
const SpatialHashScript = preload("res://src/physics/spatial_hash.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_all_true_contacts_are_candidates()
	_test_candidate_order_independent_of_insertion()
	_test_spatial_solver_matches_reference_exactly()
	_test_sparse_population_prunes_pair_search()
	if failures == 0:
		print("PASS: %d M6 spatial-index tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M6 spatial-index tests failed" % [failures, tests_run])
		quit(1)

func _config():
	var config = SimConfigScript.new()
	config.mechanical_relaxation_iterations = 40
	config.mechanical_overlap_tolerance = 1e-5
	config.mechanical_neighbor_bucket_size_grid = 2.0
	return config

func _fixture_cells() -> Array:
	return [
		CellStateScript.new(1, -1, 0, 0, Vector2(8.0, 8.0), 1.0),
		CellStateScript.new(2, -1, 0, 0, Vector2(8.35, 8.0), 1.15),
		CellStateScript.new(3, -1, 0, 0, Vector2(8.75, 8.20), 0.9),
		CellStateScript.new(4, -1, 0, 0, Vector2(9.10, 8.55), 1.2),
		CellStateScript.new(5, -1, 0, 0, Vector2(13.0, 12.0), 1.0),
		CellStateScript.new(6, -1, 0, 0, Vector2(13.4, 12.0), 1.0)
	]

func _test_all_true_contacts_are_candidates() -> void:
	var config = _config()
	var cells: Array = _fixture_cells()
	var candidate_keys: Dictionary = {}
	for key_variant in SpatialHashScript.pair_keys(cells, config, config.mechanical_neighbor_bucket_size_grid):
		candidate_keys[key_variant] = true
	var true_contacts: int = 0
	for i in range(cells.size()):
		for j in range(i + 1, cells.size()):
			var target: float = CellMechanicsScript.radius_for_cell(cells[i], config) + CellMechanicsScript.radius_for_cell(cells[j], config)
			if cells[i].position.distance_to(cells[j].position) < target:
				true_contacts += 1
				var key := Vector2i(int(cells[i].id), int(cells[j].id))
				_assert_true(candidate_keys.has(key), "spatial hash retains true contact %s" % key)
	_assert_true(true_contacts > 0, "contact-coverage fixture contains real overlaps")

func _test_candidate_order_independent_of_insertion() -> void:
	var config = _config()
	var forward: Array = _fixture_cells()
	var reverse: Array = _fixture_cells()
	reverse.reverse()
	var first_keys: Array = SpatialHashScript.pair_keys(forward, config, config.mechanical_neighbor_bucket_size_grid)
	var second_keys: Array = SpatialHashScript.pair_keys(reverse, config, config.mechanical_neighbor_bucket_size_grid)
	_assert_true(first_keys == second_keys, "spatial candidate pairs are canonical under insertion-order reversal")

func _test_spatial_solver_matches_reference_exactly() -> void:
	var config = _config()
	var world = WorldStateScript.new(24, 24, 1.0)
	var reference: Array = _fixture_cells()
	var indexed: Array = _fixture_cells()
	var reference_summary: Dictionary = CellMechanicsScript.relax(reference, world, config, false)
	var indexed_summary: Dictionary = CellMechanicsScript.relax(indexed, world, config, true)
	var reference_positions: Dictionary = _positions_by_id(reference)
	var indexed_positions: Dictionary = _positions_by_id(indexed)
	for id_variant in reference_positions.keys():
		var cell_id: int = int(id_variant)
		var a: Vector2 = reference_positions[cell_id]
		var b: Vector2 = indexed_positions[cell_id]
		_assert_close(a.x, b.x, 1e-12, "spatial broad phase preserves reference x for cell %d" % cell_id)
		_assert_close(a.y, b.y, 1e-12, "spatial broad phase preserves reference y for cell %d" % cell_id)
	_assert_close(float(indexed_summary["max_overlap"]), float(reference_summary["max_overlap"]), 1e-12, "indexed and reference solvers end with identical overlap")

func _test_sparse_population_prunes_pair_search() -> void:
	var config = _config()
	var cells: Array = []
	for i in range(12):
		cells.append(CellStateScript.new(i + 1, -1, 0, 0, Vector2(2.0 + float(i % 4) * 12.0, 2.0 + float(i / 4) * 12.0), 1.0))
	var indexed_count: int = SpatialHashScript.candidate_pairs(cells, config, config.mechanical_neighbor_bucket_size_grid).size()
	var all_pair_count: int = cells.size() * (cells.size() - 1) / 2
	_assert_true(indexed_count < all_pair_count, "spatial hash prunes distant pairs from the broad phase")

func _positions_by_id(cells: Array) -> Dictionary:
	var result: Dictionary = {}
	for cell in cells:
		result[int(cell.id)] = cell.position
	return result

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
