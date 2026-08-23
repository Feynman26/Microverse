extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const CellMechanicsScript = preload("res://src/physics/cell_mechanics.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_radius_scales_with_area()
	_test_wall_projection_contains_complete_disk()
	_test_pair_relaxation_removes_material_overlap()
	_test_symmetric_pair_preserves_midpoint()
	_test_insertion_order_does_not_change_equilibrium()
	_test_division_remains_local_and_relaxes_overlap()

	if failures == 0:
		print("PASS: %d M6 physical-cell mechanics tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M6 physical-cell mechanics tests failed" % [failures, tests_run])
		quit(1)

func _config():
	var config = SimConfigScript.new()
	config.mechanical_relaxation_iterations = 32
	config.mechanical_overlap_tolerance = 1e-5
	return config

func _test_radius_scales_with_area() -> void:
	var config = _config()
	var r1: float = CellMechanicsScript.radius_for_volume(config.ancestor_volume, config)
	var r4: float = CellMechanicsScript.radius_for_volume(config.ancestor_volume * 4.0, config)
	_assert_close(r4, 2.0 * r1, 1e-12, "2D disk radius follows sqrt(volume/area) scaling")

func _test_wall_projection_contains_complete_disk() -> void:
	var config = _config()
	var world = WorldStateScript.new(16, 16, 1.0)
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2.ZERO, config.ancestor_volume)
	CellMechanicsScript.relax([cell], world, config)
	var radius: float = CellMechanicsScript.radius_for_cell(cell, config)
	_assert_close(cell.position.x, radius, 1e-12, "left wall constrains disk edge rather than only its center")
	_assert_close(cell.position.y, radius, 1e-12, "top wall constrains disk edge rather than only its center")
	_assert_true(CellMechanicsScript.disks_within_bounds([cell], world, config), "wall projection keeps the complete cell disk inside the chamber")

func _test_pair_relaxation_removes_material_overlap() -> void:
	var config = _config()
	var world = WorldStateScript.new(16, 16, 1.0)
	var first = CellStateScript.new(1, -1, 0, 0, Vector2(7.5, 8.0), 1.0)
	var second = CellStateScript.new(2, -1, 0, 0, Vector2(7.8, 8.0), 1.0)
	var before: float = CellMechanicsScript.max_overlap([first, second], config)
	var summary: Dictionary = CellMechanicsScript.relax([first, second], world, config)
	_assert_true(before > 0.1, "controlled pair begins with material overlap")
	_assert_true(float(summary["max_overlap"]) <= config.mechanical_overlap_tolerance * 1.1, "mechanical relaxation removes pair overlap to configured tolerance")
	_assert_true(CellMechanicsScript.disks_within_bounds([first, second], world, config), "relaxed cells remain inside chamber walls")

func _test_symmetric_pair_preserves_midpoint() -> void:
	var config = _config()
	var world = WorldStateScript.new(16, 16, 1.0)
	var first = CellStateScript.new(1, -1, 0, 0, Vector2(7.6, 8.0), 1.0)
	var second = CellStateScript.new(2, -1, 0, 0, Vector2(8.0, 8.0), 1.0)
	var midpoint_before: Vector2 = 0.5 * (first.position + second.position)
	CellMechanicsScript.relax([first, second], world, config)
	var midpoint_after: Vector2 = 0.5 * (first.position + second.position)
	_assert_close(midpoint_after.x, midpoint_before.x, 1e-12, "equal/opposite contact correction preserves symmetric x midpoint")
	_assert_close(midpoint_after.y, midpoint_before.y, 1e-12, "equal/opposite contact correction preserves symmetric y midpoint")

func _test_insertion_order_does_not_change_equilibrium() -> void:
	var config = _config()
	var world = WorldStateScript.new(20, 20, 1.0)
	var forward: Array = [
		CellStateScript.new(1, -1, 0, 0, Vector2(8.0, 10.0), 1.0),
		CellStateScript.new(2, -1, 0, 0, Vector2(8.45, 10.0), 1.0),
		CellStateScript.new(3, -1, 0, 0, Vector2(8.90, 10.1), 1.0)
	]
	var reverse: Array = [
		CellStateScript.new(3, -1, 0, 0, Vector2(8.90, 10.1), 1.0),
		CellStateScript.new(2, -1, 0, 0, Vector2(8.45, 10.0), 1.0),
		CellStateScript.new(1, -1, 0, 0, Vector2(8.0, 10.0), 1.0)
	]
	CellMechanicsScript.relax(forward, world, config)
	CellMechanicsScript.relax(reverse, world, config)
	var forward_by_id: Dictionary = _positions_by_id(forward)
	var reverse_by_id: Dictionary = _positions_by_id(reverse)
	for id_variant in forward_by_id.keys():
		var cell_id: int = int(id_variant)
		var a: Vector2 = forward_by_id[cell_id]
		var b: Vector2 = reverse_by_id[cell_id]
		_assert_close(a.x, b.x, 1e-12, "reversing insertion order preserves x equilibrium for cell %d" % cell_id)
		_assert_close(a.y, b.y, 1e-12, "reversing insertion order preserves y equilibrium for cell %d" % cell_id)

func _test_division_remains_local_and_relaxes_overlap() -> void:
	var config = _config()
	var world = WorldStateScript.new(20, 20, 1.0)
	var rng = DeterministicRngScript.new(84211)
	var parent_position: Vector2 = Vector2(10.0, 10.0)
	var parent = CellStateScript.new(10, -1, 0, 0, parent_position, config.division_volume)
	parent.genome = GenomeScript.create_ancestor()
	parent.initialize_molecular_state(config)
	parent.set_pool("BIO", config.division_volume * config.biomass_units_per_volume)
	parent.volume = config.division_volume
	parent.set_pool("ATP", 5.0)
	var daughters: Array = parent.create_daughters(11, 12, 1, rng, world, config)
	CellMechanicsScript.relax(daughters, world, config)
	_assert_true(CellMechanicsScript.max_overlap(daughters, config) <= config.mechanical_overlap_tolerance * 1.1, "local daughters do not retain material overlap after mechanics")
	_assert_true(daughters[0].position.distance_to(parent_position) < 1.5 and daughters[1].position.distance_to(parent_position) < 1.5, "division daughters remain local to the parent neighborhood")
	var total_volume: float = daughters[0].volume + daughters[1].volume
	var center_of_volume: Vector2 = (daughters[0].position * daughters[0].volume + daughters[1].position * daughters[1].volume) / total_volume
	_assert_true(center_of_volume.distance_to(parent_position) < 0.05, "local division approximately conserves center of structural volume")

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
