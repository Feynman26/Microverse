extends SceneTree

const WorldStateScript = preload("res://src/world/world_state.gd")
const ChemicalFieldScript = preload("res://src/world/chemical_field.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_exact_zero_fast_path()
	_test_parallel_diffusion_matches_serial_exactly()
	if failures == 0:
		print("PASS: %d M10 performance-equivalence tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 performance-equivalence tests failed" % [failures, tests_run])
		quit(1)

func _test_exact_zero_fast_path() -> void:
	var field = ChemicalFieldScript.new(64, 64, 1.0, 1.0, 0.0)
	_assert_true(field.is_all_zero(), "new zero field carries exact-zero state")
	var before: PackedFloat64Array = field.values.duplicate()
	for _i in range(100):
		field.step_diffusion(0.1)
	_assert_true(field.values == before, "exact-zero diffusion fast path leaves every lattice value bit-identical")
	_assert_close(field.total_amount(), 0.0, 0.0, "exact-zero total remains exact")
	_assert_close(field.checksum(), 0.0, 0.0, "exact-zero checksum remains exact")
	field.add_amount(4, 4, 1e-12)
	_assert_true(not field.is_all_zero(), "any positive molecule invalidates zero fast path")

func _test_parallel_diffusion_matches_serial_exactly() -> void:
	var serial = _fixture_world()
	var parallel = _fixture_world()
	for _step in range(40):
		serial.diffuse(0.1, false)
		parallel.diffuse(0.1, true)
	_assert_close(parallel.checksum(), serial.checksum(), 0.0, "parallel independent-field diffusion preserves exact world checksum")
	for field_name in serial.field_order:
		_assert_true(
			serial.get_field(field_name).values == parallel.get_field(field_name).values,
			"parallel diffusion preserves every %s lattice value" % field_name
		)

func _fixture_world():
	var world = WorldStateScript.new(16, 16, 1.0)
	world.register_field("a", 1.0, 2.0)
	world.register_field("b", 0.7, 0.0)
	world.register_field("c", 0.35, 1.25)
	world.get_field("a").remove_amount(8, 8, 0.75)
	world.get_field("b").add_amount(3, 12, 2.5)
	world.get_field("c").add_amount(13, 2, 0.4)
	return world

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
