extends SceneTree

const ReceptorSystemScript = preload("res://src/sensing/receptor_system.gd")
const SignallingSystemScript = preload("res://src/signalling/signalling_system.gd")
const MotorSystemScript = preload("res://src/motility/motor_system.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")

const RECEPTOR_SIGNATURE: int = 0xE000
const MOTOR_SIGNATURE: int = 0xF000

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_receptor_occupancy_depends_on_concentration_and_affinity()
	_test_receptor_maintenance_has_explicit_cost()
	_test_signalling_state_persists_then_decays()
	_test_motor_cannot_move_without_atp()
	_test_same_seed_replays_heading_switches()
	_test_generic_temporal_circuit_biases_displacement_only_when_environment_is_correlated()
	if failures == 0:
		print("PASS: %d M11 sensing/signalling/motor tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M11 sensing/signalling/motor tests failed" % [failures, tests_run])
		quit(1)

func _expression_with(signature: int, abundance: float) -> Dictionary:
	return {1: {"mrna": {}, "protein": {signature: abundance}}}

func _test_receptor_occupancy_depends_on_concentration_and_affinity() -> void:
	var expression: Dictionary = _expression_with(RECEPTOR_SIGNATURE, 100.0)
	var low: Dictionary = ReceptorSystemScript.occupancy(expression, {"G": 0.1}, 1.0, 160.0, 0.5)
	var high: Dictionary = ReceptorSystemScript.occupancy(expression, {"G": 5.0}, 1.0, 160.0, 0.5)
	var mismatch: Dictionary = ReceptorSystemScript.occupancy(expression, {"P": 5.0}, 1.0, 160.0, 0.5)
	_assert_true(ReceptorSystemScript.has_receptor_localization(RECEPTOR_SIGNATURE), "E-domain protein is recognized only by generic membrane-localization motif")
	_assert_true(float(high["bound_total"]) > float(low["bound_total"]), "receptor occupancy increases with compatible extracellular concentration")
	_assert_true(float(low["bound_total"]) > 0.0, "compatible ligand produces finite occupancy")
	_assert_true(float(mismatch["bound_total"]) < float(high["bound_total"]), "digital sequence affinity changes occupancy at equal ligand concentration")
	_assert_true(not ReceptorSystemScript.has_receptor_localization(0x1357), "ancestral M10 protein is not silently granted receptor localization")

func _test_receptor_maintenance_has_explicit_cost() -> void:
	var cost: float = ReceptorSystemScript.maintenance_cost(100.0, 0.1, 0.002)
	_assert_close(cost, 0.02, 1e-12, "realized receptor inventory carries explicit per-time ATP maintenance cost")
	_assert_close(ReceptorSystemScript.maintenance_cost(0.0, 10.0, 1.0), 0.0, 1e-12, "no receptor cohort has no receptor maintenance charge")

func _test_signalling_state_persists_then_decays() -> void:
	var state: Dictionary = {}
	SignallingSystemScript.step(state, {RECEPTOR_SIGNATURE: 80.0}, 1.0, 1.5, 0.2)
	var driven: float = float(state.get(RECEPTOR_SIGNATURE, 0.0))
	_assert_true(driven > 0.0 and driven < 80.0, "occupied receptor creates a bounded activated molecular pool")
	SignallingSystemScript.step(state, {}, 0.5, 1.5, 0.2)
	var persisted: float = float(state.get(RECEPTOR_SIGNATURE, 0.0))
	_assert_true(persisted > 0.0 and persisted < driven, "activated molecular state persists after input removal while decaying")
	for _i in range(100):
		SignallingSystemScript.step(state, {}, 0.5, 1.5, 1.0)
	_assert_true(SignallingSystemScript.total_active(state) < 1e-10, "finite activation state decays away without a dedicated memory variable")

func _test_motor_cannot_move_without_atp() -> void:
	var motors: Dictionary = {MOTOR_SIGNATURE: 1.0}
	var request: Dictionary = MotorSystemScript.movement_request(motors, Vector2.RIGHT, 0.5, 2.0)
	var unfunded: Dictionary = MotorSystemScript.funded_displacement(request["requested_displacement"], 0.0, 0.5)
	var funded: Dictionary = MotorSystemScript.funded_displacement(request["requested_displacement"], 10.0, 0.5)
	_assert_true(float(request["requested_distance"]) > 0.0, "realized motor cohort proposes finite active displacement")
	_assert_close(float(unfunded["actual_distance"]), 0.0, 1e-15, "zero ATP produces exactly zero motor displacement")
	_assert_close(float(unfunded["atp_spent"]), 0.0, 1e-15, "zero ATP cannot be debited below zero")
	_assert_true(float(funded["actual_distance"]) > 0.0 and float(funded["atp_spent"]) > 0.0, "funded motor displacement consumes explicit ATP")
	_assert_close(float(funded["atp_spent"]), float(funded["actual_distance"]) * 0.5, 1e-12, "movement ATP cost scales exactly with realized distance")

func _test_same_seed_replays_heading_switches() -> void:
	var first: Array = _heading_history(81271)
	var second: Array = _heading_history(81271)
	_assert_true(first == second, "same molecular state and RNG seed replay exact stochastic heading trajectory")

func _heading_history(seed: int) -> Array:
	var rng = DeterministicRngScript.new(seed)
	var heading := Vector2.RIGHT
	var history: Array = []
	for _i in range(30):
		var update: Dictionary = MotorSystemScript.update_heading(heading, 0.0, 0.2, rng, 1.5, 1.0)
		heading = update["heading"]
		history.append([heading.x, heading.y, update["turned"]])
	return history

func _test_generic_temporal_circuit_biases_displacement_only_when_environment_is_correlated() -> void:
	var correlated_total: float = 0.0
	var homogeneous_total: float = 0.0
	for seed in [61001, 61002, 61003, 61004, 61005, 61006, 61007, 61008]:
		correlated_total += _run_temporal_control(seed, true)
		homogeneous_total += _run_temporal_control(seed, false)
	var correlated_mean: float = correlated_total / 8.0
	var homogeneous_mean: float = homogeneous_total / 8.0
	_assert_true(correlated_mean > homogeneous_mean + 0.5, "generic fast-vs-slow molecular coupling produces greater up-gradient displacement than homogeneous control")

# A compact hand-built control genotype assay. No gradient vector is provided to
# the motor: concentration is sampled only at current x. Same molecular circuit,
# RNG and starting heading are used in both environments. The only difference is
# whether current concentration is correlated with x.
func _run_temporal_control(seed: int, correlated: bool) -> float:
	var rng = DeterministicRngScript.new(seed)
	var signalling: Dictionary = {}
	var heading := Vector2.RIGHT
	var position := Vector2(2.0, 0.0)
	var receptor_expression: Dictionary = _expression_with(RECEPTOR_SIGNATURE, 160.0)
	var motor_expression: Dictionary = _expression_with(MOTOR_SIGNATURE, 160.0)
	var motors: Dictionary = MotorSystemScript.realized_motors(motor_expression, 160.0)
	for _tick in range(180):
		var concentration: float = clampf(position.x / 20.0, 0.0, 1.0) if correlated else 0.5
		var occupancy: Dictionary = ReceptorSystemScript.occupancy(
			receptor_expression, {"G": concentration}, 1.0, 160.0, 0.15
		)
		var fast: Dictionary = occupancy["by_signature"]
		var drive: float = MotorSystemScript.control_drive(motors, fast, signalling)
		var heading_update: Dictionary = MotorSystemScript.update_heading(heading, drive, 0.1, rng, 1.1, 0.08)
		heading = heading_update["heading"]
		var request: Dictionary = MotorSystemScript.movement_request(motors, heading, 0.1, 1.0)
		position += request["requested_displacement"]
		position.x = clampf(position.x, 0.0, 20.0)
		SignallingSystemScript.step(signalling, fast, 0.1, 0.35, 0.08)
	return position.x

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
