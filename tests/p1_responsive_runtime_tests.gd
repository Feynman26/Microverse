extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const InteractiveClockScript = preload("res://src/runtime/interactive_clock.gd")
const InteractiveSimulationRuntimeScript = preload("res://src/runtime/interactive_simulation_runtime.gd")

const WAIT_TIMEOUT_USEC: int = 10000000

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_clock_bounds_and_reports_unserved_demand()
	_test_threaded_exact_steps_preserve_authoritative_result()
	_test_compute_guard_pauses_on_complete_tick_boundary()
	if failures == 0:
		print("PASS: %d P1 responsive-runtime tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d P1 responsive-runtime tests failed" % [failures, tests_run])
		quit(1)

func _test_clock_bounds_and_reports_unserved_demand() -> void:
	var clock = InteractiveClockScript.new(0.10)
	clock.set_requested_multiplier(1.0)
	clock.add_wall_time(0.10)
	_assert_true(clock.can_consume_tick(), "1x clock makes one exact tick due after 0.1 real seconds")
	clock.consume_tick()
	_assert_close(clock.backlog_min, 0.0, 1e-12, "consuming a tick removes exactly one biological dt")
	_assert_close(clock.requested_ticks_per_second(), 10.0, 1e-12, "1x honestly reports ten requested ticks per second")

	clock.set_requested_multiplier(100.0)
	clock.add_wall_time(1.0)
	_assert_close(clock.backlog_min, 50.0, 1e-12, "overload backlog is bounded to the declared half-second horizon")
	_assert_close(clock.backlog_ticks(), 500.0, 1e-12, "bounded backlog remains visible in exact tick units")
	_assert_close(clock.unserved_requested_min, 50.0, 1e-12, "discarded wall-clock demand is reported instead of hidden")
	_assert_true(clock.overload_events == 1, "backpressure records an overload event")

func _test_threaded_exact_steps_preserve_authoritative_result() -> void:
	var reference_config = _test_config()
	var reference = SimulationEngineScript.new(reference_config)
	reference.seed_ancestor()
	reference.step(8)

	var runtime = InteractiveSimulationRuntimeScript.new()
	var start_error: Error = runtime.start(_test_config(), true, 1.0)
	_assert_true(start_error == OK, "interactive worker starts headlessly")
	if start_error != OK:
		return

	var initial: Dictionary = _wait_for_tick(runtime, 0)
	_assert_true(not initial.is_empty(), "worker publishes an initial immutable visual snapshot")
	if initial.is_empty():
		runtime.stop()
		return
	var initial_checksum: float = float(initial["checksum"])
	_assert_true(typeof(initial["glucose_values"]) == TYPE_PACKED_FLOAT64_ARRAY, "snapshot owns a packed field copy")
	_assert_true(typeof(initial["cell_positions"]) == TYPE_PACKED_VECTOR2_ARRAY, "snapshot owns packed cell positions")

	_assert_true(runtime.request_exact_ticks(8), "paused runtime accepts exact tick command")
	var completed: Dictionary = _wait_for_tick(runtime, 8)
	_assert_true(not completed.is_empty(), "threaded runtime completes requested exact ticks")
	if not completed.is_empty():
		_assert_true(int(completed["tick"]) == 8, "exact command stops on the requested tick boundary")
		_assert_close(float(completed["checksum"]), reference.checksum(), 0.0, "thread scheduling preserves exact M10 checksum")
		_assert_true(int(completed["event_count"]) == reference.event_log.size(), "thread scheduling preserves event count")
		_assert_true(int(completed["population"]) == reference.population_size(), "thread scheduling preserves population")
	_assert_close(float(initial["checksum"]), initial_checksum, 0.0, "previously published snapshot remains unchanged")

	var command_started_usec: int = Time.get_ticks_usec()
	runtime.set_requested_multiplier(100.0)
	var command_elapsed_usec: int = Time.get_ticks_usec() - command_started_usec
	_assert_true(command_elapsed_usec < 250000, "speed command does not wait for a biological tick")
	var accelerated: Dictionary = _wait_for_multiplier(runtime, 100.0)
	_assert_true(not accelerated.is_empty(), "worker publishes changed requested clock while paused")
	if not accelerated.is_empty():
		var runtime_state: Dictionary = accelerated["runtime"]
		_assert_close(float(runtime_state["requested_ticks_per_second"]), 1000.0, 1e-12, "100x demand is reported honestly")

	var resume_started_usec: int = Time.get_ticks_usec()
	_assert_true(runtime.set_paused(false), "runtime accepts resume below compute guard")
	var resume_elapsed_usec: int = Time.get_ticks_usec() - resume_started_usec
	_assert_true(resume_elapsed_usec < 250000, "resume command does not wait for a biological tick")
	var running: Dictionary = _wait_for_tick(runtime, 9)
	_assert_true(not running.is_empty(), "100x worker advances while the caller remains responsive")
	var snapshot_started_usec: int = Time.get_ticks_usec()
	runtime.latest_snapshot()
	var snapshot_elapsed_usec: int = Time.get_ticks_usec() - snapshot_started_usec
	_assert_true(snapshot_elapsed_usec < 250000, "visual snapshot read does not wait for biological completion")
	var pause_started_usec: int = Time.get_ticks_usec()
	runtime.set_paused(true)
	var pause_elapsed_usec: int = Time.get_ticks_usec() - pause_started_usec
	_assert_true(pause_elapsed_usec < 250000, "pause command does not wait for a biological tick")
	_assert_true(not _wait_for_paused(runtime).is_empty(), "pause becomes visible at a completed tick boundary")

	runtime.stop()

func _test_compute_guard_pauses_on_complete_tick_boundary() -> void:
	var config = _test_config()
	config.max_cells = 1
	config.validate()
	var runtime = InteractiveSimulationRuntimeScript.new()
	var start_error: Error = runtime.start(config, false, 100.0)
	_assert_true(start_error == OK, "compute-guard fixture starts")
	if start_error != OK:
		return
	var stopped: Dictionary = _wait_for_compute_limit(runtime)
	_assert_true(not stopped.is_empty(), "compute guard publishes an explicit terminal runtime state")
	if not stopped.is_empty():
		var runtime_state: Dictionary = stopped["runtime"]
		_assert_true(int(stopped["tick"]) == 1, "compute guard stops after one complete authoritative tick")
		_assert_true(bool(runtime_state["paused"]), "compute guard pauses the worker")
		_assert_true(not runtime.set_paused(false), "compute guard cannot be silently resumed")
	runtime.stop()

func _test_config():
	var config = SimConfigScript.new()
	config.seed = 910201
	config.world_width = 12
	config.world_height = 12
	config.max_cells = 32
	config.mutation_enabled = false
	config.performance_profiling_enabled = false
	config.validate()
	return config

func _wait_for_tick(runtime, target_tick: int) -> Dictionary:
	var deadline_usec: int = Time.get_ticks_usec() + WAIT_TIMEOUT_USEC
	while Time.get_ticks_usec() < deadline_usec:
		var snapshot: Dictionary = runtime.latest_snapshot()
		if not snapshot.is_empty() and int(snapshot.get("tick", -1)) >= target_tick:
			return snapshot
		OS.delay_msec(1)
	return {}

func _wait_for_multiplier(runtime, target_multiplier: float) -> Dictionary:
	var deadline_usec: int = Time.get_ticks_usec() + WAIT_TIMEOUT_USEC
	while Time.get_ticks_usec() < deadline_usec:
		var snapshot: Dictionary = runtime.latest_snapshot()
		var runtime_state: Dictionary = snapshot.get("runtime", {})
		if is_equal_approx(float(runtime_state.get("requested_multiplier", 0.0)), target_multiplier):
			return snapshot
		OS.delay_msec(1)
	return {}

func _wait_for_paused(runtime) -> Dictionary:
	var deadline_usec: int = Time.get_ticks_usec() + WAIT_TIMEOUT_USEC
	while Time.get_ticks_usec() < deadline_usec:
		var snapshot: Dictionary = runtime.latest_snapshot()
		var runtime_state: Dictionary = snapshot.get("runtime", {})
		if bool(runtime_state.get("paused", false)):
			return snapshot
		OS.delay_msec(1)
	return {}

func _wait_for_compute_limit(runtime) -> Dictionary:
	var deadline_usec: int = Time.get_ticks_usec() + WAIT_TIMEOUT_USEC
	while Time.get_ticks_usec() < deadline_usec:
		var snapshot: Dictionary = runtime.latest_snapshot()
		var runtime_state: Dictionary = snapshot.get("runtime", {})
		if bool(runtime_state.get("compute_limit_reached", false)):
			return snapshot
		OS.delay_msec(1)
	return {}

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
