extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	var control = _create_sim(false)
	var observed = _create_sim(true)
	for _tick in range(40):
		control.step(1)
		observed.step(1)
	_assert_close(observed.checksum(), control.checksum(), 0.0, "profiling preserves exact authoritative checksum")
	_assert_true(observed.event_log == control.event_log, "profiling preserves exact event history")
	_assert_true(observed.tick_index == control.tick_index, "profiling preserves exact clock")
	var report: Dictionary = observed.performance_profiler.report()
	_assert_true(report.has("tick_total"), "enabled profiler records whole-tick measurements")
	_assert_true(int(report["tick_total"]["samples"]) == 40, "profiler records one whole-tick sample per step")
	_assert_true(control.performance_profiler.report().is_empty(), "disabled profiler records no samples")
	observed.performance_profiler.reset()
	_assert_true(observed.performance_profiler.report().is_empty(), "profiler reset removes observational samples")
	if failures == 0:
		print("PASS: %d P0 profiling-equivalence tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d P0 profiling-equivalence tests failed" % [failures, tests_run])
		quit(1)

func _create_sim(profiling: bool):
	var config = SimConfigScript.new()
	config.seed = 903001
	config.world_width = 16
	config.world_height = 16
	config.max_cells = 32
	config.mutation_enabled = false
	config.performance_profiling_enabled = profiling
	config.validate()
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor(Vector2(8.0, 8.0))
	return sim

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
