extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")
const LegacySimulationBackendScript = preload("res://src/runtime/legacy_simulation_backend.gd")
const SimulationCommandScript = preload("res://src/runtime/simulation_command.gd")
const SimulationDeltaScript = preload("res://src/runtime/simulation_delta.gd")
const VisualSnapshotScript = preload("res://src/observation/visual_snapshot.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_command_delta_and_backend_schemas()
	_test_legacy_adapter_is_exact_for_seed_environment_and_steps()
	_test_serial_transfer_command_is_deterministic()
	_test_structural_telemetry_is_observational()
	_test_visual_snapshot_is_versioned_and_backend_identified()
	if failures == 0:
		print("PASS: %d P2 runtime-contract tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d P2 runtime-contract tests failed" % [failures, tests_run])
		quit(1)

func _test_command_delta_and_backend_schemas() -> void:
	var backend = LegacySimulationBackendScript.new(_config(false))
	var seed_delta: Dictionary = backend.execute(SimulationCommandScript.seed_ancestor())
	SimulationDeltaScript.validate(seed_delta)
	_assert_true(int(seed_delta["schema_version"]) == SimulationDeltaScript.SCHEMA_VERSION, "backend emits the declared delta schema")
	_assert_true(String(seed_delta["command_kind"]) == SimulationCommandScript.SEED_ANCESTOR, "delta identifies the committed command")
	_assert_true(int(seed_delta["population_before"]) == 0 and int(seed_delta["population_after"]) == 1, "seed delta describes the exact population transition")
	var metadata: Dictionary = backend.metadata()
	_assert_true(int(metadata["interface_version"]) == 1, "legacy adapter declares ISimulationBackend version 1")
	_assert_true(String(metadata["backend_id"]) == "m10-legacy-adapter", "legacy backend identity remains explicit")
	_assert_true(int(metadata["command_schema_version"]) == SimulationCommandScript.SCHEMA_VERSION, "backend advertises its command schema")

func _test_legacy_adapter_is_exact_for_seed_environment_and_steps() -> void:
	var direct_config = _config(false)
	var backend_config = _config(false)
	var direct = SimulationEngineScript.new(direct_config)
	var backend = LegacySimulationBackendScript.new(backend_config)
	var schedule: Dictionary = EnvironmentScheduleScript.square_wave(
		{"glucose": 4.0, "oxygen": 2.0, "nitrogen": 3.0, "phosphorus": 2.0},
		{"glucose": 0.5, "oxygen": 0.1, "nitrogen": 3.0, "phosphorus": 2.0},
		3
	)
	direct.seed_ancestor(Vector2(4.0, 4.0))
	backend.execute(SimulationCommandScript.seed_ancestor(Vector2(4.0, 4.0)))
	for tick in range(12):
		EnvironmentScheduleScript.apply(direct, tick, schedule)
		direct.step(1)
		backend.execute(SimulationCommandScript.apply_environment(tick, schedule))
		backend.execute(SimulationCommandScript.advance_ticks(1))
	var adapted = backend.legacy_inspection_state()
	_assert_close(adapted.checksum(), direct.checksum(), 0.0, "backend commands preserve the exact M10 checksum")
	_assert_true(adapted.event_log == direct.event_log, "backend commands preserve exact semantic event history")
	_assert_true(adapted.tick_index == direct.tick_index and adapted.population_size() == direct.population_size(), "backend and direct engine end at the same tick and population")

func _test_serial_transfer_command_is_deterministic() -> void:
	var first = LegacySimulationBackendScript.new(_config(false))
	var second = LegacySimulationBackendScript.new(_config(false))
	for position in [Vector2(2, 2), Vector2(6, 2), Vector2(2, 6), Vector2(6, 6)]:
		first.execute(SimulationCommandScript.seed_ancestor(position))
		second.execute(SimulationCommandScript.seed_ancestor(position))
	var first_delta: Dictionary = first.execute(SimulationCommandScript.serial_transfer(51502, 1, 2))
	var second_delta: Dictionary = second.execute(SimulationCommandScript.serial_transfer(51502, 1, 2))
	_assert_true(first_delta["details"] == second_delta["details"], "serial-transfer command chooses the same cells for the same seed")
	_assert_true(int(first_delta["population_after"]) == 2, "serial-transfer delta records the configured survivor count")
	_assert_close(first.legacy_inspection_state().checksum(), second.legacy_inspection_state().checksum(), 0.0, "deterministic intervention commands preserve identical authoritative state")

func _test_structural_telemetry_is_observational() -> void:
	var observed = LegacySimulationBackendScript.new(_config(true))
	var control = SimulationEngineScript.new(_config(false))
	observed.execute(SimulationCommandScript.seed_ancestor())
	control.seed_ancestor()
	observed.execute(SimulationCommandScript.advance_ticks(3))
	control.step(3)
	var report: Dictionary = observed.telemetry_report()
	var simulation_report: Dictionary = report["simulation"]
	var timers: Dictionary = simulation_report["phase_timers"]
	var work: Dictionary = simulation_report["work_counters"]
	_assert_true(int(simulation_report["schema_version"]) == 1, "phase timers and work counters share a versioned telemetry envelope")
	_assert_true(int(timers["tick_total"]["samples"]) == 3, "phase timer records one complete sample per tick")
	_assert_true(int(work["ticks"]) == 3, "work counter records authoritative tick count")
	_assert_true(int(work["intracellular_cell_steps"]) == 3, "work counter exposes cell-scale intracellular effort")
	_assert_close(observed.legacy_inspection_state().checksum(), control.checksum(), 0.0, "structural telemetry cannot change scientific state")

func _test_visual_snapshot_is_versioned_and_backend_identified() -> void:
	var backend = LegacySimulationBackendScript.new(_config(false))
	backend.execute(SimulationCommandScript.seed_ancestor())
	var snapshot: Dictionary = backend.capture_visual_snapshot()
	_assert_true(int(snapshot["schema_version"]) == VisualSnapshotScript.SCHEMA_VERSION, "visual snapshot publishes an explicit schema version")
	_assert_true(String(snapshot["snapshot_kind"]) == "visual", "presentation receives a named visual contract")
	_assert_true(String(snapshot["backend"]["backend_id"]) == "m10-legacy-adapter", "visual contract identifies its producing backend")

func _config(profiling: bool):
	var config = SimConfigScript.new()
	config.seed = 920201
	config.world_width = 10
	config.world_height = 10
	config.max_cells = 16
	config.mutation_enabled = false
	config.performance_profiling_enabled = profiling
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
