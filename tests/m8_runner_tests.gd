extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")
const ExperimentRunnerScript = preload("res://src/experiments/experiment_runner.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_environment_schedule_executes_exact_ticks()
	_test_same_spec_and_seed_replay_exactly()
	_test_sampling_is_observational_only()
	_test_extinction_is_terminal_and_explained()
	_test_anoxic_carbon_fed_environment_can_exceed_sixteen_cells()
	_test_batch_order_does_not_change_individual_runs()
	_test_run_metadata_is_self_identifying()

	if failures == 0:
		print("PASS: %d M8 runner/environment tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M8 runner/environment tests failed" % [failures, tests_run])
		quit(1)

func _test_environment_schedule_executes_exact_ticks() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	var sim = SimulationEngineScript.new(config)
	var square: Dictionary = EnvironmentScheduleScript.square_wave(
		{"glucose": 6.0, "oxygen": 5.0},
		{"glucose": 0.5, "oxygen": 0.0},
		3
	)
	EnvironmentScheduleScript.apply(sim, 2, square)
	_assert_close(sim.world.get_field("glucose").minimum_value(), 6.0, 1e-12, "square-wave high phase applies through tick 2")
	_assert_close(sim.world.get_field("oxygen").maximum_value(), 5.0, 1e-12, "high oxygen target is applied exactly")
	EnvironmentScheduleScript.apply(sim, 3, square)
	_assert_close(sim.world.get_field("glucose").maximum_value(), 0.5, 1e-12, "square-wave switches exactly at phase boundary tick")
	_assert_close(sim.world.get_field("oxygen").minimum_value(), 0.0, 1e-12, "low oxygen phase is applied exactly")

func _test_same_spec_and_seed_replay_exactly() -> void:
	var spec: Dictionary = _short_constant_spec(44021, 500, 25)
	var first: Dictionary = ExperimentRunnerScript.run(spec)
	var second: Dictionary = ExperimentRunnerScript.run(spec)
	_assert_true(first == second, "same experiment specification and seed reproduce exact sampled trajectory and final state")

func _test_sampling_is_observational_only() -> void:
	var dense: Dictionary = _short_constant_spec(55102, 500, 1)
	var sparse: Dictionary = dense.duplicate(true)
	sparse["sample_every_ticks"] = 50
	var dense_result: Dictionary = ExperimentRunnerScript.run(dense)
	var sparse_result: Dictionary = ExperimentRunnerScript.run(sparse)
	_assert_close(float(dense_result["final_checksum"]), float(sparse_result["final_checksum"]), 1e-12, "analytics sampling cadence does not feed back into simulation state")
	_assert_true(int(dense_result["final_population"]) == int(sparse_result["final_population"]), "sampling cadence leaves final population unchanged")

func _test_extinction_is_terminal_and_explained() -> void:
	var spec: Dictionary = ExperimentRunnerScript.create_spec(
		66203,
		1400,
		100,
		EnvironmentScheduleScript.closed()
	)
	spec["world_width"] = 12
	spec["world_height"] = 12
	spec["max_cells"] = 16
	spec["mutation_enabled"] = false
	spec["initial_resources"] = {
		"glucose": 0.0,
		"oxygen": 0.0,
		"nitrogen": 0.0,
		"phosphorus": 0.0
	}
	var result: Dictionary = ExperimentRunnerScript.run(spec)
	_assert_true(String(result["termination_reason"]) == "extinction", "resource-free run terminates explicitly by extinction")
	_assert_true(int(result["final_population"]) == 0, "extinct experiment remains extinct and is never silently reseeded")
	_assert_true(int(result["death_causes"].get("energy_failure", 0)) >= 1, "runner reports physiological cause of terminal population loss")
	_assert_true(int(result["realized_ticks"]) < int(result["horizon_ticks"]), "stop-on-extinction avoids meaningless post-extinction ticks")

func _test_anoxic_carbon_fed_environment_can_exceed_sixteen_cells() -> void:
	# The O2=5 and O2=0.5 characterization runs retained large external nutrient
	# inventories but accumulated ROS-driven deaths. The ancestral reaction network
	# already contains fermentative redox relief, so complete oxygen limitation is
	# a clean environmental counterfactual: no physiology constant is changed and
	# the >16 criterion remains fixed. The run stops immediately on cell 17.
	var spec: Dictionary = ExperimentRunnerScript.create_spec(
		77304,
		8000,
		200,
		EnvironmentScheduleScript.constant({
			"glucose": 4.0,
			"oxygen": 0.0,
			"nitrogen": 3.0,
			"phosphorus": 2.0
		})
	)
	spec["world_width"] = 8
	spec["world_height"] = 8
	spec["max_cells"] = 24
	spec["mutation_enabled"] = false
	spec["stop_population_at_least"] = 17
	var result: Dictionary = ExperimentRunnerScript.run(spec)
	print("M8 anoxic-reservoir diagnostic: max_population=%d final_population=%d generation=%d divisions=%d deaths=%s ticks=%d final_resources=%s" % [
		int(result["max_population"]),
		int(result["final_population"]),
		int(result["final_generation"]),
		int(result["division_events"]),
		str(result["death_causes"]),
		int(result["realized_ticks"]),
		str(result["final_resources"])
	])
	_assert_true(
		int(result["max_population"]) > 16,
		"anoxic carbon-fed reservoir allows the same core physiology to grow beyond the observed 16-cell plateau (actual max=%d, generation=%d, divisions=%d, deaths=%s)" % [
			int(result["max_population"]), int(result["final_generation"]), int(result["division_events"]), str(result["death_causes"])
		]
	)
	_assert_true(int(result["final_population"]) >= 17, "population threshold stops only after a living population crosses sixteen cells")
	_assert_true(String(result["termination_reason"]) == "population_threshold", "runner records the population-threshold stop reason explicitly")
	_assert_true(int(result["realized_ticks"]) < int(result["horizon_ticks"]), "population threshold avoids unnecessary post-gate simulation ticks")

func _test_batch_order_does_not_change_individual_runs() -> void:
	var spec: Dictionary = _short_constant_spec(1, 400, 40)
	var forward: Dictionary = ExperimentRunnerScript.run_batch(spec, [801, 802])
	var reverse: Dictionary = ExperimentRunnerScript.run_batch(spec, [802, 801])
	for seed in [801, 802]:
		var first: Dictionary = forward["by_seed"][seed]
		var second: Dictionary = reverse["by_seed"][seed]
		_assert_close(float(first["final_checksum"]), float(second["final_checksum"]), 1e-12, "batch order leaves seed %d final state unchanged" % seed)
		_assert_true(first["trajectory"] == second["trajectory"], "batch order leaves seed %d sampled history unchanged" % seed)

func _test_run_metadata_is_self_identifying() -> void:
	var spec: Dictionary = _short_constant_spec(99005, 50, 10)
	var result: Dictionary = ExperimentRunnerScript.run(spec)
	_assert_true(int(result["schema_version"]) == ExperimentRunnerScript.SCHEMA_VERSION, "run records experiment schema version")
	_assert_true(String(result["model_version"]) == ExperimentRunnerScript.MODEL_VERSION, "run records model version")
	_assert_true(int(result["seed"]) == 99005, "run records exact RNG seed")
	_assert_true(result.has("death_causes") and result.has("final_resources"), "run records diagnostic population-loss and resource context")

func _short_constant_spec(seed: int, horizon: int, cadence: int) -> Dictionary:
	var spec: Dictionary = ExperimentRunnerScript.create_spec(
		seed,
		horizon,
		cadence,
		EnvironmentScheduleScript.constant({
			"glucose": 4.0,
			"oxygen": 5.0,
			"nitrogen": 3.0,
			"phosphorus": 2.0
		})
	)
	spec["world_width"] = 16
	spec["world_height"] = 16
	spec["max_cells"] = 16
	spec["mutation_enabled"] = true
	return spec

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
