extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")
const ExperimentRunnerScript = preload("res://src/experiments/experiment_runner.gd")
const ExperimentAnalyticsScript = preload("res://src/experiments/experiment_analytics.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_seeded_pulses_are_schedule_local()
	_test_spatial_gradient_and_patch()
	_test_paired_fork_has_identical_prefix()
	_test_serial_transfer_is_deterministic_external_bottleneck()
	_test_batch_summary_and_emergence_detectors_are_observational()
	if failures == 0:
		print("PASS: %d compact M8 environment/analytics closure tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d compact M8 environment/analytics tests failed" % [failures, tests_run])
		quit(1)

func _test_seeded_pulses_are_schedule_local() -> void:
	var schedule: Dictionary = EnvironmentScheduleScript.seeded_pulses(
		{"glucose": 0.5, "oxygen": 0.5},
		{"glucose": 6.0},
		0.35,
		9191
	)
	var first: Array = []
	var second: Array = []
	for tick in range(20):
		first.append(EnvironmentScheduleScript.targets_for_tick(tick, schedule))
		second.append(EnvironmentScheduleScript.targets_for_tick(tick, schedule))
	_assert_true(first == second, "seeded nutrient-pulse schedule replays exactly")
	var has_base: bool = false
	var has_pulse: bool = false
	for target_variant in first:
		var targets: Dictionary = target_variant
		has_base = has_base or absf(float(targets["glucose"]) - 0.5) <= 1e-12
		has_pulse = has_pulse or absf(float(targets["glucose"]) - 6.0) <= 1e-12
	_assert_true(has_base and has_pulse, "fixed seeded-pulse fixture contains both baseline and pulse ticks")

	var config = SimConfigScript.new()
	config.seed = 9199
	config.world_width = 6
	config.world_height = 6
	var sim = SimulationEngineScript.new(config)
	var rng_before: int = int(sim.rng.get_state())
	EnvironmentScheduleScript.apply(sim, 7, schedule)
	_assert_true(int(sim.rng.get_state()) == rng_before, "environmental pulse scheduling does not consume biological RNG state")

func _test_spatial_gradient_and_patch() -> void:
	var config = SimConfigScript.new()
	config.world_width = 5
	config.world_height = 5
	var sim = SimulationEngineScript.new(config)
	var gradient: Dictionary = EnvironmentScheduleScript.linear_gradient(
		"glucose", 0.0, 4.0, "x", {"oxygen": 0.5}
	)
	EnvironmentScheduleScript.apply(sim, 0, gradient)
	_assert_close(sim.world.get_field("glucose").get_value(0, 2), 0.0, 1e-12, "linear gradient sets exact low endpoint")
	_assert_close(sim.world.get_field("glucose").get_value(4, 2), 4.0, 1e-12, "linear gradient sets exact high endpoint")
	_assert_close(sim.world.get_field("glucose").get_value(2, 2), 2.0, 1e-12, "linear gradient interpolates interior lattice values")
	_assert_close(sim.world.get_field("oxygen").get_value(3, 3), 0.5, 1e-12, "spatial schedule can maintain independent uniform companion resources")

	var patch: Dictionary = EnvironmentScheduleScript.patch("glucose", 0.0, 5.0, Vector2(2.0, 2.0), 1.1)
	EnvironmentScheduleScript.apply(sim, 1, patch)
	_assert_close(sim.world.get_field("glucose").get_value(2, 2), 5.0, 1e-12, "resource patch sets its physical center")
	_assert_close(sim.world.get_field("glucose").get_value(0, 0), 0.0, 1e-12, "resource patch leaves distant background at configured value")

func _test_paired_fork_has_identical_prefix() -> void:
	var base: Dictionary = ExperimentRunnerScript.create_spec(
		41401,
		80,
		20,
		EnvironmentScheduleScript.constant({"glucose": 4.0, "oxygen": 0.5, "nitrogen": 3.0, "phosphorus": 2.0})
	)
	base["world_width"] = 8
	base["world_height"] = 8
	base["max_cells"] = 8
	base["mutation_enabled"] = false
	var fork: Dictionary = ExperimentRunnerScript.run_paired_fork(base, 20, {
		"anoxic": EnvironmentScheduleScript.constant({"glucose": 4.0, "oxygen": 0.0, "nitrogen": 3.0, "phosphorus": 2.0}),
		"oxidative": EnvironmentScheduleScript.constant({"glucose": 4.0, "oxygen": 5.0, "nitrogen": 3.0, "phosphorus": 2.0})
	})
	for checksum_variant in fork["prefix_checksums"].values():
		_assert_close(float(checksum_variant), float(fork["fork_checksum"]), 1e-12, "paired environmental arms replay the exact same pre-fork authoritative state")
	var arms: Dictionary = fork["arms"]
	_assert_true(float(arms["anoxic"]["final_checksum"]) != float(arms["oxidative"]["final_checksum"]), "paired arms can diverge only after their controlled environmental fork")
	var comparison: Dictionary = ExperimentAnalyticsScript.compare_paired_fork(fork, "final_population")
	_assert_close(float(comparison["fork_checksum"]), float(fork["fork_checksum"]), 1e-12, "paired-fork analytics preserve provenance of the common prefix")

func _test_serial_transfer_is_deterministic_external_bottleneck() -> void:
	var spec: Dictionary = ExperimentRunnerScript.create_spec(
		51502,
		4,
		1,
		EnvironmentScheduleScript.constant({"glucose": 4.0, "oxygen": 0.5, "nitrogen": 3.0, "phosphorus": 2.0})
	)
	spec["world_width"] = 8
	spec["world_height"] = 8
	spec["max_cells"] = 8
	spec["mutation_enabled"] = false
	spec["initial_positions"] = [Vector2(2, 2), Vector2(5, 2), Vector2(2, 5), Vector2(5, 5)]
	spec["interventions"] = [{"kind": "serial_transfer", "tick": 1, "survivors": 2}]
	var first: Dictionary = ExperimentRunnerScript.run(spec)
	var second: Dictionary = ExperimentRunnerScript.run(spec)
	_assert_true(first == second, "serial-transfer bottleneck is exactly reproducible for fixed experiment seed")
	_assert_true(first["intervention_log"].size() == 1, "serial-transfer intervention is explicitly logged")
	var record: Dictionary = first["intervention_log"][0]
	_assert_true(int(record["before"]) == 4 and int(record["after"]) == 2, "serial transfer samples the configured survivor count without classifying cell fitness")
	_assert_true(int(first["death_causes"].get("serial_transfer", 0)) == 0, "experimental transfer removal is not mislabeled as biological death")

func _test_batch_summary_and_emergence_detectors_are_observational() -> void:
	var synthetic_run: Dictionary = {
		"seed": 1,
		"termination_reason": "horizon",
		"final_population": 5,
		"max_population": 8,
		"mutation_events": 6,
		"division_events": 9,
		"environment": {"mode": "closed"},
		"trajectory": [
			{"population": 2, "genotype_count": 1, "mutation_events": 0, "resources": {"carbon_c2": 5.0}},
			{"population": 6, "genotype_count": 2, "mutation_events": 1, "resources": {"carbon_c2": 4.0}},
			{"population": 1, "genotype_count": 2, "mutation_events": 2, "resources": {"carbon_c2": 3.0}},
			{"population": 7, "genotype_count": 2, "mutation_events": 5, "resources": {"carbon_c2": 2.0}},
			{"population": 5, "genotype_count": 2, "mutation_events": 6, "resources": {"carbon_c2": 1.0}}
		]
	}
	var before: Dictionary = synthetic_run.duplicate(true)
	var detected: Dictionary = ExperimentAnalyticsScript.detect_candidates(synthetic_run)
	_assert_true(synthetic_run == before, "emergence detector cannot mutate the run record it interprets")
	_assert_true(bool(detected["persistent_polymorphism"]["candidate"]), "detector identifies persistent sampled polymorphism as a candidate hypothesis")
	_assert_true(bool(detected["population_bottleneck"]["candidate"]), "detector identifies a severe population bottleneck")
	_assert_true(bool(detected["new_resource_utilization"]["candidate"]), "closed-world resource depletion can be flagged as a resource-use candidate")
	_assert_true(int(detected["recurring_population_cycles"]["turning_points"]) >= 3, "detector counts repeated population direction reversals without feeding them back")

	var second_run: Dictionary = synthetic_run.duplicate(true)
	second_run["seed"] = 2
	second_run["final_population"] = 0
	second_run["termination_reason"] = "extinction"
	second_run["trajectory"][-1]["population"] = 0
	var summary: Dictionary = ExperimentAnalyticsScript.summarize_batch([synthetic_run, second_run])
	_assert_true(int(summary["replicate_count"]) == 2, "batch summary retains explicit replicate count")
	_assert_close(float(summary["extinction_frequency"]), 0.5, 1e-12, "batch summary reports replicate extinction frequency rather than only a mean")
	_assert_true(summary["final_population"].has("median") and summary["final_population"].has("q25") and summary["final_population"].has("q75"), "batch summary reports median and quantiles")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
