extends SceneTree

const SpikeScript = preload("res://src/experiments/m5c_timescale_spike.gd")

const CONFIRM_PHASE_TICKS: int = 100
const CONFIRM_PHASE_MIN: float = 10.0
const REQUIRED_NEGATIVE_SEEDS: int = 31
const MATERIAL_NEGATIVE_EFFECT: float = -0.02
const TWO_SIDED_T_CRITICAL_DF47: float = -2.01

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	# Prospectively frozen after the exploratory 15001-15004 sweep nominated
	# the shortest qualifying new timescale (100 ticks = 10 min): mean D_norm
	# -0.142701, sample SD 0.334143, 3/4 negative. The 20-minute powered null
	# (14001-14024) is not reused. A two-sided one-sample t power calculation
	# using |d| ~= 0.427 gives ~82.6% power at n=48, alpha=0.05.
	#
	# Confirmation is directional because the discovery nominated a negative
	# differential. No seed, phase, threshold or sample size may be changed after
	# observing this panel.
	var seeds: Array = []
	for seed in range(16001, 16049):
		seeds.append(seed)
	var result: Dictionary = SpikeScript.run_sweep(
		seeds,
		[CONFIRM_PHASE_TICKS]
	)
	var record: Dictionary = result["by_phase"][CONFIRM_PHASE_TICKS]
	var values: Array = record["normalized_differentials"]
	var mean_value: float = float(record["mean_normalized_differential"])
	var sample_sd: float = float(record["sample_sd"])
	var standard_error: float = sample_sd / sqrt(float(values.size()))
	var signed_t: float = 0.0
	if standard_error > 0.0:
		signed_t = mean_value / standard_error
	elif mean_value < 0.0:
		signed_t = -INF
	elif mean_value > 0.0:
		signed_t = INF

	var negative_count: int = 0
	for i in range(record["seed_results"].size()):
		var seed_result: Dictionary = record["seed_results"][i]
		var value: float = float(seed_result["normalized_differential"])
		if value < 0.0:
			negative_count += 1
		print("M5-C CONFIRM10 seed=%d stable_adv=%.8f fluct_adv=%.8f normalized_delta=%.6f" % [
			int(seed_result["seed"]),
			float(seed_result["stable_advantage"]),
			float(seed_result["fluctuating_advantage"]),
			value
		])

	print("M5-C CONFIRM10 phase_min=%.1f mean=%.6f sd=%.6f se=%.6f signed_t=%.6f negative=%d/%d" % [
		CONFIRM_PHASE_MIN,
		mean_value,
		sample_sd,
		standard_error,
		signed_t,
		negative_count,
		values.size()
	])

	_assert_true(bool(result["all_divided"]), "all 48-seed M5-C 10-minute confirmation lineages reach ordinary division criteria before timeout")
	_assert_true(mean_value <= MATERIAL_NEGATIVE_EFFECT, "10-minute fluctuation reproduces the predeclared negative selection differential with at least two-percent material magnitude")
	_assert_true(signed_t <= TWO_SIDED_T_CRITICAL_DF47, "48-seed mean negative differential clears the predeclared two-sided t threshold")
	_assert_true(negative_count >= REQUIRED_NEGATIVE_SEEDS, "at least 31 of 48 independent seeds reproduce the negative environment-selection direction")

	if failures == 0:
		print("PASS: %d M5-C 10-minute confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5-C 10-minute confirmation tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
