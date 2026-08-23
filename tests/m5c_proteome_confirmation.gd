extends SceneTree

const SpikeScript = preload("res://src/experiments/m5c_timescale_spike.gd")

# Prospectively frozen after the RNG-locality correction and before any output
# from seeds 18001-18048. Earlier panels used the variable-draw Poisson sampler
# and therefore cannot serve as inference for the corrected stochastic model.
const CONFIRM_PHASE_TICKS: int = 200
const CONFIRM_PHASE_MIN: float = 20.0
const MATERIAL_EFFECT: float = 0.02
const TWO_SIDED_T_CRITICAL_DF47: float = 2.01
const REQUIRED_SAME_SIGN_SEEDS: int = 31

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	var seeds: Array = []
	for seed in range(18001, 18049):
		seeds.append(seed)

	var result: Dictionary = SpikeScript.run_sweep(seeds, [CONFIRM_PHASE_TICKS])
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

	var positive_count: int = 0
	var negative_count: int = 0
	for i in range(record["seed_results"].size()):
		var seed_result: Dictionary = record["seed_results"][i]
		var value: float = float(seed_result["normalized_differential"])
		if value > 0.0:
			positive_count += 1
		elif value < 0.0:
			negative_count += 1
		print("M5-C STREAM-STABLE seed=%d stable_adv=%.8f fluct_adv=%.8f normalized_delta=%.6f" % [
			int(seed_result["seed"]),
			float(seed_result["stable_advantage"]),
			float(seed_result["fluctuating_advantage"]),
			value
		])

	var same_sign_count: int = positive_count if mean_value > 0.0 else negative_count
	if mean_value == 0.0:
		same_sign_count = 0

	print("M5-C STREAM-STABLE phase_min=%.1f mean=%.6f sd=%.6f se=%.6f signed_t=%.6f positive=%d negative=%d same_sign=%d/%d" % [
		CONFIRM_PHASE_MIN,
		mean_value,
		sample_sd,
		standard_error,
		signed_t,
		positive_count,
		negative_count,
		same_sign_count,
		values.size()
	])

	_assert_true(bool(result["all_divided"]), "all stream-stable finite-proteome M5-C lineages reach ordinary division criteria before timeout")
	_assert_true(absf(mean_value) >= MATERIAL_EFFECT, "stream-stable finite proteome makes stable-vs-fluctuating regulatory selection differ by at least two percent of baseline growth rate")
	_assert_true(absf(signed_t) >= TWO_SIDED_T_CRITICAL_DF47, "independent 48-seed stream-stable differential clears the predeclared two-sided t threshold")
	_assert_true(same_sign_count >= REQUIRED_SAME_SIGN_SEEDS, "at least 31 of 48 independent stream-stable seeds share the observed mean selection direction")

	if failures == 0:
		print("PASS: %d M5-C stream-stable confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5-C stream-stable confirmation tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
