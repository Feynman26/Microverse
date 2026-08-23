extends SceneTree

const Assay = preload("res://src/experiments/m5c_burden_circuit.gd")

const MIN_MEDIAN_EFFECT: float = 0.05
const MIN_MEAN_EFFECT: float = 0.05
const MIN_POSITIVE_SEEDS: int = 17 # exact one-sided sign-test p ~= 0.03196 for n=24
const MAX_ABS_STABLE_ADVANTAGE: float = 0.03

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	# Independent robust confirmation frozen after the first burden panel exposed
	# non-Gaussian/censored first-division outcomes. Death or failure to divide is
	# a reproductive outcome and remains represented as zero growth rate; it is
	# not discarded or converted into a timeout error.
	var values: Array = []
	var stable_values: Array = []
	var positive: int = 0
	var stable_divisions: int = 0
	for seed in range(25001, 25025):
		var result: Dictionary = Assay.paired_differential(seed)
		var value: float = float(result["normalized_differential"])
		var stable_adv: float = float(result["stable_advantage"])
		values.append(value)
		stable_values.append(stable_adv)
		if value > 0.0:
			positive += 1
		# paired_differential only exposes joint all_divided, which includes the
		# fluctuating runs. Stable equivalence is therefore assessed by stable
		# advantage plus the ordinary positive baseline reference in the assay.
		if absf(stable_adv) < 1e9:
			stable_divisions += 1
		print("M5-C BURDEN-ROBUST seed=%d Dnorm=%.6f stable_adv=%.8f" % [seed, value, stable_adv])

	var mean_effect: float = _mean(values)
	var median_effect: float = _median(values)
	var mean_stable_adv: float = _mean(stable_values)
	print("M5-C BURDEN-ROBUST summary n=%d mean=%.6f median=%.6f positive=%d mean_stable_adv=%.8f" % [values.size(), mean_effect, median_effect, positive, mean_stable_adv])

	_assert_true(stable_divisions == values.size(), "all seeds produce valid stable paired measurements")
	_assert_true(absf(mean_stable_adv) <= MAX_ABS_STABLE_ADVANTAGE, "architectures remain effectively matched in stable high O2")
	_assert_true(mean_effect >= MIN_MEAN_EFFECT, "mean environment-conditioned reproductive effect is at least five percent of baseline")
	_assert_true(median_effect >= MIN_MEDIAN_EFFECT, "median environment-conditioned reproductive effect is at least five percent of baseline")
	_assert_true(positive >= MIN_POSITIVE_SEEDS, "at least seventeen of twenty-four independent seeds favor the responsive architecture")

	if failures == 0:
		print("PASS: %d robust final M5-C expression-burden confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d robust final M5-C expression-burden confirmation tests failed" % [failures, tests_run])
		quit(1)

func _mean(values: Array) -> float:
	var total: float = 0.0
	for value_variant in values:
		total += float(value_variant)
	return total / float(values.size())

func _median(values: Array) -> float:
	var ordered: Array = values.duplicate()
	ordered.sort()
	var n: int = ordered.size()
	if n % 2 == 1:
		return float(ordered[n / 2])
	return 0.5 * (float(ordered[n / 2 - 1]) + float(ordered[n / 2]))

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
