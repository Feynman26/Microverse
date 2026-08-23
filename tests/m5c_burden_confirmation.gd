extends SceneTree

const Assay = preload("res://src/experiments/m5c_burden_circuit.gd")

const MIN_MEAN_EFFECT: float = 0.05
const T_CRITICAL: float = 2.07 # bilateral ~95%, df=23
const MIN_POSITIVE_SEEDS: int = 16
const MAX_ABS_STABLE_ADVANTAGE: float = 0.03

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	# Final independent M5-C capability panel. Seeds frozen before observation.
	var values: Array = []
	var stable_values: Array = []
	var all_divided: bool = true
	var positive: int = 0
	for seed in range(24001, 24025):
		var result: Dictionary = Assay.paired_differential(seed)
		var value: float = float(result["normalized_differential"])
		var stable_adv: float = float(result["stable_advantage"])
		values.append(value)
		stable_values.append(stable_adv)
		all_divided = all_divided and bool(result["all_divided"])
		if value > 0.0:
			positive += 1
		print("M5-C BURDEN seed=%d Dnorm=%.6f stable_adv=%.8f" % [seed, value, stable_adv])

	var mean_effect: float = _mean(values)
	var sd: float = _sample_sd(values, mean_effect)
	var se: float = sd / sqrt(float(values.size()))
	var t_value: float = 0.0 if se <= 1e-12 else mean_effect / se
	var mean_stable_adv: float = _mean(stable_values)
	print("M5-C BURDEN summary n=%d mean=%.6f sd=%.6f t=%.6f positive=%d mean_stable_adv=%.8f" % [values.size(), mean_effect, sd, t_value, positive, mean_stable_adv])

	_assert_true(all_divided, "all stable and fluctuating burden lineages reach ordinary division")
	_assert_true(absf(mean_stable_adv) <= MAX_ABS_STABLE_ADVANTAGE, "architectures remain effectively matched in stable high O2")
	_assert_true(mean_effect >= MIN_MEAN_EFFECT, "fluctuating environment changes responsive reproductive rate by at least five percent of baseline")
	_assert_true(t_value >= T_CRITICAL, "paired reproductive-rate interaction clears the predeclared t threshold")
	_assert_true(positive >= MIN_POSITIVE_SEEDS, "at least sixteen of twenty-four independent seeds favor the responsive architecture")

	if failures == 0:
		print("PASS: %d final M5-C expression-burden confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d final M5-C expression-burden confirmation tests failed" % [failures, tests_run])
		quit(1)

func _mean(values: Array) -> float:
	var total: float = 0.0
	for value_variant in values:
		total += float(value_variant)
	return total / float(values.size())

func _sample_sd(values: Array, mean_value: float) -> float:
	var ss: float = 0.0
	for value_variant in values:
		var delta: float = float(value_variant) - mean_value
		ss += delta * delta
	return sqrt(ss / float(values.size() - 1))

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
