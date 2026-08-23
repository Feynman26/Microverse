extends SceneTree

const Experiment = preload("res://src/experiments/m5c_amplified_circuit.gd")

const MATERIAL_EFFECT: float = 0.05
const T_CRITICAL_DF31: float = 2.04
const REQUIRED_SIGN_CONSISTENCY: int = 21

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	# Prospectively frozen after earlier weaker constructs failed. These seeds have
	# not appeared in any discovery or confirmation panel. The gate is bilateral:
	# the environment may favor or penalize the responsive architecture, but the
	# direction must be reproducible and materially nonzero.
	var values: Array = []
	var all_divided: bool = true
	var positive: int = 0
	var negative: int = 0
	for seed in range(18001, 18033):
		var record: Dictionary = Experiment.paired_differential(seed)
		var value: float = float(record["normalized_differential"])
		values.append(value)
		all_divided = all_divided and bool(record["all_divided"])
		if value > 0.0:
			positive += 1
		elif value < 0.0:
			negative += 1
		print("M5-C AMPLIFIED seed=%d stable_adv=%.8f fluct_adv=%.8f normalized_delta=%.6f" % [
			seed,
			float(record["stable_advantage"]),
			float(record["fluctuating_advantage"]),
			value
		])

	var mean_value: float = _mean(values)
	var sd: float = _sample_sd(values)
	var se: float = sd / sqrt(float(values.size()))
	var signed_t: float = 0.0
	if se > 0.0:
		signed_t = mean_value / se
	elif mean_value > 0.0:
		signed_t = INF
	elif mean_value < 0.0:
		signed_t = -INF
	var same_sign: int = positive if mean_value >= 0.0 else negative

	print("M5-C AMPLIFIED mean=%.6f sd=%.6f se=%.6f signed_t=%.6f positive=%d negative=%d same_sign=%d/%d" % [
		mean_value, sd, se, signed_t, positive, negative, same_sign, values.size()
	])

	_assert_true(all_divided, "all 32 amplified-circuit lineages reach ordinary division criteria")
	_assert_true(absf(mean_value) >= MATERIAL_EFFECT, "amplified generic regulatory circuit produces at least five-percent environment-selection differential")
	_assert_true(absf(signed_t) >= T_CRITICAL_DF31, "32-seed amplified-circuit differential clears the predeclared two-sided t threshold")
	_assert_true(same_sign >= REQUIRED_SIGN_CONSISTENCY, "at least 21 of 32 independent seeds share the observed selection direction")

	if failures == 0:
		print("PASS: %d amplified M5-C confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d amplified M5-C confirmation tests failed" % [failures, tests_run])
		quit(1)

func _mean(values: Array) -> float:
	var total: float = 0.0
	for value_variant in values:
		total += float(value_variant)
	return total / float(values.size())

func _sample_sd(values: Array) -> float:
	if values.size() < 2:
		return 0.0
	var mean_value: float = _mean(values)
	var total: float = 0.0
	for value_variant in values:
		var delta: float = float(value_variant) - mean_value
		total += delta * delta
	return sqrt(total / float(values.size() - 1))

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
