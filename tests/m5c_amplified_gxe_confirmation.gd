extends SceneTree

const Assay = preload("res://src/experiments/m5c_amplified_gxe_assay.gd")

const MAX_TICKS: int = 7200
const MATERIAL_EFFECT: float = 0.05
const T_CRITICAL_DF31: float = 2.04
const REQUIRED_POSITIVE: int = 21

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	var seeds: Array = []
	for seed in range(20001, 20033):
		seeds.append(seed)
	var result: Dictionary = Assay.run_panel(seeds, MAX_TICKS)
	var values: Array = result["normalized_differentials"]
	var mean_value: float = float(result["mean_normalized_differential"])
	var sd: float = float(result["sample_sd"])
	var se: float = sd / sqrt(float(values.size()))
	var signed_t: float = 0.0
	if se > 0.0:
		signed_t = mean_value / se
	elif mean_value > 0.0:
		signed_t = INF
	var positive: int = 0
	for record_variant in result["seed_results"]:
		var record: Dictionary = record_variant
		var value: float = float(record["normalized_differential"])
		if value > 0.0:
			positive += 1
		print("M5-C AMP-GXE seed=%d stable_adv=%.8f high_anoxic_adv=%.8f normalized_delta=%.6f" % [
			int(record["seed"]), float(record["stable_advantage"]), float(record["fluctuating_advantage"]), value
		])
	print("M5-C AMP-GXE mean=%.6f sd=%.6f se=%.6f signed_t=%.6f positive=%d/%d" % [mean_value, sd, se, signed_t, positive, values.size()])

	_assert_true(bool(result["all_divided"]), "all amplified GxE lineages reach ordinary division criteria before timeout")
	_assert_true(mean_value >= MATERIAL_EFFECT, "amplified GxE interaction improves responsive relative growth by at least five percent")
	_assert_true(signed_t >= T_CRITICAL_DF31, "32-seed amplified GxE interaction clears the predeclared two-sided t threshold")
	_assert_true(positive >= REQUIRED_POSITIVE, "at least 21 of 32 independent seeds reproduce the predicted positive interaction")

	if failures == 0:
		print("PASS: %d amplified GxE M5-C confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d amplified GxE M5-C confirmation tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
