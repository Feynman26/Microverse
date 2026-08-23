extends SceneTree

const AssayScript = preload("res://src/experiments/m5c_gxe_assay.gd")

# Final prospectively frozen M5-C gate. Seeds 19001-19024 have not appeared in
# any prior M5-C analysis. The direction is predeclared positive: under stable
# O2=6 the O2-compatible repressor is largely ligand-inhibited, whereas 40-min
# anoxic phases permit the responsive genotype to reduce R03 allocation and
# redirect a finite shared proteome. No production parameter is changed by this
# test.
const MATERIAL_POSITIVE_EFFECT: float = 0.02
const TWO_SIDED_T_CRITICAL_DF23: float = 2.07
const REQUIRED_POSITIVE_SEEDS: int = 17

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	var seeds: Array = []
	for seed in range(19001, 19025):
		seeds.append(seed)

	var result: Dictionary = AssayScript.run_panel(seeds)
	var values: Array = result["normalized_differentials"]
	var mean_value: float = float(result["mean_normalized_differential"])
	var sample_sd: float = float(result["sample_sd"])
	var standard_error: float = sample_sd / sqrt(float(values.size()))
	var signed_t: float = INF if standard_error <= 0.0 and mean_value > 0.0 else 0.0
	if standard_error > 0.0:
		signed_t = mean_value / standard_error

	var positive_count: int = 0
	for seed_result_variant in result["seed_results"]:
		var seed_result: Dictionary = seed_result_variant
		var value: float = float(seed_result["normalized_differential"])
		if value > 0.0:
			positive_count += 1
		print("M5-C GXE seed=%d stable_high_adv=%.8f high_anoxic_adv=%.8f normalized_delta=%.6f" % [
			int(seed_result["seed"]),
			float(seed_result["stable_advantage"]),
			float(seed_result["fluctuating_advantage"]),
			value
		])

	print("M5-C GXE phase_min=40.0 mean=%.6f sd=%.6f se=%.6f signed_t=%.6f positive=%d/%d" % [
		mean_value,
		sample_sd,
		standard_error,
		signed_t,
		positive_count,
		values.size()
	])

	_assert_true(bool(result["all_divided"]), "all final GxE lineages reach ordinary division criteria before timeout")
	_assert_true(mean_value >= MATERIAL_POSITIVE_EFFECT, "high/anoxic fluctuation improves responsive-vs-constitutive selection by at least two percent of baseline growth rate")
	_assert_true(signed_t >= TWO_SIDED_T_CRITICAL_DF23, "24-seed GxE differential clears the predeclared two-sided t threshold in the predicted positive direction")
	_assert_true(positive_count >= REQUIRED_POSITIVE_SEEDS, "at least 17 of 24 independent seeds reproduce the predicted positive GxE direction")

	if failures == 0:
		print("PASS: %d final M5-C GxE confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d final M5-C GxE confirmation tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
