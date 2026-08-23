extends SceneTree

const SpikeScript = preload("res://src/experiments/m5c_timescale_spike.gd")
const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	# The powered 20-minute panel (seeds 14001-14024) was prospectively frozen
	# and decisively null: mean D_norm=-0.008394, SD=0.304364, SE=0.062128,
	# |t|=0.135113. All trajectories divided. We therefore do not add more seeds
	# to that hypothesis or relax its threshold.
	#
	# This is a new exploratory question: can the same generic regulatory circuit
	# matter at a different environmental timescale? The sweep and nomination rule
	# below are committed before seeds 15001-15004 are observed.
	#
	# Candidate half-periods: 10, 20, 40, 80 biological minutes. The 20-minute
	# condition is retained only as a calibration control because it already has a
	# powered null. A NEW timescale is eligible for prospective confirmation only
	# if |mean D_norm| >= 0.10 and at least 3/4 seeds share the mean sign. If more
	# than one new timescale qualifies, nominate the shortest qualifying phase.
	# If none qualifies, timescale alone is considered insufficient and M5-C must
	# next test an explicit expression/proteome resource trade-off.
	var seeds: Array = [15001, 15002, 15003, 15004]
	var phase_ticks_values: Array = [100, 200, 400, 800]
	var result: Dictionary = SpikeScript.run_sweep(seeds, phase_ticks_values)
	_assert_true(bool(result["all_divided"]), "all M5-C timescale-spike lineages reach ordinary division criteria")

	var nominated_phase: int = -1
	for phase_variant in phase_ticks_values:
		var phase_ticks: int = int(phase_variant)
		var record: Dictionary = result["by_phase"][phase_ticks]
		var values: Array = record["normalized_differentials"]
		var mean_value: float = float(record["mean_normalized_differential"])
		var direction: float = 1.0 if mean_value > 0.0 else -1.0
		var same_sign: int = 0
		for value_variant in values:
			if float(value_variant) * direction > 0.0:
				same_sign += 1
		print("M5-C TIMESCALE phase_ticks=%d phase_min=%.1f mean=%.6f sd=%.6f same_sign=%d/%d values=%s" % [
			phase_ticks,
			float(record["phase_minutes"]),
			mean_value,
			float(record["sample_sd"]),
			same_sign,
			values.size(),
			str(values)
		])
		if phase_ticks != BaseExperiment.PHASE_TICKS and nominated_phase < 0:
			if absf(mean_value) >= 0.10 and same_sign >= 3:
				nominated_phase = phase_ticks

	print("M5-C TIMESCALE nominated_phase_ticks=%d" % nominated_phase)
	# Intentional red sentinel: discovery output must be reviewed and followed by
	# a separately predeclared confirmation. A green discovery run is not M5-C.
	_assert_true(false, "M5-C timescale spike is exploratory; inspect outputs and predeclare an independent confirmation before acceptance")

	if failures == 0:
		print("PASS: %d M5-C timescale-spike tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5-C timescale-spike tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
