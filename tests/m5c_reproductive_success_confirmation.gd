extends SceneTree

const Assay = preload("res://src/experiments/m5c_ros_circuit.gd")
const BaseGxe = preload("res://src/experiments/m5c_gxe_assay.gd")

const MAX_TICKS: int = 7200
const MIN_STABLE_SUCCESS: float = 0.95
const MAX_STABLE_DIFFERENCE: float = 0.03
const MIN_INTERACTION: float = 0.04
const MIN_DISCORDANT: int = 8
const MCNEMAR_Z_CRITICAL: float = 1.96
const WIN_RATIO: int = 3

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	# Independent confirmation panel frozen before observing the ROS-responsive
	# circuit at reproductive scale. No seed from earlier M5-C panels is reused.
	var seeds: Array = []
	for seed in range(23001, 23065):
		seeds.append(seed)

	var stable_c_success: int = 0
	var stable_r_success: int = 0
	var fluct_c_success: int = 0
	var fluct_r_success: int = 0
	var responsive_only: int = 0
	var constitutive_only: int = 0
	var phase_responsive_only: Dictionary = {}
	var phase_constitutive_only: Dictionary = {}
	var offsets: Array = [0, BaseGxe.PHASE_TICKS / 2, BaseGxe.PHASE_TICKS, (3 * BaseGxe.PHASE_TICKS) / 2]
	for offset_variant in offsets:
		phase_responsive_only[int(offset_variant)] = 0
		phase_constitutive_only[int(offset_variant)] = 0

	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var stable_c: Dictionary = Assay.run_lineage(seed, "constitutive", BaseGxe.CONDITION_HIGH_STABLE, 0, MAX_TICKS)
		var stable_r: Dictionary = Assay.run_lineage(seed, "responsive", BaseGxe.CONDITION_HIGH_STABLE, 0, MAX_TICKS)
		if bool(stable_c["reached_division"]): stable_c_success += 1
		if bool(stable_r["reached_division"]): stable_r_success += 1

		for offset_variant in offsets:
			var offset: int = int(offset_variant)
			var fc: Dictionary = Assay.run_lineage(seed, "constitutive", BaseGxe.CONDITION_HIGH_ANOXIC, offset, MAX_TICKS)
			var fr: Dictionary = Assay.run_lineage(seed, "responsive", BaseGxe.CONDITION_HIGH_ANOXIC, offset, MAX_TICKS)
			var c_success: bool = bool(fc["reached_division"])
			var r_success: bool = bool(fr["reached_division"])
			if c_success: fluct_c_success += 1
			if r_success: fluct_r_success += 1
			if r_success and not c_success:
				responsive_only += 1
				phase_responsive_only[offset] = int(phase_responsive_only[offset]) + 1
			elif c_success and not r_success:
				constitutive_only += 1
				phase_constitutive_only[offset] = int(phase_constitutive_only[offset]) + 1

	var stable_n: float = float(seeds.size())
	var fluct_n: float = float(seeds.size() * offsets.size())
	var p_c_stable: float = float(stable_c_success) / stable_n
	var p_r_stable: float = float(stable_r_success) / stable_n
	var p_c_fluct: float = float(fluct_c_success) / fluct_n
	var p_r_fluct: float = float(fluct_r_success) / fluct_n
	var stable_difference: float = p_r_stable - p_c_stable
	var fluct_difference: float = p_r_fluct - p_c_fluct
	var interaction: float = fluct_difference - stable_difference
	var discordant: int = responsive_only + constitutive_only
	var signed_mcnemar_z: float = 0.0
	if discordant > 0:
		var raw_delta: int = responsive_only - constitutive_only
		var corrected_magnitude: float = maxf(0.0, float(abs(raw_delta) - 1)) / sqrt(float(discordant))
		signed_mcnemar_z = corrected_magnitude if raw_delta >= 0 else -corrected_magnitude

	print("M5-C ROS-REPRO stable C=%d/%d (%.6f) R=%d/%d (%.6f) diff=%.6f" % [stable_c_success, seeds.size(), p_c_stable, stable_r_success, seeds.size(), p_r_stable, stable_difference])
	print("M5-C ROS-REPRO fluctuating C=%d/%d (%.6f) R=%d/%d (%.6f) diff=%.6f interaction=%.6f" % [fluct_c_success, int(fluct_n), p_c_fluct, fluct_r_success, int(fluct_n), p_r_fluct, fluct_difference, interaction])
	print("M5-C ROS-REPRO discordant responsive_only=%d constitutive_only=%d total=%d corrected_z=%.6f" % [responsive_only, constitutive_only, discordant, signed_mcnemar_z])
	for offset_variant in offsets:
		var offset: int = int(offset_variant)
		print("M5-C ROS-REPRO phase offset=%d responsive_only=%d constitutive_only=%d" % [offset, int(phase_responsive_only[offset]), int(phase_constitutive_only[offset])])

	_assert_true(p_c_stable >= MIN_STABLE_SUCCESS and p_r_stable >= MIN_STABLE_SUCCESS, "both regulatory architectures reproduce in at least 95 percent of stable-high trials")
	_assert_true(absf(stable_difference) <= MAX_STABLE_DIFFERENCE, "stable-high reproductive-success difference remains within three percentage points")
	_assert_true(interaction >= MIN_INTERACTION, "environment structure changes responsive relative reproductive success by at least four percentage points")
	_assert_true(discordant >= MIN_DISCORDANT, "paired fluctuating assay contains at least eight informative discordant outcomes")
	_assert_true(signed_mcnemar_z >= MCNEMAR_Z_CRITICAL, "paired reproductive-success difference clears the predeclared continuity-corrected McNemar threshold")
	_assert_true(responsive_only >= WIN_RATIO * maxi(1, constitutive_only), "responsive-only reproductive successes outnumber constitutive-only successes by at least three to one")

	if failures == 0:
		print("PASS: %d final M5-C ROS-responsive reproductive-success confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d final M5-C ROS-responsive reproductive-success confirmation tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
