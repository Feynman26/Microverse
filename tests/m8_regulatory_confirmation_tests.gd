extends SceneTree

const ConfirmationScript = preload("res://src/experiments/m8_regulatory_confirmation.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_assert_true(ConfirmationScript.CONFIRMATORY_SEEDS.size() == 8, "confirmatory seed panel is frozen at eight untouched replicates")
	_assert_true(ConfirmationScript.HORIZON_TICKS == 1200, "confirmatory outcome uses a predeclared fixed horizon")
	_assert_true(ConfirmationScript.PRIMARY_ENDPOINT.contains("cell-count"), "primary endpoint is descendant cell-count based rather than first-division normalized")

	var replay_a: Dictionary = ConfirmationScript.run_condition(ConfirmationScript.CONFIRMATORY_SEEDS[0], "stable")
	var replay_b: Dictionary = ConfirmationScript.run_condition(ConfirmationScript.CONFIRMATORY_SEEDS[0], "stable")
	_assert_true(replay_a == replay_b, "same confirmatory seed and condition reproduce the exact fixed-horizon lineage outcome")

	var panel: Dictionary = ConfirmationScript.run_confirmatory_panel()
	print("M8 REGULATORY CONFIRMATION: interpretation=%s median=%s q25=%s q75=%s positive=%d zero=%d negative=%d stable_extinctions=%d fluctuating_extinctions=%d effects=%s" % [
		String(panel["interpretation"]),
		str(panel["median_effect"]),
		str(panel["q25_effect"]),
		str(panel["q75_effect"]),
		int(panel["positive_replicates"]),
		int(panel["zero_replicates"]),
		int(panel["negative_replicates"]),
		int(panel["stable_extinctions"]),
		int(panel["fluctuating_extinctions"]),
		str(panel["effects"])
	])
	_assert_true(String(panel["panel_status"]) == "frozen_confirmatory", "panel labels itself as confirmatory rather than exploratory")
	_assert_true(panel["rows"].size() == ConfirmationScript.CONFIRMATORY_SEEDS.size(), "no confirmatory seed is silently dropped")
	_assert_true(panel["effects"].size() == ConfirmationScript.CONFIRMATORY_SEEDS.size(), "every paired seed contributes one primary endpoint")
	_assert_true(int(panel["positive_replicates"]) + int(panel["zero_replicates"]) + int(panel["negative_replicates"]) == ConfirmationScript.CONFIRMATORY_SEEDS.size(), "positive zero and negative replicate counts partition the complete frozen panel")
	_assert_true(panel.has("median_effect") and panel.has("q25_effect") and panel.has("q75_effect"), "confirmation reports median and quantiles instead of only a mean")
	_assert_true(String(panel["interpretation"]) in [
		"supports_environment_dependent_responsive_advantage",
		"evidence_against_environment_dependent_responsive_advantage",
		"inconclusive"
	], "predeclared rule maps the complete panel to one allowed interpretation without retuning")
	for row_variant in panel["rows"]:
		var row: Dictionary = row_variant
		for condition_name in ["stable", "fluctuating"]:
			var result: Dictionary = row[condition_name]
			_assert_true(bool(result["extinct"]) or int(result["realized_ticks"]) == ConfirmationScript.HORIZON_TICKS, "non-extinct confirmatory outcomes reach the exact fixed horizon while extinctions remain terminal outcomes")

	if failures == 0:
		print("PASS: %d M8 regulatory-confirmation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M8 regulatory-confirmation tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
