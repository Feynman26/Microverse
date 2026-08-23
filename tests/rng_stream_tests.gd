extends SceneTree

const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_poisson_call_cardinality_is_lambda_independent()
	_test_zero_lambda_still_preserves_stream_alignment()
	_test_poisson_ensemble_means()
	if failures == 0:
		print("PASS: %d RNG stream-locality tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d RNG stream-locality tests failed" % [failures, tests_run])
		quit(1)

func _test_poisson_call_cardinality_is_lambda_independent() -> void:
	var low = DeterministicRngScript.new(88001)
	var high = DeterministicRngScript.new(88001)
	low.poisson(0.05)
	high.poisson(8.0)
	_assert_close(low.randf(), high.randf(), 0.0, "different Poisson lambdas consume exactly one RNG draw and do not phase-shift the next event")

func _test_zero_lambda_still_preserves_stream_alignment() -> void:
	var zero = DeterministicRngScript.new(88002)
	var nonzero = DeterministicRngScript.new(88002)
	_assert_true(zero.poisson(0.0) == 0, "zero-lambda Poisson remains deterministically zero")
	nonzero.poisson(1.0)
	_assert_close(zero.randf(), nonzero.randf(), 0.0, "zero and nonzero Poisson calls preserve equal stream cardinality")

func _test_poisson_ensemble_means() -> void:
	for lambda_value in [0.1, 2.5, 10.0]:
		var rng = DeterministicRngScript.new(99000 + int(lambda_value * 100.0))
		var trials: int = 50000
		var total: float = 0.0
		for _i in range(trials):
			total += float(rng.poisson(float(lambda_value)))
		var observed: float = total / float(trials)
		var sigma_mean: float = sqrt(float(lambda_value) / float(trials))
		_assert_true(absf(observed - float(lambda_value)) <= 5.0 * sigma_mean + 1e-4, "one-draw Poisson preserves expected ensemble mean for lambda=%.1f (observed=%.5f)" % [lambda_value, observed])

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
