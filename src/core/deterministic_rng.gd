extends RefCounted
class_name DeterministicRng

var _rng := RandomNumberGenerator.new()

func _init(seed_value: int = 1) -> void:
	_rng.seed = seed_value

func reseed(seed_value: int) -> void:
	_rng.seed = seed_value

func randf() -> float:
	return _rng.randf()

func randf_range(from: float, to: float) -> float:
	return _rng.randf_range(from, to)

func randi_range(from: int, to: int) -> int:
	return _rng.randi_range(from, to)

# Exact Knuth sampler for the small event intensities used by M5 gene
# expression. Large lambda is decomposed into independent chunks so the method
# remains numerically well behaved without introducing a second RNG or a
# platform-dependent Gaussian approximation.
func poisson(lambda_value: float) -> int:
	assert(lambda_value >= 0.0)
	if lambda_value <= 0.0:
		return 0
	var remaining: float = lambda_value
	var result: int = 0
	while remaining > 20.0:
		result += _poisson_small(20.0)
		remaining -= 20.0
	result += _poisson_small(remaining)
	return result

func _poisson_small(lambda_value: float) -> int:
	if lambda_value <= 0.0:
		return 0
	var threshold: float = exp(-lambda_value)
	var product: float = 1.0
	var k: int = 0
	while product > threshold:
		k += 1
		product *= maxf(_rng.randf(), 1e-15)
	return k - 1

func get_state() -> int:
	return _rng.state

func set_state(state_value: int) -> void:
	_rng.state = state_value
