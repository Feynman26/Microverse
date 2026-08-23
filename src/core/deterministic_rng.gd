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

# Exact inverse-CDF Poisson sampler for the small event intensities used by M5.
# Crucially, every poisson() call consumes exactly one uniform RNG draw,
# including lambda=0. The former Knuth loop consumed a variable number of
# uniforms, so a regulatory mutation that changed one transcription lambda also
# phase-shifted the random stream seen by unrelated downstream genes. That was a
# reproducible but nonlocal stochastic artifact, not molecular biology.
#
# M5 event lambdas are far below the numerical-risk regime for inverse CDF. A
# conservative upper guard makes that assumption explicit instead of silently
# switching to a second approximation with different draw cardinality.
func poisson(lambda_value: float) -> int:
	assert(lambda_value >= 0.0)
	assert(lambda_value <= 100.0, "Exact one-draw Poisson currently supports lambda <= 100")
	var u: float = clampf(_rng.randf(), 0.0, 1.0 - 1e-15)
	if lambda_value <= 0.0:
		return 0
	var probability: float = exp(-lambda_value)
	var cumulative: float = probability
	var k: int = 0
	while u > cumulative:
		k += 1
		probability *= lambda_value / float(k)
		cumulative += probability
		if probability <= 1e-16 and cumulative >= 1.0 - 1e-14:
			break
	return k

func get_state() -> int:
	return _rng.state

func set_state(state_value: int) -> void:
	_rng.state = state_value
