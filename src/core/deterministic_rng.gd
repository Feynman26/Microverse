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

func get_state() -> int:
	return _rng.state

func set_state(state_value: int) -> void:
	_rng.state = state_value
