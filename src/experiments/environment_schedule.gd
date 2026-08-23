extends RefCounted
class_name EnvironmentSchedule

const MODE_CLOSED: String = "closed"
const MODE_CONSTANT: String = "constant"
const MODE_SQUARE_WAVE: String = "square_wave"
const MODE_SEEDED_PULSES: String = "seeded_pulses"
const MODE_LINEAR_GRADIENT: String = "linear_gradient"
const MODE_PATCH: String = "patch"
const MODE_FORKED: String = "forked"

static func closed() -> Dictionary:
	return {"mode": MODE_CLOSED}

static func constant(fields: Dictionary) -> Dictionary:
	_validate_fields(fields)
	return {"mode": MODE_CONSTANT, "fields": fields.duplicate(true)}

static func square_wave(high_fields: Dictionary, low_fields: Dictionary, phase_ticks: int) -> Dictionary:
	assert(phase_ticks > 0)
	_validate_fields(high_fields)
	_validate_fields(low_fields)
	return {
		"mode": MODE_SQUARE_WAVE,
		"high_fields": high_fields.duplicate(true),
		"low_fields": low_fields.duplicate(true),
		"phase_ticks": phase_ticks
	}

# Deterministic environmental pulses use their own schedule seed and never touch
# the SimulationEngine RNG stream. A pulse is a boundary condition, not a cell
# event; neighboring simulation seeds therefore receive the same schedule when
# the schedule seed is held fixed.
static func seeded_pulses(
	base_fields: Dictionary,
	pulse_fields: Dictionary,
	probability_per_tick: float,
	schedule_seed: int
) -> Dictionary:
	assert(probability_per_tick >= 0.0 and probability_per_tick <= 1.0)
	_validate_fields(base_fields)
	_validate_fields(pulse_fields)
	return {
		"mode": MODE_SEEDED_PULSES,
		"base_fields": base_fields.duplicate(true),
		"pulse_fields": pulse_fields.duplicate(true),
		"probability_per_tick": probability_per_tick,
		"schedule_seed": schedule_seed
	}

# Spatial schedules are ordinary field boundary conditions. They do not query
# cells and therefore cannot encode a hidden fitness or resource-seeking rule.
static func linear_gradient(
	field_name: String,
	low_value: float,
	high_value: float,
	axis: String = "x",
	base_fields: Dictionary = {}
) -> Dictionary:
	assert(not field_name.is_empty())
	assert(low_value >= 0.0 and high_value >= 0.0)
	assert(axis == "x" or axis == "y")
	_validate_fields(base_fields)
	return {
		"mode": MODE_LINEAR_GRADIENT,
		"field_name": field_name,
		"low_value": low_value,
		"high_value": high_value,
		"axis": axis,
		"base_fields": base_fields.duplicate(true)
	}

static func patch(
	field_name: String,
	background_value: float,
	patch_value: float,
	center: Vector2,
	radius_grid: float,
	base_fields: Dictionary = {}
) -> Dictionary:
	assert(not field_name.is_empty())
	assert(background_value >= 0.0 and patch_value >= 0.0)
	assert(radius_grid >= 0.0)
	_validate_fields(base_fields)
	return {
		"mode": MODE_PATCH,
		"field_name": field_name,
		"background_value": background_value,
		"patch_value": patch_value,
		"center": center,
		"radius_grid": radius_grid,
		"base_fields": base_fields.duplicate(true)
	}

# Paired arms can be replayed from one deterministic prefix. The arm schedule is
# evaluated with tick zero at the fork, which makes intervention timing explicit.
static func forked(prefix: Dictionary, fork_tick: int, after_fork: Dictionary) -> Dictionary:
	assert(fork_tick >= 0)
	return {
		"mode": MODE_FORKED,
		"prefix": prefix.duplicate(true),
		"fork_tick": fork_tick,
		"after_fork": after_fork.duplicate(true)
	}

# Apply an external experimental boundary condition before a simulation tick.
# This harness changes only extracellular fields. It does not modify cell state,
# fitness, gene expression, metabolism, division, death, or SimulationEngine RNG.
static func apply(sim, tick: int, spec: Dictionary) -> Dictionary:
	assert(tick >= 0)
	var mode: String = String(spec.get("mode", MODE_CLOSED))
	match mode:
		MODE_LINEAR_GRADIENT:
			return _apply_gradient(sim, spec)
		MODE_PATCH:
			return _apply_patch(sim, spec)
		MODE_FORKED:
			var fork_tick: int = int(spec.get("fork_tick", -1))
			assert(fork_tick >= 0)
			if tick < fork_tick:
				return apply(sim, tick, spec.get("prefix", closed()))
			return apply(sim, tick - fork_tick, spec.get("after_fork", closed()))
		_:
			var targets: Dictionary = targets_for_tick(tick, spec)
			_apply_uniform_targets(sim, targets)
			return targets

static func targets_for_tick(tick: int, spec: Dictionary) -> Dictionary:
	assert(tick >= 0)
	var mode: String = String(spec.get("mode", MODE_CLOSED))
	match mode:
		MODE_CLOSED:
			return {}
		MODE_CONSTANT:
			var fields: Dictionary = spec.get("fields", {})
			_validate_fields(fields)
			return fields.duplicate(true)
		MODE_SQUARE_WAVE:
			var phase_ticks: int = int(spec.get("phase_ticks", 0))
			assert(phase_ticks > 0)
			var high_fields: Dictionary = spec.get("high_fields", {})
			var low_fields: Dictionary = spec.get("low_fields", {})
			_validate_fields(high_fields)
			_validate_fields(low_fields)
			var high_phase: bool = ((tick / phase_ticks) % 2) == 0
			return high_fields.duplicate(true) if high_phase else low_fields.duplicate(true)
		MODE_SEEDED_PULSES:
			var base_fields: Dictionary = spec.get("base_fields", {})
			var pulse_fields: Dictionary = spec.get("pulse_fields", {})
			_validate_fields(base_fields)
			_validate_fields(pulse_fields)
			var probability: float = float(spec.get("probability_per_tick", -1.0))
			assert(probability >= 0.0 and probability <= 1.0)
			var draw: float = _deterministic_unit(int(spec.get("schedule_seed", 0)), tick)
			var result: Dictionary = base_fields.duplicate(true)
			if draw < probability:
				for field_variant in pulse_fields.keys():
					result[field_variant] = pulse_fields[field_variant]
			return result
		MODE_FORKED:
			var fork_tick: int = int(spec.get("fork_tick", -1))
			assert(fork_tick >= 0)
			if tick < fork_tick:
				return targets_for_tick(tick, spec.get("prefix", closed()))
			return targets_for_tick(tick - fork_tick, spec.get("after_fork", closed()))
		MODE_LINEAR_GRADIENT, MODE_PATCH:
			# Spatial schedules do not have one scalar target per field. Return only
			# their uniform companion fields for metadata/introspection.
			var base_fields: Dictionary = spec.get("base_fields", {})
			_validate_fields(base_fields)
			return base_fields.duplicate(true)
		_:
			assert(false, "Unknown environment schedule mode: %s" % mode)
	return {}

static func _apply_gradient(sim, spec: Dictionary) -> Dictionary:
	var base_fields: Dictionary = spec.get("base_fields", {})
	_validate_fields(base_fields)
	_apply_uniform_targets(sim, base_fields)
	var field_name: String = String(spec.get("field_name", ""))
	var low_value: float = float(spec.get("low_value", -1.0))
	var high_value: float = float(spec.get("high_value", -1.0))
	var axis: String = String(spec.get("axis", "x"))
	assert(sim.world.has_field(field_name))
	assert(low_value >= 0.0 and high_value >= 0.0)
	assert(axis == "x" or axis == "y")
	var field = sim.world.get_field(field_name)
	var denominator: float = float(maxi(1, field.width - 1 if axis == "x" else field.height - 1))
	for y in range(field.height):
		for x in range(field.width):
			var coordinate: float = float(x if axis == "x" else y)
			var fraction: float = coordinate / denominator
			field.set_value(x, y, lerpf(low_value, high_value, fraction))
	return {
		"spatial_mode": MODE_LINEAR_GRADIENT,
		"field_name": field_name,
		"low_value": low_value,
		"high_value": high_value,
		"axis": axis,
		"base_fields": base_fields.duplicate(true)
	}

static func _apply_patch(sim, spec: Dictionary) -> Dictionary:
	var base_fields: Dictionary = spec.get("base_fields", {})
	_validate_fields(base_fields)
	_apply_uniform_targets(sim, base_fields)
	var field_name: String = String(spec.get("field_name", ""))
	var background: float = float(spec.get("background_value", -1.0))
	var patch_value: float = float(spec.get("patch_value", -1.0))
	var center: Vector2 = spec.get("center", Vector2.ZERO)
	var radius: float = float(spec.get("radius_grid", -1.0))
	assert(sim.world.has_field(field_name))
	assert(background >= 0.0 and patch_value >= 0.0 and radius >= 0.0)
	var field = sim.world.get_field(field_name)
	for y in range(field.height):
		for x in range(field.width):
			var distance: float = Vector2(float(x), float(y)).distance_to(center)
			field.set_value(x, y, patch_value if distance <= radius else background)
	return {
		"spatial_mode": MODE_PATCH,
		"field_name": field_name,
		"background_value": background,
		"patch_value": patch_value,
		"center": center,
		"radius_grid": radius,
		"base_fields": base_fields.duplicate(true)
	}

static func _apply_uniform_targets(sim, targets: Dictionary) -> void:
	var field_names: Array = targets.keys()
	field_names.sort()
	for field_variant in field_names:
		var field_name := String(field_variant)
		assert(sim.world.has_field(field_name), "Environment schedule references unknown field: %s" % field_name)
		_fill_field(sim.world.get_field(field_name), float(targets[field_variant]))

static func _deterministic_unit(seed: int, tick: int) -> float:
	var value: int = (seed ^ (tick * 1103515245 + 12345)) & 0x7fffffff
	value = (value * 1664525 + 1013904223) & 0x7fffffff
	return float(value) / 2147483648.0

static func _validate_fields(fields: Dictionary) -> void:
	for field_variant in fields.keys():
		assert(not String(field_variant).is_empty())
		assert(float(fields[field_variant]) >= 0.0)

static func _fill_field(field, value: float) -> void:
	assert(value >= 0.0)
	for y in range(field.height):
		for x in range(field.width):
			field.set_value(x, y, value)
