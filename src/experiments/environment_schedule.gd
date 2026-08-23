extends RefCounted
class_name EnvironmentSchedule

const MODE_CLOSED: String = "closed"
const MODE_CONSTANT: String = "constant"
const MODE_SQUARE_WAVE: String = "square_wave"

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

# Apply an external experimental boundary condition before a simulation tick.
# This harness changes only extracellular fields. It does not modify cell state,
# fitness, gene expression, metabolism, division, death, or RNG state.
static func apply(sim, tick: int, spec: Dictionary) -> Dictionary:
	assert(tick >= 0)
	var targets: Dictionary = targets_for_tick(tick, spec)
	var field_names: Array = targets.keys()
	field_names.sort()
	for field_variant in field_names:
		var field_name := String(field_variant)
		assert(sim.world.has_field(field_name), "Environment schedule references unknown field: %s" % field_name)
		_fill_field(sim.world.get_field(field_name), float(targets[field_variant]))
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
		_:
			assert(false, "Unknown environment schedule mode: %s" % mode)
	return {}

static func _validate_fields(fields: Dictionary) -> void:
	for field_variant in fields.keys():
		assert(not String(field_variant).is_empty())
		assert(float(fields[field_variant]) >= 0.0)

static func _fill_field(field, value: float) -> void:
	assert(value >= 0.0)
	for y in range(field.height):
		for x in range(field.width):
			field.set_value(x, y, value)
