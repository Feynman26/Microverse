extends RefCounted
class_name SimulationCommand

# Serializable, immutable-by-contract commands entering the authoritative
# backend. Runtime clock controls remain outside this contract because they do
# not themselves change scientific state.

const SCHEMA_VERSION: int = 1

const SEED_ANCESTOR: String = "seed_ancestor"
const ADVANCE_TICKS: String = "advance_ticks"
const APPLY_ENVIRONMENT: String = "apply_environment"
const SERIAL_TRANSFER: String = "serial_transfer"

static func seed_ancestor(position: Vector2 = Vector2(-1.0, -1.0)) -> Dictionary:
	return _create(SEED_ANCESTOR, {"position": position})

static func advance_ticks(tick_count: int = 1) -> Dictionary:
	assert(tick_count >= 0)
	return _create(ADVANCE_TICKS, {"tick_count": tick_count})

static func apply_environment(tick: int, schedule: Dictionary) -> Dictionary:
	assert(tick >= 0)
	return _create(APPLY_ENVIRONMENT, {"tick": tick, "schedule": schedule})

static func serial_transfer(seed: int, tick: int, survivors: int) -> Dictionary:
	assert(tick >= 0 and survivors >= 1)
	return _create(SERIAL_TRANSFER, {
		"seed": seed,
		"tick": tick,
		"survivors": survivors
	})

static func validate(command: Dictionary) -> void:
	assert(int(command.get("schema_version", -1)) == SCHEMA_VERSION)
	var kind: String = String(command.get("kind", ""))
	assert(kind in [SEED_ANCESTOR, ADVANCE_TICKS, APPLY_ENVIRONMENT, SERIAL_TRANSFER])
	var payload: Dictionary = command.get("payload", {})
	match kind:
		SEED_ANCESTOR:
			assert(payload.get("position", Vector2.ZERO) is Vector2)
		ADVANCE_TICKS:
			assert(int(payload.get("tick_count", -1)) >= 0)
		APPLY_ENVIRONMENT:
			assert(int(payload.get("tick", -1)) >= 0)
			assert(payload.get("schedule", null) is Dictionary)
		SERIAL_TRANSFER:
			assert(int(payload.get("tick", -1)) >= 0)
			assert(int(payload.get("survivors", 0)) >= 1)

static func _create(kind: String, payload: Dictionary) -> Dictionary:
	return {
		"schema_version": SCHEMA_VERSION,
		"kind": kind,
		"payload": payload
	}
