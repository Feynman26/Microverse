extends RefCounted
class_name SimulationDelta

# Compact receipt emitted after a command commits. It describes the affected
# authoritative interval without copying cell, field or event storage.

const SCHEMA_VERSION: int = 1

static func create(
	sequence: int,
	kind: String,
	before: Dictionary,
	after: Dictionary,
	details: Dictionary = {}
) -> Dictionary:
	var result: Dictionary = {
		"schema_version": SCHEMA_VERSION,
		"sequence": sequence,
		"command_kind": kind,
		"from_tick": int(before.get("tick", 0)),
		"to_tick": int(after.get("tick", 0)),
		"from_time_min": float(before.get("simulation_time_min", 0.0)),
		"to_time_min": float(after.get("simulation_time_min", 0.0)),
		"population_before": int(before.get("population", 0)),
		"population_after": int(after.get("population", 0)),
		"event_count_before": int(before.get("event_count", 0)),
		"event_count_after": int(after.get("event_count", 0)),
		"details": details
	}
	validate(result)
	return result

static func validate(delta: Dictionary) -> void:
	assert(int(delta.get("schema_version", -1)) == SCHEMA_VERSION)
	assert(int(delta.get("sequence", 0)) >= 1)
	assert(not String(delta.get("command_kind", "")).is_empty())
	assert(int(delta.get("to_tick", -1)) >= int(delta.get("from_tick", 0)))
	assert(float(delta.get("to_time_min", -1.0)) >= float(delta.get("from_time_min", 0.0)))
	assert(int(delta.get("event_count_after", -1)) >= int(delta.get("event_count_before", 0)))
