extends "res://src/runtime/i_simulation_backend.gd"
class_name LegacySimulationBackend

const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const SimulationCommandScript = preload("res://src/runtime/simulation_command.gd")
const SimulationDeltaScript = preload("res://src/runtime/simulation_delta.gd")
const VisualSnapshotScript = preload("res://src/observation/visual_snapshot.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")

const BACKEND_ID: String = "m10-legacy-adapter"
const MODEL_VERSION: String = "microverse-m10"

var _simulation = null
var _command_sequence: int = 0
var _command_counts: Dictionary = {}

func _init(config = null) -> void:
	if config != null:
		config.validate()
		_simulation = SimulationEngineScript.new(config)

func execute(command: Dictionary) -> Dictionary:
	_require_initialized()
	SimulationCommandScript.validate(command)
	var kind: String = String(command["kind"])
	var payload: Dictionary = command["payload"]
	var before: Dictionary = state_summary()
	var details: Dictionary = {}

	match kind:
		SimulationCommandScript.SEED_ANCESTOR:
			var cell = _simulation.seed_ancestor(payload.get("position", Vector2(-1.0, -1.0)))
			details = {"cell_id": int(cell.id)}
		SimulationCommandScript.ADVANCE_TICKS:
			_simulation.step(int(payload.get("tick_count", 1)))
		SimulationCommandScript.APPLY_ENVIRONMENT:
			details = EnvironmentScheduleScript.apply(
				_simulation,
				int(payload["tick"]),
				payload["schedule"]
			)
		SimulationCommandScript.SERIAL_TRANSFER:
			details = _apply_serial_transfer(
				int(payload["seed"]),
				int(payload["tick"]),
				int(payload["survivors"])
			)

	_command_sequence += 1
	_command_counts[kind] = int(_command_counts.get(kind, 0)) + 1
	return SimulationDeltaScript.create(
		_command_sequence,
		kind,
		before,
		state_summary(),
		details
	)

func capture_visual_snapshot() -> Dictionary:
	_require_initialized()
	var snapshot: Dictionary = VisualSnapshotScript.capture(_simulation)
	snapshot["backend"] = metadata()
	return snapshot

func state_summary() -> Dictionary:
	_require_initialized()
	return {
		"tick": int(_simulation.tick_index),
		"simulation_time_min": float(_simulation.simulation_time_min),
		"population": int(_simulation.population_size()),
		"max_cells": int(_simulation.config.max_cells),
		"event_count": int(_simulation.event_log.size())
	}

func telemetry_report() -> Dictionary:
	_require_initialized()
	return {
		"schema_version": 1,
		"backend": metadata(),
		"commands": _command_counts.duplicate(),
		"simulation": _simulation.performance_profiler.structured_report()
	}

func metadata() -> Dictionary:
	return {
		"interface_version": INTERFACE_VERSION,
		"backend_id": BACKEND_ID,
		"model_version": MODEL_VERSION,
		"command_schema_version": SimulationCommandScript.SCHEMA_VERSION,
		"delta_schema_version": SimulationDeltaScript.SCHEMA_VERSION,
		"visual_snapshot_schema_version": VisualSnapshotScript.SCHEMA_VERSION
	}

func legacy_inspection_state():
	_require_initialized()
	return _simulation

func _apply_serial_transfer(seed: int, tick: int, survivors: int) -> Dictionary:
	var before: int = _simulation.population_size()
	if before <= survivors:
		return {
			"tick": tick,
			"kind": SimulationCommandScript.SERIAL_TRANSFER,
			"before": before,
			"after": before,
			"removed_ids": []
		}

	var ranked: Array = []
	for cell in _simulation.cells:
		ranked.append({"rank": _intervention_rank(seed, tick, int(cell.id)), "cell": cell})
	ranked.sort_custom(func(a: Dictionary, b: Dictionary) -> bool:
		if int(a["rank"]) == int(b["rank"]):
			return int(a["cell"].id) < int(b["cell"].id)
		return int(a["rank"]) < int(b["rank"])
	)
	var kept: Array = []
	var removed_ids: Array = []
	for index in range(ranked.size()):
		if index < survivors:
			kept.append(ranked[index]["cell"])
		else:
			removed_ids.append(int(ranked[index]["cell"].id))
	_simulation.cells = kept
	return {
		"tick": tick,
		"kind": SimulationCommandScript.SERIAL_TRANSFER,
		"before": before,
		"after": kept.size(),
		"removed_ids": removed_ids
	}

func _intervention_rank(seed: int, tick: int, cell_id: int) -> int:
	var value: int = (seed ^ (tick * 1103515245) ^ (cell_id * 2654435761)) & 0x7fffffff
	value = (value * 1664525 + 1013904223) & 0x7fffffff
	return value

func _require_initialized() -> void:
	assert(_simulation != null, "Legacy simulation backend requires a validated config")
