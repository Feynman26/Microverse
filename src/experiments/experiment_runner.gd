extends RefCounted
class_name ExperimentRunner

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")

const SCHEMA_VERSION: int = 1
const MODEL_VERSION: String = "microverse-m8a"

static func create_spec(
	seed: int,
	horizon_ticks: int,
	sample_every_ticks: int,
	environment: Dictionary
) -> Dictionary:
	assert(horizon_ticks > 0 and sample_every_ticks > 0)
	return {
		"schema_version": SCHEMA_VERSION,
		"model_version": MODEL_VERSION,
		"seed": seed,
		"horizon_ticks": horizon_ticks,
		"sample_every_ticks": sample_every_ticks,
		"environment": environment.duplicate(true),
		"stop_on_extinction": true,
		"world_width": 64,
		"world_height": 64,
		"max_cells": 64,
		"mutation_enabled": true,
		"initial_positions": []
	}

static func run(spec: Dictionary) -> Dictionary:
	_validate_spec(spec)
	var config = _create_config(spec)
	var sim = SimulationEngineScript.new(config)
	_seed_initial_population(sim, spec)

	var horizon_ticks: int = int(spec["horizon_ticks"])
	var sample_every: int = int(spec["sample_every_ticks"])
	var stop_on_extinction: bool = bool(spec.get("stop_on_extinction", true))
	var environment: Dictionary = spec.get("environment", EnvironmentScheduleScript.closed())
	var trajectory: Array = []
	var max_population: int = sim.population_size()
	var termination_reason: String = "horizon"
	var realized_ticks: int = 0

	trajectory.append(_sample(sim))
	for tick in range(horizon_ticks):
		EnvironmentScheduleScript.apply(sim, tick, environment)
		sim.step(1)
		realized_ticks += 1
		max_population = maxi(max_population, sim.population_size())
		var extinct: bool = sim.population_size() == 0
		var should_sample: bool = (
			realized_ticks % sample_every == 0
			or realized_ticks == horizon_ticks
			or extinct
		)
		if should_sample:
			trajectory.append(_sample(sim))
		if extinct and stop_on_extinction:
			termination_reason = "extinction"
			break

	return {
		"schema_version": SCHEMA_VERSION,
		"model_version": MODEL_VERSION,
		"seed": int(spec["seed"]),
		"environment": environment.duplicate(true),
		"horizon_ticks": horizon_ticks,
		"realized_ticks": realized_ticks,
		"termination_reason": termination_reason,
		"trajectory": trajectory,
		"max_population": max_population,
		"final_population": sim.population_size(),
		"final_generation": sim.maximum_generation(),
		"final_genotype_count": sim.genotype_count(),
		"division_events": _event_count(sim, "division"),
		"mutation_events": sim.mutation_event_count(),
		"death_causes": _death_causes(sim),
		"final_resources": _field_totals(sim),
		"final_checksum": sim.checksum()
	}

static func run_batch(spec: Dictionary, seeds: Array) -> Dictionary:
	assert(not seeds.is_empty())
	var runs: Array = []
	var by_seed: Dictionary = {}
	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var run_spec: Dictionary = spec.duplicate(true)
		run_spec["seed"] = seed
		var result: Dictionary = run(run_spec)
		runs.append(result)
		by_seed[seed] = result
	return {
		"schema_version": SCHEMA_VERSION,
		"model_version": MODEL_VERSION,
		"runs": runs,
		"by_seed": by_seed
	}

static func _validate_spec(spec: Dictionary) -> void:
	assert(int(spec.get("schema_version", SCHEMA_VERSION)) == SCHEMA_VERSION)
	assert(String(spec.get("model_version", MODEL_VERSION)) == MODEL_VERSION)
	assert(spec.has("seed"))
	assert(int(spec.get("horizon_ticks", 0)) > 0)
	assert(int(spec.get("sample_every_ticks", 0)) > 0)
	assert(int(spec.get("world_width", 64)) > 2)
	assert(int(spec.get("world_height", 64)) > 2)
	assert(int(spec.get("max_cells", 64)) >= 1)

static func _create_config(spec: Dictionary):
	var config = SimConfigScript.new()
	config.seed = int(spec["seed"])
	config.world_width = int(spec.get("world_width", config.world_width))
	config.world_height = int(spec.get("world_height", config.world_height))
	config.max_cells = int(spec.get("max_cells", config.max_cells))
	config.mutation_enabled = bool(spec.get("mutation_enabled", config.mutation_enabled))
	var initial: Dictionary = spec.get("initial_resources", {})
	config.initial_glucose = float(initial.get("glucose", config.initial_glucose))
	config.initial_oxygen = float(initial.get("oxygen", config.initial_oxygen))
	config.initial_nitrogen = float(initial.get("nitrogen", config.initial_nitrogen))
	config.initial_phosphorus = float(initial.get("phosphorus", config.initial_phosphorus))
	config.validate()
	return config

static func _seed_initial_population(sim, spec: Dictionary) -> void:
	var positions: Array = spec.get("initial_positions", [])
	if positions.is_empty():
		sim.seed_ancestor()
		return
	for position_variant in positions:
		assert(position_variant is Vector2)
		sim.seed_ancestor(position_variant)

static func _sample(sim) -> Dictionary:
	return {
		"tick": sim.tick_index,
		"time_min": sim.simulation_time_min,
		"population": sim.population_size(),
		"total_biomass": sim.total_cell_volume(),
		"max_generation": sim.maximum_generation(),
		"genotype_count": sim.genotype_count(),
		"mutation_events": sim.mutation_event_count(),
		"division_events": _event_count(sim, "division"),
		"death_causes": _death_causes(sim),
		"resources": _field_totals(sim),
		"checksum": sim.checksum()
	}

static func _field_totals(sim) -> Dictionary:
	var result: Dictionary = {}
	for field_name in sim.world.field_order:
		result[field_name] = sim.world.get_field(field_name).total_amount()
	return result

static func _death_causes(sim) -> Dictionary:
	var result: Dictionary = {}
	for event in sim.event_log:
		if String(event.get("kind", "")) != "death":
			continue
		var reason: String = String(event.get("reason", "unknown"))
		result[reason] = int(result.get(reason, 0)) + 1
	return result

static func _event_count(sim, kind: String) -> int:
	var result: int = 0
	for event in sim.event_log:
		if String(event.get("kind", "")) == kind:
			result += 1
	return result
