extends RefCounted
class_name ExperimentRunner

const SimConfigScript = preload("res://src/core/sim_config.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")
const MutationDynamicsAnalyticsScript = preload("res://src/experiments/mutation_dynamics_analytics.gd")
const LegacySimulationBackendScript = preload("res://src/runtime/legacy_simulation_backend.gd")
const SimulationCommandScript = preload("res://src/runtime/simulation_command.gd")

const SCHEMA_VERSION: int = 2
const MODEL_VERSION: String = "microverse-m10"
const DEFAULT_MAX_CELLS: int = 256

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
		"interventions": [],
		"stop_on_extinction": true,
		"stop_population_at_least": 0,
		"world_width": 64,
		"world_height": 64,
		# This is a computational guard, never an ecological carrying capacity.
		# A run terminates explicitly when it reaches the guard so suppressed
		# divisions cannot be misread as stationary-phase biology.
		"max_cells": DEFAULT_MAX_CELLS,
		"mutation_enabled": true,
		"initial_positions": []
	}

static func run(spec: Dictionary) -> Dictionary:
	_validate_spec(spec)
	var config = _create_config(spec)
	var backend = LegacySimulationBackendScript.new(config)
	_seed_initial_population(backend, spec)
	# P2 keeps existing analytics on an explicitly read-only adapter view. Every
	# authoritative mutation below enters through a versioned backend command.
	var sim = backend.legacy_inspection_state()

	var horizon_ticks: int = int(spec["horizon_ticks"])
	var sample_every: int = int(spec["sample_every_ticks"])
	var stop_on_extinction: bool = bool(spec.get("stop_on_extinction", true))
	var stop_population_at_least: int = int(spec.get("stop_population_at_least", 0))
	var environment: Dictionary = spec.get("environment", EnvironmentScheduleScript.closed())
	var interventions: Array = spec.get("interventions", [])
	var intervention_log: Array = []
	var trajectory: Array = []
	var max_population: int = sim.population_size()
	var termination_reason: String = "horizon"
	var realized_ticks: int = 0

	trajectory.append(_sample(sim))
	for tick in range(horizon_ticks):
		_apply_interventions(backend, tick, interventions, int(spec["seed"]), intervention_log)
		backend.execute(SimulationCommandScript.apply_environment(tick, environment))
		backend.execute(SimulationCommandScript.advance_ticks(1))
		realized_ticks += 1
		max_population = maxi(max_population, sim.population_size())
		var extinct: bool = sim.population_size() == 0
		var population_threshold_reached: bool = (
			stop_population_at_least > 0
			and sim.population_size() >= stop_population_at_least
		)
		var computational_limit_reached: bool = sim.population_size() >= int(config.max_cells)
		var should_sample: bool = (
			realized_ticks % sample_every == 0
			or realized_ticks == horizon_ticks
			or extinct
			or population_threshold_reached
			or computational_limit_reached
		)
		if should_sample:
			trajectory.append(_sample(sim))
		if extinct and stop_on_extinction:
			termination_reason = "extinction"
			break
		if population_threshold_reached:
			termination_reason = "population_threshold"
			break
		if computational_limit_reached:
			termination_reason = "computational_population_limit"
			break

	return {
		"schema_version": SCHEMA_VERSION,
		"model_version": MODEL_VERSION,
		"seed": int(spec["seed"]),
		"environment": environment.duplicate(true),
		"interventions": interventions.duplicate(true),
		"intervention_log": intervention_log,
		"horizon_ticks": horizon_ticks,
		"realized_ticks": realized_ticks,
		"termination_reason": termination_reason,
		"computational_population_limit": int(config.max_cells),
		"computational_limit_reached": termination_reason == "computational_population_limit",
		"trajectory": trajectory,
		"max_population": max_population,
		"final_population": sim.population_size(),
		"final_generation": sim.maximum_generation(),
		"final_genotype_count": sim.genotype_count(),
		"final_genotype_frequencies": _genotype_frequencies(sim),
		"final_generation_distribution": _generation_distribution(sim),
		"division_events": _event_count(sim, "division"),
		"mutation_events": sim.mutation_event_count(),
		"mutation_event_summary": MutationDynamicsAnalyticsScript.summarize_event_log(sim.event_log),
		"death_causes": _death_causes(sim),
		"final_resources": _field_totals(sim),
		"final_cell_diagnostics": _cell_diagnostics(sim),
		"final_mutation_dynamics": MutationDynamicsAnalyticsScript.sample_population(sim),
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

# Paired forks use deterministic prefix replay. Every arm has the same seed,
# initial state, prefix schedule and interventions through fork_tick; prefix
# checksums prove identical authoritative state before environmental divergence.
static func run_paired_fork(
	base_spec: Dictionary,
	fork_tick: int,
	arm_environments: Dictionary
) -> Dictionary:
	_validate_spec(base_spec)
	assert(fork_tick > 0 and fork_tick < int(base_spec["horizon_ticks"]))
	assert(arm_environments.size() >= 2)
	var prefix_spec: Dictionary = base_spec.duplicate(true)
	prefix_spec["horizon_ticks"] = fork_tick
	prefix_spec["sample_every_ticks"] = fork_tick
	prefix_spec["stop_population_at_least"] = 0
	var canonical_prefix: Dictionary = run(prefix_spec)
	var fork_checksum: float = float(canonical_prefix["final_checksum"])

	var arms: Dictionary = {}
	var prefix_checksums: Dictionary = {}
	var arm_names: Array = arm_environments.keys()
	arm_names.sort()
	for name_variant in arm_names:
		var name: String = String(name_variant)
		var replay_prefix: Dictionary = run(prefix_spec)
		prefix_checksums[name] = float(replay_prefix["final_checksum"])
		assert(absf(float(replay_prefix["final_checksum"]) - fork_checksum) <= 1e-12)
		var arm_spec: Dictionary = base_spec.duplicate(true)
		arm_spec["environment"] = EnvironmentScheduleScript.forked(
			base_spec.get("environment", EnvironmentScheduleScript.closed()),
			fork_tick,
			arm_environments[name_variant]
		)
		arms[name] = run(arm_spec)
	return {
		"schema_version": SCHEMA_VERSION,
		"model_version": MODEL_VERSION,
		"seed": int(base_spec["seed"]),
		"fork_tick": fork_tick,
		"fork_checksum": fork_checksum,
		"prefix_checksums": prefix_checksums,
		"arms": arms
	}

static func _validate_spec(spec: Dictionary) -> void:
	assert(int(spec.get("schema_version", SCHEMA_VERSION)) == SCHEMA_VERSION)
	assert(String(spec.get("model_version", MODEL_VERSION)) == MODEL_VERSION)
	assert(spec.has("seed"))
	assert(int(spec.get("horizon_ticks", 0)) > 0)
	assert(int(spec.get("sample_every_ticks", 0)) > 0)
	assert(int(spec.get("stop_population_at_least", 0)) >= 0)
	assert(int(spec.get("world_width", 64)) > 2)
	assert(int(spec.get("world_height", 64)) > 2)
	assert(int(spec.get("max_cells", DEFAULT_MAX_CELLS)) >= 1)
	if int(spec.get("stop_population_at_least", 0)) > 0:
		assert(int(spec.get("stop_population_at_least", 0)) <= int(spec.get("max_cells", DEFAULT_MAX_CELLS)))
	for intervention_variant in spec.get("interventions", []):
		var intervention: Dictionary = intervention_variant
		assert(String(intervention.get("kind", "")) == "serial_transfer")
		assert(int(intervention.get("tick", -1)) >= 0)
		assert(int(intervention.get("survivors", 0)) >= 1)

static func _create_config(spec: Dictionary):
	var config = SimConfigScript.new()
	config.seed = int(spec["seed"])
	config.world_width = int(spec.get("world_width", config.world_width))
	config.world_height = int(spec.get("world_height", config.world_height))
	config.max_cells = int(spec.get("max_cells", DEFAULT_MAX_CELLS))
	config.mutation_enabled = bool(spec.get("mutation_enabled", config.mutation_enabled))
	var initial: Dictionary = spec.get("initial_resources", {})
	config.initial_glucose = float(initial.get("glucose", config.initial_glucose))
	config.initial_oxygen = float(initial.get("oxygen", config.initial_oxygen))
	config.initial_nitrogen = float(initial.get("nitrogen", config.initial_nitrogen))
	config.initial_phosphorus = float(initial.get("phosphorus", config.initial_phosphorus))
	config.validate()
	return config

static func _seed_initial_population(backend, spec: Dictionary) -> void:
	var positions: Array = spec.get("initial_positions", [])
	if positions.is_empty():
		backend.execute(SimulationCommandScript.seed_ancestor())
		return
	for position_variant in positions:
		assert(position_variant is Vector2)
		backend.execute(SimulationCommandScript.seed_ancestor(position_variant))

static func _sample(sim) -> Dictionary:
	return {
		"tick": sim.tick_index,
		"time_min": sim.simulation_time_min,
		"population": sim.population_size(),
		"total_biomass": sim.total_cell_volume(),
		"max_generation": sim.maximum_generation(),
		"generation_distribution": _generation_distribution(sim),
		"genotype_count": sim.genotype_count(),
		"genotype_frequencies": _genotype_frequencies(sim),
		"mutation_events": sim.mutation_event_count(),
		"division_events": _event_count(sim, "division"),
		"death_causes": _death_causes(sim),
		"resources": _field_totals(sim),
		"cell_diagnostics": _cell_diagnostics(sim),
		"mutation_dynamics": MutationDynamicsAnalyticsScript.sample_population(sim),
		"checksum": sim.checksum()
	}

static func _apply_interventions(backend, tick: int, interventions: Array, seed: int, log: Array) -> void:
	for intervention_variant in interventions:
		var intervention: Dictionary = intervention_variant
		if int(intervention.get("tick", -1)) != tick:
			continue
		match String(intervention.get("kind", "")):
			"serial_transfer":
				var survivors: int = int(intervention.get("survivors", 1))
				var delta: Dictionary = backend.execute(
					SimulationCommandScript.serial_transfer(seed, tick, survivors)
				)
				log.append(delta["details"].duplicate(true))

static func _genotype_frequencies(sim) -> Dictionary:
	var counts: Dictionary = {}
	for cell in sim.cells:
		var key: String = str(cell.genome.fingerprint())
		counts[key] = int(counts.get(key, 0)) + 1
	return counts

static func _generation_distribution(sim) -> Dictionary:
	var counts: Dictionary = {}
	for cell in sim.cells:
		var key: String = str(int(cell.generation))
		counts[key] = int(counts.get(key, 0)) + 1
	return counts

static func _cell_diagnostics(sim) -> Dictionary:
	var population: int = sim.population_size()
	if population == 0:
		return {
			"mean_damage": 0.0,
			"max_damage": 0.0,
			"mean_energy_debt": 0.0,
			"max_energy_debt": 0.0,
			"total_intracellular_ros": 0.0
		}
	var total_damage: float = 0.0
	var max_damage: float = 0.0
	var total_energy_debt: float = 0.0
	var max_energy_debt: float = 0.0
	var total_ros: float = 0.0
	for cell in sim.cells:
		total_damage += float(cell.damage)
		max_damage = maxf(max_damage, float(cell.damage))
		total_energy_debt += float(cell.energy_debt)
		max_energy_debt = maxf(max_energy_debt, float(cell.energy_debt))
		total_ros += float(cell.pool("ROS"))
	return {
		"mean_damage": total_damage / float(population),
		"max_damage": max_damage,
		"mean_energy_debt": total_energy_debt / float(population),
		"max_energy_debt": max_energy_debt,
		"total_intracellular_ros": total_ros
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
