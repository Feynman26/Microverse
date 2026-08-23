extends RefCounted
class_name TrajectoryDiagnostics

# Post-hoc/observational diagnostics only. This module reads authoritative state
# and never changes physiology, RNG, resources, mutation, division or survival.

static func sample(sim) -> Dictionary:
	var population: int = sim.population_size()
	var volumes: Array = []
	var atp_values: Array = []
	var damage_values: Array = []
	var debt_values: Array = []
	var ros_values: Array = []
	var ages: Array = []
	var local_glucose: Array = []
	var local_oxygen: Array = []
	var local_nitrogen: Array = []
	var local_phosphorus: Array = []
	var division_ready: int = 0
	var volume_ready: int = 0
	var atp_ready: int = 0
	var volume_ready_atp_blocked: int = 0
	var atp_ready_volume_blocked: int = 0
	var oversized_ids: Array = []

	for cell in sim.cells:
		if not cell.alive:
			continue
		var volume: float = float(cell.volume)
		var atp: float = float(cell.pool("ATP"))
		var has_volume: bool = volume >= float(sim.config.division_volume)
		var has_atp: bool = atp >= float(sim.config.division_atp_cost)
		volumes.append(volume)
		atp_values.append(atp)
		damage_values.append(float(cell.damage))
		debt_values.append(float(cell.energy_debt))
		ros_values.append(float(cell.pool("ROS")))
		ages.append(float(sim.tick_index - int(cell.birth_tick)))
		local_glucose.append(float(sim.world.sample("glucose", cell.position)))
		local_oxygen.append(float(sim.world.sample("oxygen", cell.position)))
		local_nitrogen.append(float(sim.world.sample("nitrogen", cell.position)))
		local_phosphorus.append(float(sim.world.sample("phosphorus", cell.position)))
		if has_volume: volume_ready += 1
		if has_atp: atp_ready += 1
		if has_volume and has_atp:
			division_ready += 1
		elif has_volume:
			volume_ready_atp_blocked += 1
			if oversized_ids.size() < 12:
				oversized_ids.append(int(cell.id))
		elif has_atp:
			atp_ready_volume_blocked += 1

	return {
		"tick": int(sim.tick_index),
		"time_min": float(sim.simulation_time_min),
		"population": population,
		"max_generation": int(sim.maximum_generation()),
		"genotype_count": int(sim.genotype_count()),
		"mutation_events": int(sim.mutation_event_count()),
		"division_events": _event_count(sim, "division"),
		"death_events": _event_count(sim, "death"),
		"death_causes": _death_causes(sim),
		"total_biomass": float(sim.total_cell_volume()),
		"volume": _stats(volumes),
		"atp": _stats(atp_values),
		"damage": _stats(damage_values),
		"energy_debt": _stats(debt_values),
		"intracellular_ros": _stats(ros_values),
		"age_ticks": _stats(ages),
		"division": {
			"division_volume": float(sim.config.division_volume),
			"division_atp_cost": float(sim.config.division_atp_cost),
			"ready": division_ready,
			"volume_ready": volume_ready,
			"atp_ready": atp_ready,
			"volume_ready_atp_blocked": volume_ready_atp_blocked,
			"atp_ready_volume_blocked": atp_ready_volume_blocked,
			"oversized_blocked_ids": oversized_ids
		},
		"local_environment": {
			"glucose": _stats(local_glucose),
			"oxygen": _stats(local_oxygen),
			"nitrogen": _stats(local_nitrogen),
			"phosphorus": _stats(local_phosphorus)
		},
		"global_environment": _field_totals(sim),
		"mechanics": sim.last_mechanics_summary.duplicate(true)
	}

static func trace(sim, horizon_ticks: int, sample_every_ticks: int) -> Array:
	assert(horizon_ticks >= 0 and sample_every_ticks > 0)
	var result: Array = [sample(sim)]
	for _tick in range(horizon_ticks):
		sim.step(1)
		if sim.tick_index % sample_every_ticks == 0 or sim.population_size() == 0:
			result.append(sample(sim))
		if sim.population_size() == 0:
			break
	return result

static func compact_line(record: Dictionary) -> String:
	var volume: Dictionary = record["volume"]
	var atp: Dictionary = record["atp"]
	var damage: Dictionary = record["damage"]
	var debt: Dictionary = record["energy_debt"]
	var division: Dictionary = record["division"]
	var local: Dictionary = record["local_environment"]
	return (
		"TRACE tick=%d t=%.1f pop=%d gen=%d div=%d mut=%d deaths=%d "
		+ "BIO=%.3f vol_mean=%.3f vol_max=%.3f ATP_mean=%.3f ATP_min=%.3f ATP_max=%.3f "
		+ "ready=%d vol_ready=%d vol_ATP_block=%d ATP_vol_block=%d "
		+ "damage_mean=%.3f damage_max=%.3f debt_mean=%.3f debt_max=%.3f "
		+ "local_G=%.3f local_O2=%.3f local_N=%.3f local_P=%.3f"
	) % [
		int(record["tick"]), float(record["time_min"]), int(record["population"]), int(record["max_generation"]),
		int(record["division_events"]), int(record["mutation_events"]), int(record["death_events"]),
		float(record["total_biomass"]), float(volume["mean"]), float(volume["max"]),
		float(atp["mean"]), float(atp["min"]), float(atp["max"]),
		int(division["ready"]), int(division["volume_ready"]), int(division["volume_ready_atp_blocked"]),
		int(division["atp_ready_volume_blocked"]), float(damage["mean"]), float(damage["max"]),
		float(debt["mean"]), float(debt["max"]),
		float(local["glucose"]["mean"]), float(local["oxygen"]["mean"]),
		float(local["nitrogen"]["mean"]), float(local["phosphorus"]["mean"])
	]

static func diagnose(trace_records: Array) -> Dictionary:
	assert(not trace_records.is_empty())
	var max_population: int = 0
	var first_sustained_atp_block_tick: int = -1
	var max_volume: float = 0.0
	var max_blocked: int = 0
	var first_population_peak_tick: int = 0
	for record_variant in trace_records:
		var record: Dictionary = record_variant
		var population: int = int(record["population"])
		if population > max_population:
			max_population = population
			first_population_peak_tick = int(record["tick"])
		var blocked: int = int(record["division"]["volume_ready_atp_blocked"])
		max_blocked = maxi(max_blocked, blocked)
		max_volume = maxf(max_volume, float(record["volume"]["max"]))
		if first_sustained_atp_block_tick < 0 and blocked > 0:
			first_sustained_atp_block_tick = int(record["tick"])
	var final_record: Dictionary = trace_records[-1]
	return {
		"max_population": max_population,
		"population_peak_tick": first_population_peak_tick,
		"max_volume": max_volume,
		"max_volume_ready_atp_blocked": max_blocked,
		"first_volume_ready_atp_block_tick": first_sustained_atp_block_tick,
		"final_population": int(final_record["population"]),
		"final_death_causes": final_record["death_causes"].duplicate(true),
		"final_division_state": final_record["division"].duplicate(true)
	}

static func _stats(values: Array) -> Dictionary:
	if values.is_empty():
		return {"min": 0.0, "mean": 0.0, "max": 0.0, "total": 0.0}
	var minimum: float = INF
	var maximum: float = -INF
	var total: float = 0.0
	for value_variant in values:
		var value: float = float(value_variant)
		minimum = minf(minimum, value)
		maximum = maxf(maximum, value)
		total += value
	return {"min": minimum, "mean": total / float(values.size()), "max": maximum, "total": total}

static func _field_totals(sim) -> Dictionary:
	var result: Dictionary = {}
	for field_name in sim.world.field_order:
		result[field_name] = float(sim.world.get_field(field_name).total_amount())
	return result

static func _event_count(sim, kind: String) -> int:
	var result: int = 0
	for event_variant in sim.event_log:
		if String(event_variant.get("kind", "")) == kind:
			result += 1
	return result

static func _death_causes(sim) -> Dictionary:
	var result: Dictionary = {}
	for event_variant in sim.event_log:
		if String(event_variant.get("kind", "")) != "death":
			continue
		var reason: String = String(event_variant.get("reason", "unknown"))
		result[reason] = int(result.get(reason, 0)) + 1
	return result
