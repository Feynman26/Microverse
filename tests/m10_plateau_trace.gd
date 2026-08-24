extends SceneTree

const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const TrajectoryDiagnosticsScript = preload("res://src/experiments/trajectory_diagnostics.gd")

# One-off M10 characterization of the exact default visual configuration. This
# is evidence, not an acceptance test: it never changes configuration or cell
# physiology and intentionally allows extinction/plateau to occur naturally.
func _init() -> void:
	var sim = SimulationEngineScript.new()
	sim.seed_ancestor()
	var records: Array = [TrajectoryDiagnosticsScript.sample(sim)]
	var max_ticks: int = 9000
	var sample_every: int = 100
	var peak_population: int = sim.population_size()
	var peak_tick: int = 0
	var first_sixteen_tick: int = -1
	var divisions_at_sixteen: int = -1
	var deaths_at_sixteen: int = -1

	for _unused in range(max_ticks):
		sim.step(1)
		if sim.population_size() > peak_population:
			peak_population = sim.population_size()
			peak_tick = sim.tick_index
		if first_sixteen_tick < 0 and sim.population_size() >= 16:
			first_sixteen_tick = sim.tick_index
			var snap: Dictionary = TrajectoryDiagnosticsScript.sample(sim)
			divisions_at_sixteen = int(snap["division_events"])
			deaths_at_sixteen = int(snap["death_events"])
			records.append(snap)
		var should_sample: bool = (
			sim.tick_index % sample_every == 0
			or sim.population_size() == 0
		)
		if should_sample:
			records.append(TrajectoryDiagnosticsScript.sample(sim))
		# Observe 120 virtual minutes after first reaching 16. That is long enough
		# to distinguish a hard division ceiling from continued turnover/resumed
		# births while bounding this characterization's CPU cost.
		if first_sixteen_tick >= 0 and sim.tick_index >= first_sixteen_tick + 1200:
			break
		if sim.population_size() == 0:
			break

	# Print the complete sampled history so the plateau can be inspected rather
	# than reduced to one endpoint. Duplicate tick samples are harmless evidence.
	for record_variant in records:
		print(TrajectoryDiagnosticsScript.compact_line(record_variant))
	var diagnosis: Dictionary = TrajectoryDiagnosticsScript.diagnose(records)
	var final_record: Dictionary = records[-1]
	print("PLATEAU_RESULT=" + JSON.stringify({
		"seed": int(sim.config.seed),
		"max_ticks": max_ticks,
		"realized_ticks": int(sim.tick_index),
		"peak_population": peak_population,
		"peak_tick": peak_tick,
		"first_sixteen_tick": first_sixteen_tick,
		"divisions_at_sixteen": divisions_at_sixteen,
		"deaths_at_sixteen": deaths_at_sixteen,
		"final_population": int(sim.population_size()),
		"final_generation": int(sim.maximum_generation()),
		"final_mutations": int(sim.mutation_event_count()),
		"diagnosis": diagnosis,
		"final_division": final_record["division"],
		"final_replication": final_record["replication"],
		"final_damage": final_record["damage"],
		"final_debt": final_record["energy_debt"],
		"final_local_environment": final_record["local_environment"],
		"final_global_environment": final_record["global_environment"],
		"final_mechanics": final_record["mechanics"]
	}))
	quit()
