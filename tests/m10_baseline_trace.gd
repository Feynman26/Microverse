extends SceneTree

const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const TrajectoryDiagnosticsScript = preload("res://src/experiments/trajectory_diagnostics.gd")

func _init() -> void:
	var sim = SimulationEngineScript.new()
	sim.seed_ancestor()
	var records: Array = [TrajectoryDiagnosticsScript.sample(sim)]
	var max_ticks: int = 4000
	var sample_every: int = 100
	var peak_population: int = sim.population_size()
	var peak_tick: int = 0

	for _unused in range(max_ticks):
		sim.step(1)
		if sim.population_size() > peak_population:
			peak_population = sim.population_size()
			peak_tick = sim.tick_index
		var should_sample: bool = (
			sim.tick_index % sample_every == 0
			or sim.population_size() == 0
		)
		if should_sample:
			records.append(TrajectoryDiagnosticsScript.sample(sim))
		# Once the observed ~16-cell plateau has occurred, collect enough history
		# to see whether division resumes or deaths dominate, without wasting CI.
		if peak_population >= 16 and sim.tick_index >= peak_tick + 800:
			break
		if sim.population_size() == 0:
			break

	for record_variant in records:
		print(TrajectoryDiagnosticsScript.compact_line(record_variant))
	var diagnosis: Dictionary = TrajectoryDiagnosticsScript.diagnose(records)
	print("BASELINE_DIAGNOSIS=" + JSON.stringify(diagnosis))
	print("BASELINE_FINAL_GLOBAL=" + JSON.stringify(records[-1]["global_environment"]))
	print("BASELINE_FINAL_DIVISION=" + JSON.stringify(records[-1]["division"]))
	print("BASELINE_FINAL_MECHANICS=" + JSON.stringify(records[-1]["mechanics"]))
	assert(records.size() >= 2)
	quit()
