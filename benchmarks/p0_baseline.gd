extends SceneTree

const SCHEMA_VERSION: int = 1
const MODEL_VERSION: String = "microverse-m10"
const DEFAULT_POPULATIONS: Array[int] = [1, 16, 64, 256, 1000]
const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var _warmup_ticks: int = 10
var _measured_ticks: int = 30
var _populations: Array[int] = DEFAULT_POPULATIONS.duplicate()
var _output_path: String = ""

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_parse_args(OS.get_cmdline_user_args())
	var manifest: Dictionary = _machine_manifest()
	var results: Array = []
	for population in _populations:
		print("P0: benchmarking %d initial cells (%d warmup + %d measured ticks)" % [population, _warmup_ticks, _measured_ticks])
		results.append(_run_scenario(population))
	var report: Dictionary = {
		"schema_version": SCHEMA_VERSION,
		"model_version": MODEL_VERSION,
		"profiler_version": "p0-v1",
		"purpose": "observational_performance_baseline",
		"machine": manifest,
		"parameters": {
			"warmup_ticks": _warmup_ticks,
			"measured_ticks": _measured_ticks,
			"populations": _populations
		},
		"results": results
	}
	var json: String = JSON.stringify(report, "  ")
	print("P0_RESULT_JSON_BEGIN")
	print(json)
	print("P0_RESULT_JSON_END")
	if not _output_path.is_empty():
		_write_report(_output_path, json)
	quit(0)

func _run_scenario(initial_population: int) -> Dictionary:
	var config = SimConfigScript.new()
	config.seed = 903000 + initial_population
	config.max_cells = maxi(initial_population * 2, initial_population + 2)
	config.mutation_enabled = false
	config.performance_profiling_enabled = true
	config.validate()
	var sim = SimulationEngineScript.new(config)
	_seed_grid(sim, initial_population)
	for _tick in range(_warmup_ticks):
		sim.step(1)
	sim.performance_profiler.reset()
	var memory_before: Dictionary = _memory_snapshot()
	var started_usec: int = Time.get_ticks_usec()
	for _tick in range(_measured_ticks):
		sim.step(1)
	var elapsed_usec: int = maxi(1, Time.get_ticks_usec() - started_usec)
	var memory_after: Dictionary = _memory_snapshot()
	return {
		"scenario": "fixed_64x64_m10_cells_%d" % initial_population,
		"initial_population": initial_population,
		"final_population": sim.population_size(),
		"world_width": config.world_width,
		"world_height": config.world_height,
		"registered_fields": sim.world.field_order.size(),
		"dynamic_protein_fields": sim.world.protein_fields.size(),
		"measured_ticks": _measured_ticks,
		"elapsed_usec": elapsed_usec,
		"ticks_per_second": float(_measured_ticks) * 1000000.0 / float(elapsed_usec),
		"mean_usec_per_tick": float(elapsed_usec) / float(_measured_ticks),
		"phase_profile": sim.performance_profiler.report(),
		"memory_before": memory_before,
		"memory_after": memory_after,
		"final_checksum": sim.checksum(),
		"event_count": sim.event_log.size(),
		"mechanics": sim.last_mechanics_summary.duplicate(true)
	}

func _seed_grid(sim, population: int) -> void:
	var spacing: float = 1.05
	var margin: float = 0.55
	var columns: int = int(floor((float(sim.world.width - 1) - 2.0 * margin) / spacing)) + 1
	assert(columns > 0)
	for i in range(population):
		var x: float = margin + float(i % columns) * spacing
		var y: float = margin + float(i / columns) * spacing
		assert(y <= float(sim.world.height - 1) - margin, "P0 population does not fit benchmark chamber")
		sim.seed_ancestor(Vector2(x, y))

func _machine_manifest() -> Dictionary:
	var version: Dictionary = Engine.get_version_info()
	return {
		"godot_version": String(version.get("string", "unknown")),
		"os_name": OS.get_name(),
		"os_distribution": OS.get_distribution_name(),
		"os_version": OS.get_version(),
		"processor_name": OS.get_processor_name(),
		"processor_count": OS.get_processor_count(),
		"display_server": DisplayServer.get_name(),
		"command_line": OS.get_cmdline_args()
	}

func _memory_snapshot() -> Dictionary:
	return {
		"static_bytes": int(Performance.get_monitor(Performance.MEMORY_STATIC)),
		"static_peak_bytes": int(Performance.get_monitor(Performance.MEMORY_STATIC_MAX)),
		"object_count": int(Performance.get_monitor(Performance.OBJECT_COUNT)),
		"resource_count": int(Performance.get_monitor(Performance.OBJECT_RESOURCE_COUNT)),
		"node_count": int(Performance.get_monitor(Performance.OBJECT_NODE_COUNT)),
		"orphan_node_count": int(Performance.get_monitor(Performance.OBJECT_ORPHAN_NODE_COUNT))
	}

func _parse_args(args: PackedStringArray) -> void:
	for argument in args:
		if argument.begins_with("--warmup="):
			_warmup_ticks = int(argument.trim_prefix("--warmup="))
		elif argument.begins_with("--ticks="):
			_measured_ticks = int(argument.trim_prefix("--ticks="))
		elif argument.begins_with("--populations="):
			_populations = []
			for value in argument.trim_prefix("--populations=").split(",", false):
				_populations.append(int(value))
		elif argument.begins_with("--output="):
			_output_path = argument.trim_prefix("--output=")
		else:
			push_error("Unknown P0 benchmark argument: %s" % argument)
			quit(2)
			return
	assert(_warmup_ticks >= 0)
	assert(_measured_ticks > 0)
	assert(not _populations.is_empty())
	for population in _populations:
		assert(population >= 1 and population <= 1000)

func _write_report(path: String, json: String) -> void:
	var file := FileAccess.open(path, FileAccess.WRITE)
	if file == null:
		push_error("Could not write P0 report to %s (error %s)" % [path, FileAccess.get_open_error()])
		quit(3)
		return
	file.store_string(json + "\n")
	file.close()
	print("P0: wrote report to %s" % path)
