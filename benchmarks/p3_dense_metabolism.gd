extends SceneTree

const SCHEMA_VERSION: int = 1
const DEFAULT_POPULATIONS: Array[int] = [1, 16, 64, 256, 1000]
const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var _warmup_ticks: int = 3
var _measured_ticks: int = 10
var _populations: Array[int] = DEFAULT_POPULATIONS.duplicate()
var _output_path: String = ""
var _failure_message: String = ""

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	if not _parse_args(OS.get_cmdline_user_args()):
		return
	var results: Array = []
	for population in _populations:
		print("P3: comparing dense and legacy metabolism at %d cells" % population)
		var comparison: Dictionary = _compare(population)
		if comparison.is_empty():
			push_error("P3 benchmark aborted: %s" % _failure_message)
			quit(4)
			return
		results.append(comparison)
	var report: Dictionary = {
		"schema_version": SCHEMA_VERSION,
		"model_version": "microverse-m10",
		"purpose": "p3_dense_metabolism_paired_benchmark",
		"machine": _machine_manifest(),
		"parameters": {
			"warmup_ticks": _warmup_ticks,
			"measured_ticks": _measured_ticks,
			"populations": _populations
		},
		"results": results
	}
	var json: String = JSON.stringify(report, "  ")
	print("P3_RESULT_JSON_BEGIN")
	print(json)
	print("P3_RESULT_JSON_END")
	if not _output_path.is_empty() and not _write_report(_output_path, json):
		quit(3)
		return
	quit(0)

func _compare(population: int) -> Dictionary:
	var dense = _create_simulation(population, true)
	var legacy = _create_simulation(population, false)
	if dense.compiled_reactions == null or legacy.compiled_reactions == null:
		return _fail_comparison(
			population,
			"compiled reaction network did not load; close Godot and rebuild the project cache"
		)
	for _tick in range(_warmup_ticks):
		dense.step(1)
		legacy.step(1)
		if dense.checksum() != legacy.checksum():
			return _fail_comparison(population, "dense and M10 checksums diverged during warmup")

	var dense_usec: int = 0
	var legacy_usec: int = 0
	# Alternate which solver runs first on every measured tick so CPU boost,
	# thermal drift and neighboring system load do not consistently favor one.
	for tick in range(_measured_ticks):
		if tick % 2 == 0:
			dense_usec += _measure_step(dense)
			legacy_usec += _measure_step(legacy)
		else:
			legacy_usec += _measure_step(legacy)
			dense_usec += _measure_step(dense)
		if dense.checksum() != legacy.checksum():
			return _fail_comparison(population, "dense and M10 checksums diverged at measured tick %d" % (tick + 1))
	if dense.event_log != legacy.event_log:
		return _fail_comparison(population, "dense and M10 event histories diverged")
	return {
		"initial_population": population,
		"final_population": dense.population_size(),
		"measured_ticks": _measured_ticks,
		"dense_elapsed_usec": dense_usec,
		"legacy_elapsed_usec": legacy_usec,
		"dense_ticks_per_second": float(_measured_ticks) * 1000000.0 / float(dense_usec),
		"legacy_ticks_per_second": float(_measured_ticks) * 1000000.0 / float(legacy_usec),
		"speedup": float(legacy_usec) / float(dense_usec),
		"exact_checksum": dense.checksum(),
		"event_count": dense.event_log.size()
	}

func _fail_comparison(population: int, message: String) -> Dictionary:
	_failure_message = "%d cells: %s" % [population, message]
	return {}

func _measure_step(sim) -> int:
	var started: int = Time.get_ticks_usec()
	sim.step(1)
	return maxi(1, Time.get_ticks_usec() - started)

func _create_simulation(population: int, dense_solver: bool):
	var config = SimConfigScript.new()
	config.seed = 930000 + population
	config.max_cells = maxi(population * 2, population + 2)
	config.mutation_enabled = false
	config.performance_profiling_enabled = false
	config.metabolic_use_dense_solver = dense_solver
	config.validate()
	var sim = SimulationEngineScript.new(config)
	_seed_grid(sim, population)
	return sim

func _seed_grid(sim, population: int) -> void:
	var spacing: float = 1.05
	var margin: float = 0.55
	var columns: int = int(floor((float(sim.world.width - 1) - 2.0 * margin) / spacing)) + 1
	for index in range(population):
		var x: float = margin + float(index % columns) * spacing
		var y: float = margin + float(index / columns) * spacing
		sim.seed_ancestor(Vector2(x, y))

func _machine_manifest() -> Dictionary:
	var version: Dictionary = Engine.get_version_info()
	return {
		"godot_version": String(version.get("string", "unknown")),
		"os_name": OS.get_name(),
		"os_version": OS.get_version(),
		"processor_name": OS.get_processor_name(),
		"processor_count": OS.get_processor_count()
	}

func _parse_args(args: PackedStringArray) -> bool:
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
			push_error("Unknown P3 benchmark argument: %s" % argument)
			quit(2)
			return false
	if _warmup_ticks < 0 or _measured_ticks <= 0 or _populations.is_empty():
		push_error("Invalid P3 benchmark parameters")
		quit(2)
		return false
	return true

func _write_report(path: String, json: String) -> bool:
	var file := FileAccess.open(path, FileAccess.WRITE)
	if file == null:
		push_error("Could not write P3 report to %s" % path)
		return false
	file.store_string(json + "\n")
	file.close()
	return true
