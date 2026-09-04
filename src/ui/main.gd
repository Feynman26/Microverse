extends Node2D

const InteractiveSimulationRuntimeScript = preload("res://src/runtime/interactive_simulation_runtime.gd")
const SimConfigScript = preload("res://src/core/sim_config.gd")

const TILE_PX: float = 8.0
const WORLD_ORIGIN: Vector2 = Vector2(20.0, 70.0)
# Explicit workstation safety guard, never a biological carrying capacity.
const INTERACTIVE_MAX_CELLS: int = 256

var _runtime
var _snapshot: Dictionary = {}
var _snapshot_sequence: int = -1
var _status_label: Label
var _help_label: Label

func _ready() -> void:
	_create_labels()
	var config = SimConfigScript.new()
	config.max_cells = INTERACTIVE_MAX_CELLS
	config.validate()
	_runtime = InteractiveSimulationRuntimeScript.new()
	var error: Error = _runtime.start(config, false, 1.0)
	if error != OK:
		_status_label.text = "STATE: RUNTIME START FAILED\n\nError: %s" % error_string(error)

func _exit_tree() -> void:
	if _runtime != null:
		_runtime.stop()

func _process(_delta: float) -> void:
	if _runtime == null:
		return
	var candidate: Dictionary = _runtime.latest_snapshot()
	if candidate.is_empty():
		return
	var runtime_state: Dictionary = candidate.get("runtime", {})
	var sequence: int = int(runtime_state.get("sequence", -1))
	if sequence == _snapshot_sequence:
		return
	_snapshot = candidate
	_snapshot_sequence = sequence
	_update_status()
	queue_redraw()

func _unhandled_input(event: InputEvent) -> void:
	if _runtime == null or not event is InputEventKey or not event.pressed or event.echo:
		return
	match event.keycode:
		KEY_SPACE:
			_runtime.toggle_paused()
		KEY_1:
			_runtime.set_requested_multiplier(1.0)
		KEY_2:
			_runtime.set_requested_multiplier(10.0)
		KEY_3:
			_runtime.set_requested_multiplier(100.0)
		KEY_4:
			_runtime.set_requested_multiplier(1000.0)
		KEY_5:
			_runtime.set_requested_multiplier(5000.0)

func _draw() -> void:
	if _snapshot.is_empty():
		return
	_draw_glucose_field()
	_draw_cells()

func _draw_glucose_field() -> void:
	var width: int = int(_snapshot.get("world_width", 0))
	var height: int = int(_snapshot.get("world_height", 0))
	var values: PackedFloat64Array = _snapshot.get("glucose_values", PackedFloat64Array())
	if values.size() != width * height:
		return
	var reference: float = maxf(float(_snapshot.get("initial_glucose", 0.0)), 1e-9)
	for y in range(height):
		for x in range(width):
			var normalized: float = clampf(float(values[y * width + x]) / reference, 0.0, 1.0)
			var intensity: float = 0.06 + 0.34 * normalized
			var rect := Rect2(
				WORLD_ORIGIN + Vector2(float(x), float(y)) * TILE_PX,
				Vector2(TILE_PX, TILE_PX)
			)
			draw_rect(rect, Color(0.04, intensity, 0.10 + 0.30 * normalized, 1.0), true)

func _draw_cells() -> void:
	var positions: PackedVector2Array = _snapshot.get("cell_positions", PackedVector2Array())
	var volumes: PackedFloat32Array = _snapshot.get("cell_volumes", PackedFloat32Array())
	for index in range(mini(positions.size(), volumes.size())):
		var screen_position: Vector2 = (
			WORLD_ORIGIN
			+ positions[index] * TILE_PX
			+ Vector2(TILE_PX * 0.5, TILE_PX * 0.5)
		)
		var radius: float = maxf(2.5, sqrt(float(volumes[index])) * 4.0)
		draw_circle(screen_position, radius, Color(0.88, 0.93, 0.96, 1.0))
		draw_circle(screen_position, maxf(1.0, radius - 1.5), Color(0.16, 0.30, 0.34, 1.0))

func _create_labels() -> void:
	_status_label = Label.new()
	_status_label.position = Vector2(555.0, 45.0)
	_status_label.size = Vector2(390.0, 690.0)
	_status_label.text = "STATE: STARTING"
	add_child(_status_label)

	_help_label = Label.new()
	_help_label.position = Vector2(20.0, 15.0)
	_help_label.text = "Microverse M10 / P1 async  |  SPACE pause  |  1: 1x  2: 10x  3: 100x  4: 1000x  5: 5000x"
	add_child(_help_label)

func _update_status() -> void:
	var runtime_state: Dictionary = _snapshot.get("runtime", {})
	var environment: Dictionary = _snapshot.get("environment", {})
	var division: Dictionary = _snapshot.get("division", {})
	var worker_active: bool = bool(runtime_state.get("worker_active", false))
	var compute_limit: bool = bool(runtime_state.get("compute_limit_reached", false))
	var paused: bool = bool(runtime_state.get("paused", false))
	var overloaded: bool = bool(runtime_state.get("overloaded", false))
	var state: String = "RUNNING"
	if not worker_active:
		state = "STOPPED"
	elif compute_limit:
		state = "COMPUTE LIMIT"
	elif paused:
		state = "PAUSED"
	elif overloaded:
		state = "RUNNING / OVERLOADED"

	_status_label.text = (
		"STATE: %s\n\n" % state
		+ "Requested clock: %.0fx (%.1f ticks/s)\n" % [
			float(runtime_state.get("requested_multiplier", 1.0)),
			float(runtime_state.get("requested_ticks_per_second", 0.0))
		]
		+ "Achieved clock: %.1fx (%.1f ticks/s)\n" % [
			float(runtime_state.get("achieved_multiplier", 0.0)),
			float(runtime_state.get("actual_ticks_per_second", 0.0))
		]
		+ "Backlog: %.1f ticks / %.2f sim min\n" % [
			float(runtime_state.get("backlog_ticks", 0.0)),
			float(runtime_state.get("backlog_min", 0.0))
		]
		+ "Unserved clock demand: %.1f sim min\n" % float(runtime_state.get("unserved_requested_min", 0.0))
		+ "Tick: %d\n" % int(_snapshot.get("tick", 0))
		+ "Virtual time: %.1f min\n\n" % float(_snapshot.get("simulation_time_min", 0.0))
		+ "Population: %d / %d COMPUTE CAP\n" % [
			int(_snapshot.get("population", 0)),
			int(_snapshot.get("max_cells", 0))
		]
		+ "Max generation: %d\n" % int(_snapshot.get("maximum_generation", 0))
		+ "Genotypes alive: %d\n" % int(_snapshot.get("genotype_count", 0))
		+ "Mutation events: %d\n" % int(_snapshot.get("mutation_event_count", 0))
		+ "Total BIO-volume: %.3f\n" % float(_snapshot.get("total_cell_volume", 0.0))
		+ "Division ready: %d\n" % int(division.get("ready", 0))
		+ "  size-ready / ATP-blocked: %d / %d\n" % [
			int(division.get("volume_ready", 0)),
			int(division.get("volume_ready_atp_blocked", 0))
		]
		+ "  size+ATP / DNA-blocked: %d\n\n" % int(division.get("volume_atp_ready_replication_blocked", 0))
		+ "Environment\n"
		+ "  G: %.2f  O2: %.2f\n" % [
			float(environment.get("glucose", 0.0)),
			float(environment.get("oxygen", 0.0))
		]
		+ "  N: %.2f  P: %.2f\n\n" % [
			float(environment.get("nitrogen", 0.0)),
			float(environment.get("phosphorus", 0.0))
		]
		+ _focal_cell_status()
		+ "\nEvents recorded: %d\n" % int(_snapshot.get("event_count", 0))
		+ "Seed: %d\n" % int(_snapshot.get("seed", 0))
		+ "Checksum: %.6f" % float(_snapshot.get("checksum", 0.0))
	)

func _focal_cell_status() -> String:
	var focal: Dictionary = _snapshot.get("focal_cell", {})
	if focal.is_empty():
		return "Focal cell: none (extinction)\n"
	var dominant: Dictionary = focal.get("dominant_flux", {})
	return (
		"Focal cell #%d  gen %d\n" % [int(focal.get("id", 0)), int(focal.get("generation", 0))]
		+ "  genotype: %d\n" % int(focal.get("genotype", 0))
		+ "  volume/BIO: %.3f / %.3f\n" % [float(focal.get("volume", 0.0)), float(focal.get("bio", 0.0))]
		+ "  ATP/ADP: %.3f / %.3f  total %.3f\n" % [
			float(focal.get("atp", 0.0)),
			float(focal.get("adp", 0.0)),
			float(focal.get("adenylate", 0.0))
		]
		+ "  NAD/NADH total: %.3f\n" % float(focal.get("redox_currency", 0.0))
		+ "  DNA copy progress: %.3f\n" % float(focal.get("replication_progress", 0.0))
		+ "  G/C3/C2: %.3f / %.3f / %.3f\n" % [
			float(focal.get("g", 0.0)),
			float(focal.get("c3", 0.0)),
			float(focal.get("c2", 0.0))
		]
		+ "  W1/W2: %.3f / %.3f\n" % [float(focal.get("w1", 0.0)), float(focal.get("w2", 0.0))]
		+ "  ROS/damage: %.3f / %.3f\n" % [float(focal.get("ros", 0.0)), float(focal.get("damage", 0.0))]
		+ "  mRNA/protein: %.2f / %.2f\n" % [float(focal.get("mrna", 0.0)), float(focal.get("protein", 0.0))]
		+ "  protein cohorts: %d\n" % int(focal.get("protein_cohorts", 0))
		+ "  this tick tx/tl: %.2f / %.2f\n" % [
			float(focal.get("transcribed", 0.0)),
			float(focal.get("translated", 0.0))
		]
		+ "  expression ATP: %.4f\n" % float(focal.get("expression_atp", 0.0))
		+ "  dominant flux: %s %.5f\n\n" % [
			String(dominant.get("reaction_id", "none")),
			float(dominant.get("flux", 0.0))
		]
	)
