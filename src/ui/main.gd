extends Node2D

const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

const TILE_PX: float = 8.0
const WORLD_ORIGIN: Vector2 = Vector2(20.0, 70.0)
const BASE_SIM_MIN_PER_REAL_SEC: float = 1.0

var simulation
var paused: bool = false
var time_multiplier: float = 10.0
var _sim_minute_budget: float = 0.0
var _status_label: Label
var _help_label: Label

func _ready() -> void:
	simulation = SimulationEngineScript.new()
	simulation.seed_ancestor()
	_create_labels()
	queue_redraw()

func _process(delta: float) -> void:
	if not paused:
		_sim_minute_budget += delta * BASE_SIM_MIN_PER_REAL_SEC * time_multiplier
		var ticks_to_run: int = int(floor(_sim_minute_budget / simulation.config.tick_dt_min))
		ticks_to_run = mini(ticks_to_run, 5000)
		if ticks_to_run > 0:
			simulation.step(ticks_to_run)
			_sim_minute_budget -= float(ticks_to_run) * simulation.config.tick_dt_min
	_update_status()
	queue_redraw()

func _unhandled_input(event: InputEvent) -> void:
	if not event is InputEventKey or not event.pressed or event.echo:
		return
	match event.keycode:
		KEY_SPACE:
			paused = not paused
		KEY_1:
			time_multiplier = 1.0
		KEY_2:
			time_multiplier = 10.0
		KEY_3:
			time_multiplier = 100.0
		KEY_4:
			time_multiplier = 1000.0
		KEY_5:
			time_multiplier = 5000.0

func _draw() -> void:
	_draw_glucose_field()
	_draw_cells()

func _draw_glucose_field() -> void:
	var field = simulation.world.get_field("glucose")
	var reference: float = maxf(float(simulation.config.initial_glucose), 1e-9)
	for y in range(field.height):
		for x in range(field.width):
			var normalized: float = clampf(float(field.get_value(x, y)) / reference, 0.0, 1.0)
			var intensity: float = 0.06 + 0.34 * normalized
			var rect: Rect2 = Rect2(WORLD_ORIGIN + Vector2(float(x), float(y)) * TILE_PX, Vector2(TILE_PX, TILE_PX))
			draw_rect(rect, Color(0.04, intensity, 0.10 + 0.30 * normalized, 1.0), true)

func _draw_cells() -> void:
	for cell in simulation.cells:
		var cell_position: Vector2 = cell.position
		var cell_volume: float = float(cell.volume)
		var screen_position: Vector2 = WORLD_ORIGIN + cell_position * TILE_PX + Vector2(TILE_PX * 0.5, TILE_PX * 0.5)
		var radius: float = maxf(2.5, sqrt(cell_volume) * 4.0)
		draw_circle(screen_position, radius, Color(0.88, 0.93, 0.96, 1.0))
		draw_circle(screen_position, maxf(1.0, radius - 1.5), Color(0.16, 0.30, 0.34, 1.0))

func _create_labels() -> void:
	_status_label = Label.new()
	_status_label.position = Vector2(555.0, 45.0)
	_status_label.size = Vector2(390.0, 690.0)
	add_child(_status_label)

	_help_label = Label.new()
	_help_label.position = Vector2(20.0, 15.0)
	_help_label.text = "Microverse M5-A  |  SPACE pause  |  1: 1x  2: 10x  3: 100x  4: 1000x  5: 5000x"
	add_child(_help_label)

func _update_status() -> void:
	var glucose_total: float = simulation.world.get_field("glucose").total_amount()
	var oxygen_total: float = simulation.world.get_field("oxygen").total_amount()
	var nitrogen_total: float = simulation.world.get_field("nitrogen").total_amount()
	var phosphorus_total: float = simulation.world.get_field("phosphorus").total_amount()
	var state: String = "PAUSED" if paused else "RUNNING"
	_status_label.text = (
		"STATE: %s\n\n" % state
		+ "Experimental clock: %.0fx\n" % time_multiplier
		+ "Tick: %d\n" % simulation.tick_index
		+ "Virtual time: %.1f min\n\n" % simulation.simulation_time_min
		+ "Population: %d / %d\n" % [simulation.population_size(), simulation.config.max_cells]
		+ "Max generation: %d\n" % simulation.maximum_generation()
		+ "Genotypes alive: %d\n" % simulation.genotype_count()
		+ "Mutation events: %d\n" % simulation.mutation_event_count()
		+ "Total BIO-volume: %.3f\n\n" % simulation.total_cell_volume()
		+ "Environment\n"
		+ "  G: %.2f  O2: %.2f\n" % [glucose_total, oxygen_total]
		+ "  N: %.2f  P: %.2f\n\n" % [nitrogen_total, phosphorus_total]
		+ _focal_cell_status()
		+ "\nEvents recorded: %d\n" % simulation.event_log.size()
		+ "Seed: %d\n" % simulation.config.seed
		+ "Checksum: %.6f" % simulation.checksum()
	)

func _focal_cell_status() -> String:
	if simulation.cells.is_empty():
		return "Focal cell: none (extinction)\n"
	var cell = simulation.cells[0]
	var dominant: Dictionary = _dominant_flux(cell)
	var expression: Dictionary = cell.last_expression_summary
	var transcribed: float = float(expression.get("transcribed", 0.0))
	var translated: float = float(expression.get("translated", 0.0))
	var expression_atp: float = float(expression.get("atp_spent", 0.0))
	var cohort_count: int = _protein_cohort_count(cell)
	return (
		"Focal cell #%d  gen %d\n" % [cell.id, cell.generation]
		+ "  genotype: %d\n" % cell.genome.fingerprint()
		+ "  volume/BIO: %.3f / %.3f\n" % [cell.volume, cell.pool("BIO")]
		+ "  ATP/ADP: %.3f / %.3f\n" % [cell.pool("ATP"), cell.pool("ADP")]
		+ "  G/C3/C2: %.3f / %.3f / %.3f\n" % [cell.pool("G"), cell.pool("C3"), cell.pool("C2")]
		+ "  W1/W2: %.3f / %.3f\n" % [cell.pool("W1"), cell.pool("W2")]
		+ "  ROS/damage: %.3f / %.3f\n" % [cell.pool("ROS"), cell.damage]
		+ "  mRNA/protein: %.2f / %.2f\n" % [cell.total_mrna(), cell.total_protein()]
		+ "  protein cohorts: %d\n" % cohort_count
		+ "  this tick tx/tl: %.2f / %.2f\n" % [transcribed, translated]
		+ "  expression ATP: %.4f\n" % expression_atp
		+ "  dominant flux: %s %.5f\n\n" % [dominant["reaction_id"], dominant["flux"]]
	)

func _protein_cohort_count(cell) -> int:
	var count: int = 0
	for locus_state in cell.expression_state.values():
		count += locus_state["protein"].size()
	return count

func _dominant_flux(cell) -> Dictionary:
	var best_id: String = "none"
	var best_flux: float = 0.0
	for reaction_id in cell.last_fluxes.keys():
		var flux: float = float(cell.last_fluxes[reaction_id])
		if flux > best_flux:
			best_flux = flux
			best_id = String(reaction_id)
	return {"reaction_id": best_id, "flux": best_flux}
