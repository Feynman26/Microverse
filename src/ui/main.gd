extends Node2D

const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const SimConfigScript = preload("res://src/core/sim_config.gd")
const DNAReplicationScript = preload("res://src/genetics/dna_replication.gd")

const TILE_PX: float = 8.0
const WORLD_ORIGIN: Vector2 = Vector2(20.0, 70.0)
const BASE_SIM_MIN_PER_REAL_SEC: float = 1.0
# Rendering and status are observers. They do not need to execute at the
# biological tick rate. A bounded simulation CPU slice keeps input/rendering
# responsive even when the requested virtual-time multiplier exceeds hardware.
const MAX_SIM_CPU_USEC_PER_FRAME: int = 12000
const MAX_TICKS_PER_FRAME: int = 256
const STATUS_REFRESH_SEC: float = 0.25
const DRAW_REFRESH_SEC: float = 0.10
const SPEED_MEASURE_SEC: float = 0.50
# This is explicitly a workstation safety guard, not a carrying capacity. The
# UI pauses the instant it is reached so the engine never continues an
# artificial "stationary phase" in which otherwise-ready divisions are blocked.
const INTERACTIVE_MAX_CELLS: int = 256

var simulation
var paused: bool = false
var compute_limit_reached: bool = false
var time_multiplier: float = 10.0
var _sim_minute_budget: float = 0.0
var _status_refresh_budget: float = 0.0
var _draw_refresh_budget: float = 0.0
var _speed_elapsed: float = 0.0
var _speed_ticks: int = 0
var _actual_ticks_per_sec: float = 0.0
var _status_label: Label
var _help_label: Label

func _ready() -> void:
	var config = SimConfigScript.new()
	config.max_cells = INTERACTIVE_MAX_CELLS
	config.validate()
	simulation = SimulationEngineScript.new(config)
	simulation.seed_ancestor()
	_create_labels()
	_update_status()
	queue_redraw()

func _process(delta: float) -> void:
	var ticks_ran: int = 0
	if not paused:
		_sim_minute_budget += delta * BASE_SIM_MIN_PER_REAL_SEC * time_multiplier
		# Do not accumulate minutes of wall-clock lag when hardware cannot meet a
		# requested multiplier. No biological tick is skipped: the achieved clock
		# simply runs slower and is reported explicitly to the user.
		var max_backlog: float = maxf(1.0, BASE_SIM_MIN_PER_REAL_SEC * time_multiplier * 0.5)
		_sim_minute_budget = minf(_sim_minute_budget, max_backlog)
		var deadline_usec: int = Time.get_ticks_usec() + MAX_SIM_CPU_USEC_PER_FRAME
		while (
			_sim_minute_budget + 1e-12 >= float(simulation.config.tick_dt_min)
			and ticks_ran < MAX_TICKS_PER_FRAME
		):
			simulation.step(1)
			_sim_minute_budget -= float(simulation.config.tick_dt_min)
			ticks_ran += 1
			if simulation.population_size() >= int(simulation.config.max_cells):
				compute_limit_reached = true
				paused = true
				_sim_minute_budget = 0.0
				break
			if Time.get_ticks_usec() >= deadline_usec:
				break

	_speed_ticks += ticks_ran
	_speed_elapsed += delta
	if _speed_elapsed >= SPEED_MEASURE_SEC:
		_actual_ticks_per_sec = float(_speed_ticks) / maxf(_speed_elapsed, 1e-9)
		_speed_elapsed = 0.0
		_speed_ticks = 0

	_status_refresh_budget += delta
	_draw_refresh_budget += delta
	if _status_refresh_budget >= STATUS_REFRESH_SEC:
		_status_refresh_budget = fmod(_status_refresh_budget, STATUS_REFRESH_SEC)
		_update_status()
	var draw_refresh: float = _current_draw_refresh_sec()
	if _draw_refresh_budget >= draw_refresh:
		_draw_refresh_budget = fmod(_draw_refresh_budget, draw_refresh)
		queue_redraw()

func _current_draw_refresh_sec() -> float:
	# At high requested clocks, rendering is deliberately sampled more sparsely.
	# This changes only observer cadence; every biological tick remains intact.
	if time_multiplier >= 1000.0:
		return 0.50
	if time_multiplier >= 100.0:
		return 0.25
	return DRAW_REFRESH_SEC

func _unhandled_input(event: InputEvent) -> void:
	if not event is InputEventKey or not event.pressed or event.echo:
		return
	match event.keycode:
		KEY_SPACE:
			# A reached computational guard cannot be resumed without changing the
			# experiment's configured guard. Restart with a larger guard instead of
			# silently allowing cap-induced biology.
			if not compute_limit_reached:
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
	_help_label.text = "Microverse M10  |  SPACE pause  |  1: 1x  2: 10x  3: 100x  4: 1000x  5: 5000x"
	add_child(_help_label)

func _update_status() -> void:
	var glucose_total: float = simulation.world.get_field("glucose").total_amount()
	var oxygen_total: float = simulation.world.get_field("oxygen").total_amount()
	var nitrogen_total: float = simulation.world.get_field("nitrogen").total_amount()
	var phosphorus_total: float = simulation.world.get_field("phosphorus").total_amount()
	var state: String = "COMPUTE LIMIT" if compute_limit_reached else ("PAUSED" if paused else "RUNNING")
	var achieved_multiplier: float = (
		_actual_ticks_per_sec * float(simulation.config.tick_dt_min) / BASE_SIM_MIN_PER_REAL_SEC
	)
	var division: Dictionary = _division_status()
	_status_label.text = (
		"STATE: %s\n\n" % state
		+ "Requested clock: %.0fx\n" % time_multiplier
		+ "Achieved clock: %.1fx (%.1f ticks/s)\n" % [achieved_multiplier, _actual_ticks_per_sec]
		+ "Tick: %d\n" % simulation.tick_index
		+ "Virtual time: %.1f min\n\n" % simulation.simulation_time_min
		+ "Population: %d / %d COMPUTE CAP\n" % [simulation.population_size(), simulation.config.max_cells]
		+ "Max generation: %d\n" % simulation.maximum_generation()
		+ "Genotypes alive: %d\n" % simulation.genotype_count()
		+ "Mutation events: %d\n" % simulation.mutation_event_count()
		+ "Total BIO-volume: %.3f\n" % simulation.total_cell_volume()
		+ "Division ready: %d\n" % int(division["ready"])
		+ "  size-ready / ATP-blocked: %d / %d\n" % [int(division["volume_ready"]), int(division["volume_ready_atp_blocked"])]
		+ "  size+ATP / DNA-blocked: %d\n\n" % int(division["volume_atp_ready_replication_blocked"])
		+ "Environment\n"
		+ "  G: %.2f  O2: %.2f\n" % [glucose_total, oxygen_total]
		+ "  N: %.2f  P: %.2f\n\n" % [nitrogen_total, phosphorus_total]
		+ _focal_cell_status()
		+ "\nEvents recorded: %d\n" % simulation.event_log.size()
		+ "Seed: %d\n" % simulation.config.seed
		+ "Checksum: %.6f" % simulation.checksum()
	)

func _division_status() -> Dictionary:
	var ready: int = 0
	var volume_ready: int = 0
	var volume_ready_atp_blocked: int = 0
	var volume_atp_ready_replication_blocked: int = 0
	for cell in simulation.cells:
		if not cell.alive:
			continue
		var has_volume: bool = float(cell.volume) >= float(simulation.config.division_volume)
		var has_atp: bool = float(cell.pool("ATP")) >= float(simulation.config.division_atp_cost)
		var has_replication: bool = (
			not bool(simulation.config.evolvable_replication_enabled)
			or DNAReplicationScript.replication_complete(cell)
		)
		if has_volume:
			volume_ready += 1
			if not has_atp:
				volume_ready_atp_blocked += 1
			elif not has_replication:
				volume_atp_ready_replication_blocked += 1
			else:
				ready += 1
	return {
		"ready": ready,
		"volume_ready": volume_ready,
		"volume_ready_atp_blocked": volume_ready_atp_blocked,
		"volume_atp_ready_replication_blocked": volume_atp_ready_replication_blocked
	}

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
		+ "  ATP/ADP: %.3f / %.3f  total %.3f\n" % [cell.pool("ATP"), cell.pool("ADP"), cell.total_adenylate()]
		+ "  NAD/NADH total: %.3f\n" % cell.total_redox_currency()
		+ "  DNA copy progress: %.3f\n" % cell.replication_progress
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
