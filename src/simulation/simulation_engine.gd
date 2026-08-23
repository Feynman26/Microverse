extends RefCounted
class_name SimulationEngine

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")

var config
var rng
var world
var cells: Array = []
var tick_index: int = 0
var simulation_time_min: float = 0.0
var next_cell_id: int = 1
var event_log: Array = []

func _init(p_config = null) -> void:
	config = p_config if p_config != null else SimConfigScript.new()
	config.validate()
	rng = DeterministicRngScript.new(config.seed)
	world = WorldStateScript.new(config.world_width, config.world_height, config.grid_cell_size_um)
	world.register_field("glucose", config.glucose_diffusion, config.initial_glucose)
	world.register_field("oxygen", config.oxygen_diffusion, config.initial_oxygen)

func seed_ancestor(position: Vector2 = Vector2(-1.0, -1.0)):
	if position.x < 0.0 or position.y < 0.0:
		position = Vector2(float(config.world_width - 1) * 0.5, float(config.world_height - 1) * 0.5)
	var cell = CellStateScript.new(_allocate_cell_id(), -1, 0, tick_index, world.clamp_position(position), config.ancestor_volume)
	cells.append(cell)
	_record_event("birth", {"cell_id": cell.id, "parent_id": -1, "generation": 0})
	return cell

func step(tick_count: int = 1) -> void:
	assert(tick_count >= 0)
	for _unused in range(tick_count):
		_step_once()

func _step_once() -> void:
	var dt: float = float(config.tick_dt_min)
	world.diffuse(dt)
	_allocate_membrane_transport(dt)

	for cell in cells:
		if cell.alive:
			cell.step_intracellular(dt, config)

	_process_deaths()
	_process_divisions()
	world.assert_nonnegative()
	tick_index += 1
	simulation_time_min += dt

# Resource allocation deliberately has three distinct phases:
# 1) all cells request from one immutable post-diffusion snapshot;
# 2) one scale factor is computed per grid site and chemical from that snapshot;
# 3) the precomputed allocations are applied.
# Recomputing the scale after each removal would make later records receive less
# solely because of iteration order, an artificial selective advantage.
func _allocate_membrane_transport(dt: float) -> void:
	var records: Array = []
	var glucose_totals: Dictionary = {}
	var oxygen_totals: Dictionary = {}

	for cell in cells:
		if not cell.alive:
			continue
		var request: Dictionary = cell.transport_requests(dt, world, config)
		var key: Vector2i = _grid_key(cell.position)
		records.append({"cell": cell, "key": key, "request": request})
		glucose_totals[key] = float(glucose_totals.get(key, 0.0)) + float(request["glucose"])
		oxygen_totals[key] = float(oxygen_totals.get(key, 0.0)) + float(request["oxygen"])

	var glucose_scales: Dictionary = {}
	var oxygen_scales: Dictionary = {}
	for key_variant in glucose_totals.keys():
		var key: Vector2i = key_variant
		var available: float = float(world.get_field("glucose").get_value(key.x, key.y))
		var total_request: float = float(glucose_totals[key])
		glucose_scales[key] = 1.0 if total_request <= available or total_request <= 0.0 else available / total_request
	for key_variant in oxygen_totals.keys():
		var key: Vector2i = key_variant
		var available: float = float(world.get_field("oxygen").get_value(key.x, key.y))
		var total_request: float = float(oxygen_totals[key])
		oxygen_scales[key] = 1.0 if total_request <= available or total_request <= 0.0 else available / total_request

	for record in records:
		var cell = record["cell"]
		var key: Vector2i = record["key"]
		var request: Dictionary = record["request"]
		var glucose_allocated: float = float(request["glucose"]) * float(glucose_scales[key])
		var oxygen_allocated: float = float(request["oxygen"]) * float(oxygen_scales[key])
		var glucose_removed: float = float(world.get_field("glucose").remove_amount(key.x, key.y, glucose_allocated))
		var oxygen_removed: float = float(world.get_field("oxygen").remove_amount(key.x, key.y, oxygen_allocated))
		cell.apply_uptake(glucose_removed, oxygen_removed)

func _process_deaths() -> void:
	var survivors: Array = []
	for cell in cells:
		if cell.alive:
			survivors.append(cell)
			continue
		var pools: Dictionary = cell.releasable_pools()
		world.release("glucose", cell.position, float(pools["glucose"]))
		world.release("oxygen", cell.position, float(pools["oxygen"]))
		_record_event("death", {
			"cell_id": cell.id,
			"generation": cell.generation,
			"reason": cell.death_reason
		})
	cells = survivors

func _process_divisions() -> void:
	var next_population: Array = []
	var projected_population: int = cells.size()

	for cell in cells:
		if not cell.alive:
			continue
		if cell.ready_to_divide(config) and projected_population < int(config.max_cells):
			var daughters: Array = cell.create_daughters(_allocate_cell_id(), _allocate_cell_id(), tick_index, rng, world, config)
			projected_population += 1
			_record_event("division", {
				"parent_id": cell.id,
				"daughter_ids": [daughters[0].id, daughters[1].id],
				"generation": cell.generation + 1
			})
			for daughter in daughters:
				_record_event("birth", {
					"cell_id": daughter.id,
					"parent_id": cell.id,
					"generation": daughter.generation
				})
			next_population.append_array(daughters)
		else:
			next_population.append(cell)

	cells = next_population

func _grid_key(position: Vector2) -> Vector2i:
	return Vector2i(
		clampi(roundi(position.x), 0, int(config.world_width) - 1),
		clampi(roundi(position.y), 0, int(config.world_height) - 1)
	)

func _allocate_cell_id() -> int:
	var result: int = next_cell_id
	next_cell_id += 1
	return result

func _record_event(kind: String, payload: Dictionary) -> void:
	var event: Dictionary = payload.duplicate(true)
	event["kind"] = kind
	event["tick"] = tick_index
	event["time_min"] = simulation_time_min
	event_log.append(event)

func population_size() -> int:
	return cells.size()

func maximum_generation() -> int:
	var result: int = 0
	for cell in cells:
		result = maxi(result, int(cell.generation))
	return result

func total_cell_volume() -> float:
	var result: float = 0.0
	for cell in cells:
		result += float(cell.volume)
	return result

func checksum() -> float:
	var result: float = float(world.checksum()) + float(tick_index) * 37.0 + simulation_time_min * 41.0
	for cell in cells:
		result += float(cell.checksum())
	result += float(rng.get_state() % 1000003) * 1e-6
	return result
