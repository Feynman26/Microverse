extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const CellMechanicsScript = preload("res://src/physics/cell_mechanics.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_growing_colony_generates_local_nutrient_gradient()
	_test_position_alone_changes_ecological_outcome()
	_test_local_reproduction_supports_clonal_sectors()
	_test_same_seed_replays_spatial_trajectory()

	if failures == 0:
		print("PASS: %d M6 spatial-ecology tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M6 spatial-ecology tests failed" % [failures, tests_run])
		quit(1)

func _gradient_config():
	var config = SimConfigScript.new()
	config.world_width = 20
	config.world_height = 20
	config.max_cells = 8
	config.mutation_enabled = false
	config.initial_glucose = 20.0
	config.initial_oxygen = 20.0
	config.initial_nitrogen = 12.0
	config.initial_phosphorus = 8.0
	config.glucose_diffusion = 0.08
	config.oxygen_diffusion = 0.08
	config.nitrogen_diffusion = 0.08
	config.phosphorus_diffusion = 0.08
	return config

func _test_growing_colony_generates_local_nutrient_gradient() -> void:
	var config = _gradient_config()
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor(Vector2(10.0, 10.0))
	sim.step(1600)

	_assert_true(sim.population_size() > 1, "resource-fed M6 colony produces descendants before gradient measurement")
	_assert_true(sim.maximum_generation() >= 1, "gradient assay contains a genuinely growing colony")

	var stats: Dictionary = _near_far_field_means(sim, "glucose", 2.5, 6.0)
	_assert_true(int(stats["near_count"]) > 0 and int(stats["far_count"]) > 0, "gradient assay contains both colony-near and colony-far lattice sites")
	_assert_true(float(stats["far_mean"]) > float(stats["near_mean"]) + 0.01, "local colony consumption leaves less glucose near cells than far from the colony")
	_assert_true(sim.world.get_field("glucose").maximum_value() > sim.world.get_field("glucose").minimum_value() + 0.01, "colony-driven consumption creates a nonuniform extracellular glucose field")

func _near_far_field_means(sim, field_name: String, near_distance: float, far_distance: float) -> Dictionary:
	var field = sim.world.get_field(field_name)
	var near_sum: float = 0.0
	var far_sum: float = 0.0
	var near_count: int = 0
	var far_count: int = 0
	for y in range(field.height):
		for x in range(field.width):
			var point := Vector2(float(x), float(y))
			var nearest: float = INF
			for cell in sim.cells:
				nearest = minf(nearest, point.distance_to(cell.position))
			if nearest <= near_distance:
				near_sum += float(field.get_value(x, y))
				near_count += 1
			elif nearest >= far_distance:
				far_sum += float(field.get_value(x, y))
				far_count += 1
	return {
		"near_mean": near_sum / maxf(1.0, float(near_count)),
		"far_mean": far_sum / maxf(1.0, float(far_count)),
		"near_count": near_count,
		"far_count": far_count
	}

func _localized_patch_config():
	var config = SimConfigScript.new()
	config.seed = 611903
	config.world_width = 16
	config.world_height = 16
	config.max_cells = 6
	config.mutation_enabled = false
	config.initial_glucose = 0.0
	config.initial_oxygen = 0.0
	config.initial_nitrogen = 0.0
	config.initial_phosphorus = 0.0
	config.glucose_diffusion = 0.0
	config.oxygen_diffusion = 0.0
	config.nitrogen_diffusion = 0.0
	config.phosphorus_diffusion = 0.0
	return config

func _build_patch_sim(position: Vector2):
	var sim = SimulationEngineScript.new(_localized_patch_config())
	for field_name in ["glucose", "oxygen", "nitrogen", "phosphorus"]:
		var field = sim.world.get_field(field_name)
		for y in range(2, 7):
			for x in range(2, 7):
				field.set_value(x, y, 100.0)
	sim.seed_ancestor(position)
	return sim

func _test_position_alone_changes_ecological_outcome() -> void:
	# The two worlds are otherwise identical, including seed and molecular state.
	# Only the ancestor's position differs relative to one fixed resource patch.
	var rich_sim = _build_patch_sim(Vector2(4.0, 4.0))
	var poor_sim = _build_patch_sim(Vector2(12.0, 12.0))
	rich_sim.step(1000)
	poor_sim.step(1000)

	_assert_true(rich_sim.population_size() > 0, "cell placed inside the resource patch remains viable")
	_assert_true(rich_sim.maximum_generation() >= 1, "local resource access can propagate through chemistry to reproduction")
	_assert_true(poor_sim.population_size() == 0, "same ancestor outside the nondiffusing resource patch reaches ordinary extinction")
	_assert_true(rich_sim.population_size() > poor_sim.population_size(), "changing position alone changes the fixed-horizon ecological outcome")

func _test_local_reproduction_supports_clonal_sectors() -> void:
	var config = SimConfigScript.new()
	config.world_width = 24
	config.world_height = 24
	config.mechanical_relaxation_iterations = 48
	config.mechanical_overlap_tolerance = 1e-5
	config.mutation_enabled = false
	var world = WorldStateScript.new(config.world_width, config.world_height, config.grid_cell_size_um)
	var rng = DeterministicRngScript.new(481151)
	var left_origin := Vector2(8.0, 12.0)
	var right_origin := Vector2(14.0, 12.0)
	var cells: Array = [
		_make_divide_ready_cell(1, left_origin, 1001, config),
		_make_divide_ready_cell(2, right_origin, 2001, config)
	]
	var next_id: int = 3

	for _generation in range(3):
		var next_cells: Array = []
		for cell in cells:
			cell.set_pool("BIO", config.division_volume * config.biomass_units_per_volume)
			cell.volume = config.division_volume
			cell.set_pool("ATP", 10.0)
			var daughters: Array = cell.create_daughters(next_id, next_id + 1, _generation + 1, rng, world, config)
			next_id += 2
			next_cells.append_array(daughters)
		cells = next_cells
		CellMechanicsScript.relax(cells, world, config, true)

	var left_cells: Array = []
	var right_cells: Array = []
	for cell in cells:
		var marker: int = int(cell.genome.get_gene_by_locus(1).neutral_marker)
		if marker == 1001:
			left_cells.append(cell)
		elif marker == 2001:
			right_cells.append(cell)

	_assert_true(left_cells.size() == 8 and right_cells.size() == 8, "neutral founder markers remain inherited across three local generations")
	_assert_true(CellMechanicsScript.max_overlap(cells, config) <= config.mechanical_overlap_tolerance * 1.1, "multi-lineage local reproduction remains mechanically nonoverlapping")

	var left_center: Vector2 = _centroid(left_cells)
	var right_center: Vector2 = _centroid(right_cells)
	var separation: float = left_center.distance_to(right_center)
	var left_radius: float = _maximum_distance(left_cells, left_center)
	var right_radius: float = _maximum_distance(right_cells, right_center)
	_assert_true(left_center.x < right_center.x, "founder ordering remains spatial rather than being erased by local birth mechanics")
	_assert_true(separation > maxf(left_radius, right_radius), "descendants form two coherent clonal spatial sectors without a kin-grouping rule")
	_assert_true(left_center.distance_to(left_origin) < left_center.distance_to(right_origin), "left lineage remains associated with its local founder neighborhood")
	_assert_true(right_center.distance_to(right_origin) < right_center.distance_to(left_origin), "right lineage remains associated with its local founder neighborhood")

func _make_divide_ready_cell(cell_id: int, position: Vector2, marker: int, config):
	var cell = CellStateScript.new(cell_id, -1, 0, 0, position, config.ancestor_volume)
	cell.genome = GenomeScript.create_ancestor()
	cell.genome.get_gene_by_locus(1).neutral_marker = marker
	cell.initialize_molecular_state(config)
	cell.set_pool("BIO", config.division_volume * config.biomass_units_per_volume)
	cell.volume = config.division_volume
	cell.set_pool("ATP", 10.0)
	return cell

func _centroid(cells: Array) -> Vector2:
	var result := Vector2.ZERO
	for cell in cells:
		result += cell.position
	return result / maxf(1.0, float(cells.size()))

func _maximum_distance(cells: Array, center: Vector2) -> float:
	var maximum: float = 0.0
	for cell in cells:
		maximum = maxf(maximum, cell.position.distance_to(center))
	return maximum

func _test_same_seed_replays_spatial_trajectory() -> void:
	var first = _build_patch_sim(Vector2(4.0, 4.0))
	var second = _build_patch_sim(Vector2(4.0, 4.0))
	first.step(1000)
	second.step(1000)

	_assert_close(first.checksum(), second.checksum(), 1e-12, "same seed and same spatial environment reproduce the complete world checksum")
	_assert_true(first.event_log == second.event_log, "same seed reproduces exact birth/division/death history under M6 mechanics")
	_assert_true(first.population_size() == second.population_size(), "same seed reproduces final population size")
	if first.population_size() == second.population_size():
		var first_positions: Dictionary = _positions_by_id(first.cells)
		var second_positions: Dictionary = _positions_by_id(second.cells)
		_assert_true(first_positions.keys() == second_positions.keys(), "same seed reproduces the same surviving cell identities")
		for id_variant in first_positions.keys():
			var cell_id: int = int(id_variant)
			var a: Vector2 = first_positions[cell_id]
			var b: Vector2 = second_positions[cell_id]
			_assert_close(a.x, b.x, 1e-12, "same seed reproduces x trajectory for cell %d" % cell_id)
			_assert_close(a.y, b.y, 1e-12, "same seed reproduces y trajectory for cell %d" % cell_id)

func _positions_by_id(cells: Array) -> Dictionary:
	var result: Dictionary = {}
	for cell in cells:
		result[int(cell.id)] = cell.position
	return result

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
