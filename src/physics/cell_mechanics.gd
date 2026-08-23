extends RefCounted
class_name CellMechanics

const SpatialHashScript = preload("res://src/physics/spatial_hash.gd")

# M6 reference mechanics. Cells are finite disks whose radius follows area
# scaling in the 2D chamber. Contact corrections are computed from one common
# position snapshot and accumulated before any cell moves (Jacobi semantics),
# so iteration order cannot grant a cell mechanical priority.

const EPSILON: float = 1e-12

static func radius_for_volume(volume: float, config) -> float:
	assert(volume > 0.0)
	return float(config.ancestor_radius_grid) * sqrt(volume / float(config.ancestor_volume))

static func radius_for_cell(cell, config) -> float:
	return radius_for_volume(float(cell.volume), config)

static func clamp_position(position: Vector2, radius: float, world) -> Vector2:
	assert(radius > 0.0)
	var min_x: float = radius
	var min_y: float = radius
	var max_x: float = float(world.width - 1) - radius
	var max_y: float = float(world.height - 1) - radius
	assert(max_x >= min_x and max_y >= min_y, "Cell disk does not fit inside chamber")
	return Vector2(
		clampf(position.x, min_x, max_x),
		clampf(position.y, min_y, max_y)
	)

static func relax(cells: Array, world, config, use_spatial_index: bool = false) -> Dictionary:
	var ordered: Array = []
	for cell in cells:
		if cell.alive:
			ordered.append(cell)
	# Canonicalizing arithmetic order by immutable ID is only a reproducibility
	# device. Pair corrections remain equal/opposite and never use ID to choose
	# magnitude, direction, survival, or priority.
	ordered.sort_custom(func(a, b): return int(a.id) < int(b.id))
	if ordered.is_empty():
		return {"iterations": 0, "max_overlap": 0.0, "contacts": 0, "candidate_pairs": 0}

	for cell in ordered:
		cell.position = clamp_position(cell.position, radius_for_cell(cell, config), world)

	var iterations_used: int = 0
	var last_contacts: int = 0
	var last_candidate_pairs: int = 0
	for iteration in range(int(config.mechanical_relaxation_iterations)):
		var snapshot: Dictionary = {}
		var displacement: Dictionary = {}
		for cell in ordered:
			snapshot[int(cell.id)] = cell.position
			displacement[int(cell.id)] = Vector2.ZERO

		var pairs: Array = _candidate_pairs(ordered, config, use_spatial_index)
		last_candidate_pairs = pairs.size()
		var contacts: int = 0
		var iteration_max_overlap: float = 0.0
		for pair in pairs:
			var first = pair[0]
			var second = pair[1]
			var first_radius: float = radius_for_cell(first, config)
			var second_radius: float = radius_for_cell(second, config)
			var first_position: Vector2 = snapshot[int(first.id)]
			var second_position: Vector2 = snapshot[int(second.id)]
			var delta: Vector2 = second_position - first_position
			var distance: float = delta.length()
			var target_distance: float = first_radius + second_radius
			var overlap: float = target_distance - distance
			if overlap <= float(config.mechanical_overlap_tolerance):
				continue
			contacts += 1
			iteration_max_overlap = maxf(iteration_max_overlap, overlap)
			var direction: Vector2
			if distance > EPSILON:
				direction = delta / distance
			else:
				direction = _degenerate_direction(first_position, target_distance)
			var correction: Vector2 = (
				direction
				* 0.5
				* overlap
				* float(config.mechanical_relaxation_fraction)
			)
			var first_displacement: Vector2 = displacement[int(first.id)]
			var second_displacement: Vector2 = displacement[int(second.id)]
			displacement[int(first.id)] = first_displacement - correction
			displacement[int(second.id)] = second_displacement + correction

		last_contacts = contacts
		if contacts == 0:
			break
		for cell in ordered:
			var cell_displacement: Vector2 = displacement[int(cell.id)]
			var next_position: Vector2 = cell.position + cell_displacement
			cell.position = clamp_position(next_position, radius_for_cell(cell, config), world)
		iterations_used = iteration + 1
		if iteration_max_overlap <= float(config.mechanical_overlap_tolerance):
			break

	return {
		"iterations": iterations_used,
		"max_overlap": max_overlap(ordered, config),
		"contacts": last_contacts,
		"candidate_pairs": last_candidate_pairs
	}

static func _candidate_pairs(ordered: Array, config, use_spatial_index: bool) -> Array:
	if use_spatial_index:
		return SpatialHashScript.candidate_pairs(
			ordered,
			config,
			float(config.mechanical_neighbor_bucket_size_grid)
		)
	var result: Array = []
	for i in range(ordered.size()):
		for j in range(i + 1, ordered.size()):
			result.append([ordered[i], ordered[j]])
	return result

static func max_overlap(cells: Array, config) -> float:
	var maximum: float = 0.0
	for i in range(cells.size()):
		var first = cells[i]
		if not first.alive:
			continue
		for j in range(i + 1, cells.size()):
			var second = cells[j]
			if not second.alive:
				continue
			var target: float = radius_for_cell(first, config) + radius_for_cell(second, config)
			maximum = maxf(maximum, target - first.position.distance_to(second.position))
	return maxf(0.0, maximum)

static func disks_within_bounds(cells: Array, world, config, tolerance: float = 1e-6) -> bool:
	for cell in cells:
		if not cell.alive:
			continue
		var radius: float = radius_for_cell(cell, config)
		if cell.position.x < radius - tolerance:
			return false
		if cell.position.y < radius - tolerance:
			return false
		if cell.position.x > float(world.width - 1) - radius + tolerance:
			return false
		if cell.position.y > float(world.height - 1) - radius + tolerance:
			return false
	return true

# Exact coincident centers have no geometric normal. The fallback is derived
# only from pair geometry/chamber coordinates, never from cell ID or insertion
# order. Production division already gives daughters a seeded nonzero axis; this
# path mainly makes externally constructed degenerate states well-defined.
static func _degenerate_direction(position: Vector2, target_distance: float) -> Vector2:
	var phase: float = sin(
		(position.x + 0.5) * 12.9898
		+ (position.y + 0.5) * 78.233
		+ target_distance * 37.719
	) * 43758.5453
	var fractional: float = phase - floor(phase)
	var angle: float = fractional * TAU
	return Vector2(cos(angle), sin(angle))
