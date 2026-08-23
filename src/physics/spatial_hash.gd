extends RefCounted
class_name CellSpatialHash

# Broad-phase index for finite cell disks. A disk is inserted into every bucket
# touched by its axis-aligned bounding box. Therefore any two overlapping disks
# must share at least one bucket regardless of the chosen positive bucket size.
# The returned pair order is canonical by immutable cell ID so replacing the
# O(N^2) broad phase cannot change floating-point accumulation order.

static func candidate_pairs(cells: Array, config, bucket_size: float) -> Array:
	assert(bucket_size > 0.0)
	var buckets: Dictionary = {}
	var by_id: Dictionary = {}

	var ordered: Array = []
	for cell in cells:
		if cell.alive:
			ordered.append(cell)
	ordered.sort_custom(func(a, b): return int(a.id) < int(b.id))

	for cell in ordered:
		var cell_id: int = int(cell.id)
		by_id[cell_id] = cell
		var radius: float = float(config.ancestor_radius_grid) * sqrt(float(cell.volume) / float(config.ancestor_volume))
		var min_bucket_x: int = floori((cell.position.x - radius) / bucket_size)
		var max_bucket_x: int = floori((cell.position.x + radius) / bucket_size)
		var min_bucket_y: int = floori((cell.position.y - radius) / bucket_size)
		var max_bucket_y: int = floori((cell.position.y + radius) / bucket_size)
		for bucket_y in range(min_bucket_y, max_bucket_y + 1):
			for bucket_x in range(min_bucket_x, max_bucket_x + 1):
				var key := Vector2i(bucket_x, bucket_y)
				var ids: Array = buckets.get(key, [])
				ids.append(cell_id)
				buckets[key] = ids

	var unique_pairs: Dictionary = {}
	var bucket_keys: Array = buckets.keys()
	bucket_keys.sort_custom(func(a, b): return a.y < b.y or (a.y == b.y and a.x < b.x))
	for bucket_key in bucket_keys:
		var ids: Array = buckets[bucket_key]
		ids.sort()
		for i in range(ids.size()):
			for j in range(i + 1, ids.size()):
				var pair_key := Vector2i(int(ids[i]), int(ids[j]))
				unique_pairs[pair_key] = true

	var pair_keys: Array = unique_pairs.keys()
	pair_keys.sort_custom(func(a, b): return a.x < b.x or (a.x == b.x and a.y < b.y))
	var result: Array = []
	for pair_key in pair_keys:
		result.append([by_id[int(pair_key.x)], by_id[int(pair_key.y)]])
	return result

static func pair_keys(cells: Array, config, bucket_size: float) -> Array:
	var result: Array = []
	for pair in candidate_pairs(cells, config, bucket_size):
		result.append(Vector2i(int(pair[0].id), int(pair[1].id)))
	return result
