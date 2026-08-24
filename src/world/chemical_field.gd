extends RefCounted
class_name ChemicalField

var width: int
var height: int
var cell_size: float
var diffusion_coefficient: float
var values := PackedFloat64Array()
var _buffer := PackedFloat64Array()
# Exact state cache only: true means every stored value is bit-exact 0.0.
# It is never inferred approximately. Any positive local write invalidates it.
var _all_zero: bool = true
# M10 exact minimum cache. Diffusion already visits every lattice value, so it
# records the resulting minimum in that same pass. Pure removal can only lower
# a local value and therefore updates the minimum exactly in O(1). A local
# increase of the site that held the known minimum marks the cache dirty; only
# that uncommon case falls back to a full scan on the next invariant check.
var _minimum_cache: float = 0.0
var _minimum_dirty: bool = false

func _init(p_width: int = 1, p_height: int = 1, p_cell_size: float = 1.0, p_diffusion: float = 0.0, initial_value: float = 0.0) -> void:
	width = p_width
	height = p_height
	cell_size = p_cell_size
	diffusion_coefficient = p_diffusion
	assert(width > 0 and height > 0)
	assert(cell_size > 0.0)
	assert(diffusion_coefficient >= 0.0)
	assert(initial_value >= 0.0)
	values.resize(width * height)
	_buffer.resize(width * height)
	values.fill(initial_value)
	_buffer.fill(initial_value)
	_all_zero = initial_value == 0.0
	_minimum_cache = initial_value
	_minimum_dirty = false

func _index(x: int, y: int) -> int:
	return y * width + x

func get_value(x: int, y: int) -> float:
	return values[_index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))]

func set_value(x: int, y: int, value: float) -> void:
	assert(value >= 0.0)
	var i: int = _index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))
	var prior: float = values[i]
	values[i] = value
	if not _minimum_dirty:
		if value < _minimum_cache:
			_minimum_cache = value
		elif prior == _minimum_cache and value > prior:
			# There may or may not be another site with the old minimum.
			_minimum_dirty = true
	if value > 0.0:
		_all_zero = false
	elif prior > 0.0:
		# A local zero write cannot prove that every other lattice site is zero.
		_all_zero = false

func fill_uniform(value: float) -> void:
	assert(value >= 0.0)
	values.fill(value)
	_buffer.fill(value)
	_all_zero = value == 0.0
	_minimum_cache = value
	_minimum_dirty = false

func replace_values(new_values: PackedFloat64Array) -> void:
	assert(new_values.size() == width * height)
	values = new_values.duplicate()
	_buffer = values.duplicate()
	_all_zero = true
	_minimum_cache = INF
	for value in values:
		assert(value >= -1e-10)
		_minimum_cache = minf(_minimum_cache, value)
		if value != 0.0:
			_all_zero = false
	if values.is_empty():
		_minimum_cache = 0.0
	_minimum_dirty = false

func is_all_zero() -> bool:
	return _all_zero

func add_amount(x: int, y: int, amount: float) -> void:
	assert(amount >= 0.0)
	if amount <= 0.0:
		return
	var i := _index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))
	var prior: float = values[i]
	values[i] += amount
	if not _minimum_dirty and prior == _minimum_cache:
		# Increasing a minimum site may expose a different minimum elsewhere.
		_minimum_dirty = true
	_all_zero = false

func remove_amount(x: int, y: int, requested: float) -> float:
	assert(requested >= 0.0)
	if _all_zero or requested <= 0.0:
		return 0.0
	var i := _index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))
	var removed := minf(values[i], requested)
	values[i] -= removed
	if values[i] < 1e-12:
		values[i] = 0.0
	if not _minimum_dirty:
		_minimum_cache = minf(_minimum_cache, values[i])
	# Do not rescan to rediscover all-zero state after a local removal. Keeping
	# false is conservative and changes performance only, never chemistry.
	return removed

func total_amount() -> float:
	if _all_zero:
		return 0.0
	var total := 0.0
	for value in values:
		total += value
	return total

func minimum_value() -> float:
	if _all_zero:
		return 0.0
	if not _minimum_dirty:
		return _minimum_cache
	var result := INF
	for value in values:
		result = minf(result, value)
	_minimum_cache = result
	_minimum_dirty = false
	return result

func maximum_value() -> float:
	if _all_zero:
		return 0.0
	var result := 0.0
	for value in values:
		result = maxf(result, value)
	return result

# Explicit 5-point finite-difference diffusion with reflecting (no-flux)
# boundaries. Reflecting boundaries are represented by substituting the center
# value for the missing outside neighbor, preserving the closed chamber's mass.
func step_diffusion(dt: float) -> void:
	# Diffusion of an exactly zero field is exactly a no-op. This skips only work
	# whose mathematical result is known bit-for-bit, never a small concentration.
	if _all_zero or diffusion_coefficient == 0.0 or dt == 0.0:
		return
	assert(dt > 0.0)
	var alpha := diffusion_coefficient * dt / (cell_size * cell_size)
	assert(alpha <= 0.25 + 1e-12, "Unstable explicit 2D diffusion step")
	var next_minimum: float = INF

	for y in range(height):
		for x in range(width):
			var i := _index(x, y)
			var center := values[i]
			var left := values[_index(x - 1, y)] if x > 0 else center
			var right := values[_index(x + 1, y)] if x < width - 1 else center
			var up := values[_index(x, y - 1)] if y > 0 else center
			var down := values[_index(x, y + 1)] if y < height - 1 else center
			var next_value := center + alpha * (left + right + up + down - 4.0 * center)
			assert(next_value >= -1e-10, "Diffusion produced a materially negative concentration")
			var stored_value: float = maxf(0.0, next_value)
			_buffer[i] = stored_value
			next_minimum = minf(next_minimum, stored_value)

	var old_values := values
	values = _buffer
	_buffer = old_values
	_minimum_cache = next_minimum
	_minimum_dirty = false

func checksum() -> float:
	# Weighted checksum catches spatial differences that total mass cannot.
	if _all_zero:
		return 0.0
	var result := 0.0
	for i in range(values.size()):
		result += values[i] * float(i + 1)
	return result