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

func _index(x: int, y: int) -> int:
	return y * width + x

func get_value(x: int, y: int) -> float:
	return values[_index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))]

func set_value(x: int, y: int, value: float) -> void:
	assert(value >= 0.0)
	var i: int = _index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))
	var prior: float = values[i]
	values[i] = value
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

func replace_values(new_values: PackedFloat64Array) -> void:
	assert(new_values.size() == width * height)
	values = new_values.duplicate()
	_buffer = values.duplicate()
	_all_zero = true
	for value in values:
		assert(value >= -1e-10)
		if value != 0.0:
			_all_zero = false

func is_all_zero() -> bool:
	return _all_zero

func add_amount(x: int, y: int, amount: float) -> void:
	assert(amount >= 0.0)
	if amount <= 0.0:
		return
	var i := _index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))
	values[i] += amount
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
	var result := INF
	for value in values:
		result = minf(result, value)
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
			_buffer[i] = maxf(0.0, next_value)

	var old_values := values
	values = _buffer
	_buffer = old_values

func checksum() -> float:
	# Weighted checksum catches spatial differences that total mass cannot.
	if _all_zero:
		return 0.0
	var result := 0.0
	for i in range(values.size()):
		result += values[i] * float(i + 1)
	return result
