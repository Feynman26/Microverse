extends RefCounted
class_name ChemicalField

var width: int
var height: int
var cell_size: float
var diffusion_coefficient: float
var values := PackedFloat64Array()
var _buffer := PackedFloat64Array()

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

func _index(x: int, y: int) -> int:
	return y * width + x

func get_value(x: int, y: int) -> float:
	return values[_index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))]

func set_value(x: int, y: int, value: float) -> void:
	assert(value >= 0.0)
	values[_index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))] = value

func add_amount(x: int, y: int, amount: float) -> void:
	assert(amount >= 0.0)
	var i := _index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))
	values[i] += amount

func remove_amount(x: int, y: int, requested: float) -> float:
	assert(requested >= 0.0)
	var i := _index(clampi(x, 0, width - 1), clampi(y, 0, height - 1))
	var removed := minf(values[i], requested)
	values[i] -= removed
	if values[i] < 1e-12:
		values[i] = 0.0
	return removed

func total_amount() -> float:
	var total := 0.0
	for value in values:
		total += value
	return total

func minimum_value() -> float:
	var result := INF
	for value in values:
		result = minf(result, value)
	return result

func maximum_value() -> float:
	var result := 0.0
	for value in values:
		result = maxf(result, value)
	return result

# Explicit 5-point finite-difference diffusion with reflecting (no-flux)
# boundaries. Reflecting boundaries are represented by substituting the center
# value for the missing outside neighbor, preserving the closed chamber's mass.
func step_diffusion(dt: float) -> void:
	if diffusion_coefficient == 0.0 or dt == 0.0:
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
	var result := 0.0
	for i in range(values.size()):
		result += values[i] * float(i + 1)
	return result
