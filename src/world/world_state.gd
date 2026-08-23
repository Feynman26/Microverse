extends RefCounted
class_name WorldState

const ChemicalFieldScript = preload("res://src/world/chemical_field.gd")

var width: int
var height: int
var cell_size: float
var fields: Dictionary = {}
var field_order: Array[String] = []

func _init(p_width: int = 1, p_height: int = 1, p_cell_size: float = 1.0) -> void:
	width = p_width
	height = p_height
	cell_size = p_cell_size
	assert(width > 0 and height > 0 and cell_size > 0.0)

func register_field(name: String, diffusion: float, initial_value: float = 0.0) -> void:
	assert(not fields.has(name), "Chemical field already registered: %s" % name)
	fields[name] = ChemicalFieldScript.new(width, height, cell_size, diffusion, initial_value)
	field_order.append(name)

func has_field(name: String) -> bool:
	return fields.has(name)

func get_field(name: String):
	assert(fields.has(name), "Unknown chemical field: %s" % name)
	return fields[name]

func sample(name: String, position: Vector2) -> float:
	var field = get_field(name)
	return field.get_value(roundi(position.x), roundi(position.y))

func consume(name: String, position: Vector2, amount: float) -> float:
	var field = get_field(name)
	return field.remove_amount(roundi(position.x), roundi(position.y), amount)

func release(name: String, position: Vector2, amount: float) -> void:
	var field = get_field(name)
	field.add_amount(roundi(position.x), roundi(position.y), amount)

func diffuse(dt: float) -> void:
	for name in field_order:
		fields[name].step_diffusion(dt)

func clamp_position(position: Vector2) -> Vector2:
	return Vector2(
		clampf(position.x, 0.0, float(width - 1)),
		clampf(position.y, 0.0, float(height - 1))
	)

func assert_nonnegative() -> void:
	for name in field_order:
		assert(fields[name].minimum_value() >= -1e-10, "Negative concentration in %s" % name)

func checksum() -> float:
	var result := 0.0
	for i in range(field_order.size()):
		result += fields[field_order[i]].checksum() * float(i + 1)
	return result
