extends RefCounted
class_name WorldState

const ChemicalFieldScript = preload("res://src/world/chemical_field.gd")

var width: int
var height: int
var cell_size: float
var fields: Dictionary = {}
var field_order: Array[String] = []

# M7-E extracellular protein material is stored by exact 16-bit sequence. These
# are physical protein molecules, not abstract activity fields. Dynamic fields
# are iterated in sorted signature order so diffusion/checksum semantics remain
# deterministic regardless of the order in which secretory lineages appear.
var protein_fields: Dictionary = {}

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

func ensure_protein_field(protein_signature: int, diffusion: float):
	var signature: int = protein_signature & 0xFFFF
	assert(diffusion >= 0.0)
	if not protein_fields.has(signature):
		protein_fields[signature] = ChemicalFieldScript.new(width, height, cell_size, diffusion, 0.0)
	return protein_fields[signature]

func has_protein_field(protein_signature: int) -> bool:
	return protein_fields.has(protein_signature & 0xFFFF)

func get_protein_field(protein_signature: int):
	var signature: int = protein_signature & 0xFFFF
	assert(protein_fields.has(signature), "Unknown extracellular protein signature: %s" % signature)
	return protein_fields[signature]

func release_protein(protein_signature: int, position: Vector2, amount: float, diffusion: float) -> void:
	assert(amount >= 0.0)
	if amount <= 0.0:
		return
	var field = ensure_protein_field(protein_signature, diffusion)
	field.add_amount(roundi(position.x), roundi(position.y), amount)

func protein_signatures() -> Array:
	var signatures: Array = protein_fields.keys()
	signatures.sort()
	return signatures

func total_extracellular_protein() -> float:
	var result: float = 0.0
	for signature_variant in protein_signatures():
		result += float(protein_fields[int(signature_variant)].total_amount())
	return result

func diffuse(dt: float) -> void:
	for name in field_order:
		fields[name].step_diffusion(dt)
	for signature_variant in protein_signatures():
		protein_fields[int(signature_variant)].step_diffusion(dt)

func clamp_position(position: Vector2) -> Vector2:
	return Vector2(
		clampf(position.x, 0.0, float(width - 1)),
		clampf(position.y, 0.0, float(height - 1))
	)

func assert_nonnegative() -> void:
	for name in field_order:
		assert(fields[name].minimum_value() >= -1e-10, "Negative concentration in %s" % name)
	for signature_variant in protein_signatures():
		var signature: int = int(signature_variant)
		assert(protein_fields[signature].minimum_value() >= -1e-10, "Negative extracellular protein %s" % signature)

func checksum() -> float:
	var result := 0.0
	for i in range(field_order.size()):
		result += fields[field_order[i]].checksum() * float(i + 1)
	var signatures: Array = protein_signatures()
	for i in range(signatures.size()):
		var signature: int = int(signatures[i])
		result += protein_fields[signature].checksum() * float(field_order.size() + i + 1) * 0.731
		result += float(signature) * float(i + 1) * 0.00000031
	return result
