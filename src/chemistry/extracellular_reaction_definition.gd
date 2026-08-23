extends RefCounted
class_name ExtracellularReactionDefinition

var reaction_id: String
var name: String
var signature: int
var substrates: Dictionary
var products: Dictionary
var catalytic_ceiling: float

func _init(
	p_reaction_id: String,
	p_name: String,
	p_signature: int,
	p_substrates: Dictionary,
	p_products: Dictionary,
	p_catalytic_ceiling: float = 1.0
) -> void:
	reaction_id = p_reaction_id
	name = p_name
	signature = p_signature & 0xFFFF
	substrates = p_substrates.duplicate(true)
	products = p_products.duplicate(true)
	catalytic_ceiling = p_catalytic_ceiling
	validate()

func validate() -> void:
	assert(not reaction_id.is_empty(), "Extracellular reaction ID cannot be empty")
	assert(signature >= 0 and signature <= 0xFFFF)
	assert(not substrates.is_empty(), "Extracellular reaction requires at least one substrate")
	assert(not products.is_empty(), "Extracellular reaction requires at least one product")
	assert(catalytic_ceiling > 0.0)
	for field_name in substrates.keys():
		assert(not String(field_name).is_empty())
		assert(float(substrates[field_name]) > 0.0)
	for field_name in products.keys():
		assert(not String(field_name).is_empty())
		assert(float(products[field_name]) > 0.0)
