extends RefCounted
class_name ReactionDefinition

const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

var reaction_id: String
var reaction_name: String
var signature: int
var substrates: Dictionary
var products: Dictionary
var catalytic_ceiling: float
var compartment: String

func _init(
	p_reaction_id: String = "",
	p_reaction_name: String = "",
	p_signature: int = 0,
	p_substrates: Dictionary = {},
	p_products: Dictionary = {},
	p_catalytic_ceiling: float = 1.0,
	p_compartment: String = "cytosol"
) -> void:
	reaction_id = p_reaction_id
	reaction_name = p_reaction_name
	signature = p_signature & 0xFFFF
	substrates = p_substrates.duplicate(true)
	products = p_products.duplicate(true)
	catalytic_ceiling = p_catalytic_ceiling
	compartment = p_compartment
	validate()

func validate() -> void:
	assert(not reaction_id.is_empty())
	assert(not reaction_name.is_empty())
	assert(signature >= 0 and signature <= 0xFFFF)
	assert(not substrates.is_empty())
	assert(not products.is_empty())
	assert(catalytic_ceiling > 0.0)
	assert(compartment in ["cytosol", "membrane", "extracellular"])
	for metabolite_id in substrates.keys():
		assert(MetaboliteCatalogScript.has(String(metabolite_id)), "Unknown substrate: %s" % metabolite_id)
		assert(float(substrates[metabolite_id]) > 0.0)
	for metabolite_id in products.keys():
		assert(MetaboliteCatalogScript.has(String(metabolite_id)), "Unknown product: %s" % metabolite_id)
		assert(float(products[metabolite_id]) > 0.0)
	assert(is_structurally_balanced(), "Reaction violates digital C/N/P balance: %s %s" % [reaction_id, structural_balance()])

func structural_balance() -> Dictionary:
	var balance: Dictionary = {"C": 0.0, "N": 0.0, "P": 0.0}
	for metabolite_id in products.keys():
		var units: Dictionary = MetaboliteCatalogScript.structural_units(String(metabolite_id))
		var coefficient: float = float(products[metabolite_id])
		for element in balance.keys():
			balance[element] = float(balance[element]) + float(units[element]) * coefficient
	for metabolite_id in substrates.keys():
		var units: Dictionary = MetaboliteCatalogScript.structural_units(String(metabolite_id))
		var coefficient: float = float(substrates[metabolite_id])
		for element in balance.keys():
			balance[element] = float(balance[element]) - float(units[element]) * coefficient
	return balance

func is_structurally_balanced(tolerance: float = 1e-12) -> bool:
	var balance: Dictionary = structural_balance()
	return absf(float(balance["C"])) <= tolerance and absf(float(balance["N"])) <= tolerance and absf(float(balance["P"])) <= tolerance

func canonical_key() -> String:
	return "%s:%d:%s>%s" % [reaction_id, signature, _canonical_stoichiometry(substrates), _canonical_stoichiometry(products)]

func _canonical_stoichiometry(stoichiometry: Dictionary) -> String:
	var keys: Array = stoichiometry.keys()
	keys.sort()
	var parts: PackedStringArray = PackedStringArray()
	for key in keys:
		parts.append("%s=%.6f" % [String(key), float(stoichiometry[key])])
	return ",".join(parts)
