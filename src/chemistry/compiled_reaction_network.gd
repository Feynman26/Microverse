extends RefCounted
class_name CompiledReactionNetwork

const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

# Immutable, integer-indexed view of one ordered reaction list. Substrate and
# product key order is preserved so dense accumulation follows the exact M10
# floating-point operation order.

var reaction_ids: Array[String] = []
var substrate_indices: Array = []
var substrate_coefficients: Array = []
var product_indices: Array = []
var product_coefficients: Array = []

static func compile(reactions: Array):
	var result = CompiledReactionNetwork.new()
	for reaction in reactions:
		result.reaction_ids.append(String(reaction.reaction_id))
		var substrate_index_buffer := PackedInt32Array()
		var substrate_coefficient_buffer := PackedFloat64Array()
		for metabolite_variant in reaction.substrates.keys():
			var metabolite_id: String = String(metabolite_variant)
			substrate_index_buffer.append(MetaboliteCatalogScript.index_of(metabolite_id))
			substrate_coefficient_buffer.append(float(reaction.substrates[metabolite_variant]))
		result.substrate_indices.append(substrate_index_buffer)
		result.substrate_coefficients.append(substrate_coefficient_buffer)

		var product_index_buffer := PackedInt32Array()
		var product_coefficient_buffer := PackedFloat64Array()
		for metabolite_variant in reaction.products.keys():
			var metabolite_id: String = String(metabolite_variant)
			product_index_buffer.append(MetaboliteCatalogScript.index_of(metabolite_id))
			product_coefficient_buffer.append(float(reaction.products[metabolite_variant]))
		result.product_indices.append(product_index_buffer)
		result.product_coefficients.append(product_coefficient_buffer)
	result.validate(reactions)
	return result

func validate(reactions: Array) -> void:
	assert(reaction_ids.size() == reactions.size())
	assert(substrate_indices.size() == reactions.size())
	assert(substrate_coefficients.size() == reactions.size())
	assert(product_indices.size() == reactions.size())
	assert(product_coefficients.size() == reactions.size())
	for reaction_index in range(reactions.size()):
		assert(reaction_ids[reaction_index] == String(reactions[reaction_index].reaction_id))
		assert(substrate_indices[reaction_index].size() == substrate_coefficients[reaction_index].size())
		assert(product_indices[reaction_index].size() == product_coefficients[reaction_index].size())
