extends RefCounted
class_name ExtracellularReactionCatalog

const ExtracellularReactionDefinitionScript = preload("res://src/chemistry/extracellular_reaction_definition.gd")

# M7-E begins with one materially balanced extracellular reaction. Extracellular
# lipid precursor (C4) is hydrolyzed into two C2 units. The product is an
# ordinary M7 secondary metabolite and may be imported by any cell whose
# realized proteome happens to support C2 transport. Nothing in this reaction
# encodes producer, beneficiary, cooperation, cheating, or ownership.
#
# The reaction signature is distance 4 from the one-bit-accessible secreted
# D136 protein used in the controlled capability assay, while remaining at least
# distance 5 from every ancestral protein. Thus the ancestral proteome is not
# accidentally an extracellular hydrolase even if protein material were ever
# placed outside by an intervention.
const E01_SIGNATURE: int = 0xD0F2

static func create_m7_candidate() -> Array:
	return [
		ExtracellularReactionDefinitionScript.new(
			"E01",
			"extracellular_lipid_hydrolysis",
			E01_SIGNATURE,
			{"lipids": 1.0},
			{"carbon_c2": 2.0},
			1.0
		)
	]

static func by_id(reactions: Array, reaction_id: String):
	for reaction in reactions:
		if reaction.reaction_id == reaction_id:
			return reaction
	return null

static func validate_unique(reactions: Array) -> void:
	var ids: Dictionary = {}
	for reaction in reactions:
		reaction.validate()
		assert(not ids.has(reaction.reaction_id), "Duplicate extracellular reaction ID: %s" % reaction.reaction_id)
		ids[reaction.reaction_id] = true
