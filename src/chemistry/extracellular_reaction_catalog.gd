extends RefCounted
class_name ExtracellularReactionCatalog

const ExtracellularReactionDefinitionScript = preload("res://src/chemistry/extracellular_reaction_definition.gd")

# M7 extracellular reactions remain ordinary chemistry. E01 hydrolyzes one C4
# lipid precursor into two C2 units. E02 converts tracked extracellular oxidant
# into the otherwise semantically neutral X species. Neither reaction contains
# producer, beneficiary, cooperation, cheating, detoxification-role, or
# ownership state; those interpretations are measured only by counterfactuals.
#
# E01 is distance 4 from the one-bit-accessible secreted D136 protein used in
# M7-E while remaining outside the ancestral catalytic radius.
# E02 is distance 4 from DACE. DACE itself is one ordinary coding-bit mutation
# from ancestral locus 5 (5ACE -> DACE), carries the same generic Dxxx secretion
# motif, and E02 remains distance >= 5 from every ancestral protein. D136 is
# distance 6 from E02, preventing the M7-E public-resource enzyme from receiving
# this detox activity for free.
const E01_SIGNATURE: int = 0xD0F2
const E02_SIGNATURE: int = 0xC0DE

static func create_m7_candidate() -> Array:
	return [
		ExtracellularReactionDefinitionScript.new(
			"E01",
			"extracellular_lipid_hydrolysis",
			E01_SIGNATURE,
			{"lipids": 1.0},
			{"carbon_c2": 2.0},
			1.0
		),
		ExtracellularReactionDefinitionScript.new(
			"E02",
			"extracellular_oxidant_neutralization",
			E02_SIGNATURE,
			{"oxidant": 1.0},
			{"neutral_x": 1.0},
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
