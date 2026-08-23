extends RefCounted
class_name ExtracellularReactionCatalog

const ExtracellularReactionDefinitionScript = preload("res://src/chemistry/extracellular_reaction_definition.gd")

# M7 extracellular reactions remain ordinary chemistry. E01 hydrolyzes one C4
# lipid precursor into two C2 units. E02 converts tracked extracellular oxidant
# into the otherwise semantically neutral X species. E03 performs the reverse
# tracked conversion, allowing damaging extracellular chemistry to arise from a
# secreted catalyst without adding any toxin/combat behavior state.
#
# Each capability is sequence-separated: D136 -> E01, DACE -> E02, DE03 -> E03
# are distance 4 pairs, while the matching secreted proteins are outside the
# active radius of the other two extracellular reactions. Each secreted protein
# is itself one ordinary coding-bit mutation from an ancestral locus, and every
# reaction remains outside the catalytic radius of the unmodified ancestor.
const E01_SIGNATURE: int = 0xD0F2
const E02_SIGNATURE: int = 0xC0DE
const E03_SIGNATURE: int = 0xC202

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
		),
		ExtracellularReactionDefinitionScript.new(
			"E03",
			"extracellular_oxidant_generation",
			E03_SIGNATURE,
			{"neutral_x": 1.0},
			{"oxidant": 1.0},
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
