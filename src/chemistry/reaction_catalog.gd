extends RefCounted
class_name ReactionCatalog

const ReactionDefinitionScript = preload("res://src/chemistry/reaction_definition.gd")

# Candidate M4 network. Signatures are deliberately positioned relative to the
# M3 ancestor so core reactions begin strongly accessible while alternative
# waste-recovery routes sit just outside the active affinity radius. The
# landscape probe must validate this arrangement before physiology depends on it.

static func create_m4_candidate() -> Array:
	return [
		ReactionDefinitionScript.new(
			"R01", "carbon_activation", 0x1357,
			{"G": 1.0, "ADP": 2.0, "NAD": 2.0},
			{"C3": 2.0, "ATP": 2.0, "NADH": 2.0}, 1.0
		),
		ReactionDefinitionScript.new(
			"R02", "oxidative_carbon_processing", 0x2468,
			{"C3": 1.0, "NAD": 1.0, "ADP": 1.0},
			{"C2": 1.0, "CO2": 1.0, "NADH": 1.0, "ATP": 1.0}, 1.0
		),
		ReactionDefinitionScript.new(
			"R03", "oxidative_phosphorylation", 0x369C,
			{"NADH": 1.0, "O2": 1.0, "ADP": 3.0},
			{"NAD": 1.0, "ATP": 3.0, "ROS": 1.0}, 1.0
		),
		ReactionDefinitionScript.new(
			"R04", "fermentative_redox_relief", 0x48AA,
			{"C3": 1.0, "NADH": 1.0, "ADP": 1.0},
			{"W1": 1.0, "NAD": 1.0, "ATP": 1.0}, 0.8
		),
		ReactionDefinitionScript.new(
			"R05", "waste_1_recovery", 0x48B2,
			{"W1": 1.0, "NAD": 1.0, "ADP": 1.0},
			{"C2": 1.0, "CO2": 1.0, "NADH": 1.0, "ATP": 1.0}, 0.8
		),
		ReactionDefinitionScript.new(
			"R06", "amino_precursor_synthesis", 0x5ACE,
			{"C2": 1.0, "NH4": 1.0, "ATP": 1.0},
			{"AA": 1.0, "ADP": 1.0}, 0.8
		),
		ReactionDefinitionScript.new(
			"R07", "lipid_precursor_synthesis", 0x6BDF,
			{"C2": 2.0, "ATP": 2.0},
			{"LIP": 1.0, "ADP": 2.0}, 0.7
		),
		ReactionDefinitionScript.new(
			"R08", "nucleotide_precursor_synthesis", 0x7CE1,
			{"C2": 1.0, "NH4": 1.0, "P": 1.0, "ATP": 2.0},
			{"NUC": 1.0, "ADP": 2.0}, 0.7
		),
		ReactionDefinitionScript.new(
			"R09", "oxidative_damage_control", 0x8DF2,
			{"ROS": 1.0, "ATP": 1.0},
			{"ADP": 1.0}, 0.6
		),
		ReactionDefinitionScript.new(
			"R10", "overflow_waste_2", 0x9E04,
			{"C3": 1.0, "NADH": 1.0},
			{"W2": 1.0, "NAD": 1.0}, 0.6
		),
		ReactionDefinitionScript.new(
			"R11", "waste_2_oxidative_recovery", 0x9E1C,
			{"W2": 1.0, "O2": 1.0, "ADP": 1.0},
			{"C2": 1.0, "CO2": 1.0, "ATP": 1.0, "ROS": 1.0}, 0.7
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
		assert(not ids.has(reaction.reaction_id), "Duplicate reaction ID: %s" % reaction.reaction_id)
		ids[reaction.reaction_id] = true
