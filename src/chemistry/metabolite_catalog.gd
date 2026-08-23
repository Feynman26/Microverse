extends RefCounted
class_name MetaboliteCatalog

# Structural C/N/P bookkeeping is separate from energetic/redox currencies.
# M5 additionally assigns every small-molecule pool a 16-bit ligand signature.
# The signature has no semantic meaning: it is only a molecular-binding target
# used by the generic allosteric sensing grammar.

const DEFINITIONS: Dictionary = {
	"G": {"name": "carbon_source", "C": 6, "N": 0, "P": 0, "ligand": 0x0000},
	"C3": {"name": "carbon_intermediate_3", "C": 3, "N": 0, "P": 0, "ligand": 0x2222},
	"C2": {"name": "carbon_intermediate_2", "C": 2, "N": 0, "P": 0, "ligand": 0x4444},
	"W1": {"name": "waste_intermediate_1", "C": 3, "N": 0, "P": 0, "ligand": 0xAAAA},
	"W2": {"name": "waste_intermediate_2", "C": 3, "N": 0, "P": 0, "ligand": 0xBEEF},
	"CO2": {"name": "oxidized_carbon_waste", "C": 1, "N": 0, "P": 0, "ligand": 0x0101},
	"NH4": {"name": "reduced_nitrogen", "C": 0, "N": 1, "P": 0, "ligand": 0x8080},
	"P": {"name": "phosphorus_resource", "C": 0, "N": 0, "P": 1, "ligand": 0x1111},
	"AA": {"name": "amino_acid_precursor", "C": 2, "N": 1, "P": 0, "ligand": 0x5555},
	"LIP": {"name": "lipid_precursor", "C": 4, "N": 0, "P": 0, "ligand": 0xDEAD},
	"NUC": {"name": "nucleotide_precursor", "C": 2, "N": 1, "P": 1, "ligand": 0xFFFF},
	"BIO": {"name": "assembled_structural_biomass", "C": 12, "N": 4, "P": 2, "ligand": 0x3333},
	"ATP": {"name": "energy_currency_high", "C": 0, "N": 0, "P": 0, "ligand": 0xA5A5},
	"ADP": {"name": "energy_currency_low", "C": 0, "N": 0, "P": 0, "ligand": 0x5A5A},
	"NAD": {"name": "redox_currency_oxidized", "C": 0, "N": 0, "P": 0, "ligand": 0x0F0F},
	"NADH": {"name": "redox_currency_reduced", "C": 0, "N": 0, "P": 0, "ligand": 0xF0F0},
	"O2": {"name": "electron_acceptor", "C": 0, "N": 0, "P": 0, "ligand": 0xCCCC},
	"ROS": {"name": "oxidative_stress", "C": 0, "N": 0, "P": 0, "ligand": 0x3C3C},
	"X": {"name": "semantically_neutral_diffusible_compound", "C": 0, "N": 0, "P": 0, "ligand": 0x9696}
}

static func has(metabolite_id: String) -> bool:
	return DEFINITIONS.has(metabolite_id)

static func definition(metabolite_id: String) -> Dictionary:
	assert(DEFINITIONS.has(metabolite_id), "Unknown metabolite: %s" % metabolite_id)
	return DEFINITIONS[metabolite_id]

static func structural_units(metabolite_id: String) -> Dictionary:
	var item: Dictionary = definition(metabolite_id)
	return {"C": int(item["C"]), "N": int(item["N"]), "P": int(item["P"])}

static func ligand_signature(metabolite_id: String) -> int:
	return int(definition(metabolite_id)["ligand"]) & 0xFFFF

static func ids() -> Array[String]:
	var result: Array[String] = []
	for key in DEFINITIONS.keys(): result.append(String(key))
	result.sort()
	return result
