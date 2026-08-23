extends RefCounted
class_name MetaboliteCatalog

# C/N/P bookkeeping, ligand identity, and extracellular availability are
# independent molecular properties. An extracellular field only means the
# molecule can physically exist in the chamber; it does not grant any cell the
# machinery to import, export, sense, or benefit from it.
const DEFINITIONS: Dictionary = {
	"G": {"name": "carbon_source", "C": 6, "N": 0, "P": 0, "ligand": 0x0000, "field": "glucose"},
	"C3": {"name": "carbon_intermediate_3", "C": 3, "N": 0, "P": 0, "ligand": 0x2222, "field": "carbon_c3"},
	"C2": {"name": "carbon_intermediate_2", "C": 2, "N": 0, "P": 0, "ligand": 0x4444, "field": "carbon_c2"},
	"W1": {"name": "waste_intermediate_1", "C": 3, "N": 0, "P": 0, "ligand": 0xAAAA, "field": "waste_1"},
	"W2": {"name": "waste_intermediate_2", "C": 3, "N": 0, "P": 0, "ligand": 0xBEEF, "field": "waste_2"},
	"CO2": {"name": "oxidized_carbon_waste", "C": 1, "N": 0, "P": 0, "ligand": 0x0101, "field": "carbon_dioxide"},
	"NH4": {"name": "reduced_nitrogen", "C": 0, "N": 1, "P": 0, "ligand": 0x8080, "field": "nitrogen"},
	"P": {"name": "phosphorus_resource", "C": 0, "N": 0, "P": 1, "ligand": 0x1111, "field": "phosphorus"},
	"AA": {"name": "amino_acid_precursor", "C": 2, "N": 1, "P": 0, "ligand": 0x5555, "field": "amino_acids"},
	"LIP": {"name": "lipid_precursor", "C": 4, "N": 0, "P": 0, "ligand": 0xDEAD, "field": "lipids"},
	"NUC": {"name": "nucleotide_precursor", "C": 2, "N": 1, "P": 1, "ligand": 0xFFFF, "field": "nucleotides"},
	"BIO": {"name": "assembled_structural_biomass", "C": 12, "N": 4, "P": 2, "ligand": 0x3333},
	"ATP": {"name": "energy_currency_high", "C": 0, "N": 0, "P": 0, "ligand": 0xA5A5},
	"ADP": {"name": "energy_currency_low", "C": 0, "N": 0, "P": 0, "ligand": 0x5A5A},
	"NAD": {"name": "redox_currency_oxidized", "C": 0, "N": 0, "P": 0, "ligand": 0x0F0F},
	"NADH": {"name": "redox_currency_reduced", "C": 0, "N": 0, "P": 0, "ligand": 0xF0F0},
	"O2": {"name": "electron_acceptor", "C": 0, "N": 0, "P": 0, "ligand": 0xCCCC, "field": "oxygen"},
	"ROS": {"name": "oxidative_stress", "C": 0, "N": 0, "P": 0, "ligand": 0x3C3C, "field": "oxidant"},
	"X": {"name": "semantically_neutral_diffusible_compound", "C": 0, "N": 0, "P": 0, "ligand": 0x9696, "field": "neutral_x"}
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

static func has_extracellular_field(metabolite_id: String) -> bool:
	return has(metabolite_id) and definition(metabolite_id).has("field")

static func extracellular_field(metabolite_id: String) -> String:
	assert(has_extracellular_field(metabolite_id), "Metabolite has no extracellular field: %s" % metabolite_id)
	return String(definition(metabolite_id)["field"])

static func extracellular_metabolite_for_field(field_name: String) -> String:
	for metabolite_id in extracellular_ids():
		if extracellular_field(metabolite_id) == field_name:
			return metabolite_id
	return ""

static func extracellular_ids() -> Array[String]:
	var result: Array[String] = []
	for key in DEFINITIONS.keys():
		var metabolite_id := String(key)
		if has_extracellular_field(metabolite_id):
			result.append(metabolite_id)
	result.sort()
	return result

static func ids() -> Array[String]:
	var result: Array[String] = []
	for key in DEFINITIONS.keys(): result.append(String(key))
	result.sort()
	return result
