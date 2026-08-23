extends RefCounted
class_name MetaboliteCatalog

# Structural C/N/P bookkeeping is deliberately separate from ATP/ADP and
# NAD/NADH energetic/redox currencies. The values are digital structural units,
# not claims about exact molecular formulas. BIO is assembled cell material: it
# remains inside the authoritative pool ledger and therefore lets physical cell
# volume be derived from chemistry instead of incremented by a growth knob.

const DEFINITIONS: Dictionary = {
	"G": {"name": "carbon_source", "C": 6, "N": 0, "P": 0},
	"C3": {"name": "carbon_intermediate_3", "C": 3, "N": 0, "P": 0},
	"C2": {"name": "carbon_intermediate_2", "C": 2, "N": 0, "P": 0},
	"W1": {"name": "waste_intermediate_1", "C": 3, "N": 0, "P": 0},
	"W2": {"name": "waste_intermediate_2", "C": 3, "N": 0, "P": 0},
	"CO2": {"name": "oxidized_carbon_waste", "C": 1, "N": 0, "P": 0},
	"NH4": {"name": "reduced_nitrogen", "C": 0, "N": 1, "P": 0},
	"P": {"name": "phosphorus_resource", "C": 0, "N": 0, "P": 1},
	"AA": {"name": "amino_acid_precursor", "C": 2, "N": 1, "P": 0},
	"LIP": {"name": "lipid_precursor", "C": 4, "N": 0, "P": 0},
	"NUC": {"name": "nucleotide_precursor", "C": 2, "N": 1, "P": 1},
	"BIO": {"name": "assembled_structural_biomass", "C": 12, "N": 4, "P": 2},
	"ATP": {"name": "energy_currency_high", "C": 0, "N": 0, "P": 0},
	"ADP": {"name": "energy_currency_low", "C": 0, "N": 0, "P": 0},
	"NAD": {"name": "redox_currency_oxidized", "C": 0, "N": 0, "P": 0},
	"NADH": {"name": "redox_currency_reduced", "C": 0, "N": 0, "P": 0},
	"O2": {"name": "electron_acceptor", "C": 0, "N": 0, "P": 0},
	"ROS": {"name": "oxidative_stress", "C": 0, "N": 0, "P": 0},
	"X": {"name": "semantically_neutral_diffusible_compound", "C": 0, "N": 0, "P": 0}
}

static func has(metabolite_id: String) -> bool:
	return DEFINITIONS.has(metabolite_id)

static func definition(metabolite_id: String) -> Dictionary:
	assert(DEFINITIONS.has(metabolite_id), "Unknown metabolite: %s" % metabolite_id)
	return DEFINITIONS[metabolite_id]

static func structural_units(metabolite_id: String) -> Dictionary:
	var item: Dictionary = definition(metabolite_id)
	return {"C": int(item["C"]), "N": int(item["N"]), "P": int(item["P"])}

static func ids() -> Array[String]:
	var result: Array[String] = []
	for key in DEFINITIONS.keys():
		result.append(String(key))
	result.sort()
	return result
