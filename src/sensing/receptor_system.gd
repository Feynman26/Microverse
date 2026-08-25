extends RefCounted
class_name ReceptorSystem

const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

# Upper sequence nibble is a coarse membrane-localization motif. No ancestral
# M10 protein has Exxx, so receptor capability is not gifted. Binding itself uses
# the remaining 12 bits and ordinary Hamming affinity against molecular ligand
# signatures. Molecule names never determine whether something is a "signal".
const RECEPTOR_LOCALIZATION_MASK: int = 0xF000
const RECEPTOR_LOCALIZATION_VALUE: int = 0xE000
const BINDING_MASK: int = 0x0FFF
const DEFAULT_MAX_DISTANCE: int = 4
const DEFAULT_DISTANCE_DECAY: float = 0.70

static var _affinity_cache: Dictionary = {}

static func has_receptor_localization(protein_signature: int) -> bool:
	return (protein_signature & RECEPTOR_LOCALIZATION_MASK) == RECEPTOR_LOCALIZATION_VALUE

static func binding_signature(protein_signature: int) -> int:
	return protein_signature & BINDING_MASK

static func ligand_binding_signature(metabolite_id: String) -> int:
	assert(MetaboliteCatalogScript.has_extracellular_field(metabolite_id))
	return MetaboliteCatalogScript.ligand_signature(metabolite_id) & BINDING_MASK

static func hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & BINDING_MASK
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

static func affinity(
	protein_signature: int,
	metabolite_id: String,
	max_distance: int = DEFAULT_MAX_DISTANCE,
	distance_decay: float = DEFAULT_DISTANCE_DECAY
) -> float:
	if not has_receptor_localization(protein_signature):
		return 0.0
	var protein_binding: int = binding_signature(protein_signature)
	var target: int = ligand_binding_signature(metabolite_id)
	var key: String = "%d:%d:%d:%.12f" % [protein_binding, target, max_distance, distance_decay]
	if _affinity_cache.has(key):
		return float(_affinity_cache[key])
	var distance: int = hamming_distance(protein_binding, target)
	var result: float = 0.0 if distance > max_distance else exp(-distance_decay * float(distance))
	_affinity_cache[key] = result
	return result

# Competitive receptor occupancy from a common extracellular snapshot. For each
# receptor sequence, all compatible extracellular compounds contribute binding
# drive and compete for the finite receptor cohort. The authoritative output is
# activated molecular amount by receptor sequence, not ligand labels.
static func occupancy(
	expression_state: Dictionary,
	ligand_concentrations: Dictionary,
	cell_volume: float,
	reference_protein_count: float,
	binding_km: float,
	max_distance: int = DEFAULT_MAX_DISTANCE,
	distance_decay: float = DEFAULT_DISTANCE_DECAY
) -> Dictionary:
	assert(cell_volume > 0.0)
	assert(reference_protein_count > 0.0)
	assert(binding_km > 0.0)
	var by_signature: Dictionary = {}
	var receptor_total: float = 0.0
	var bound_total: float = 0.0
	var loci: Array = expression_state.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			if not has_receptor_localization(signature):
				continue
			var abundance: float = maxf(0.0, float(cohorts[signature_variant]))
			if abundance <= 0.0:
				continue
			receptor_total += abundance
			var drive: float = 0.0
			var ligand_ids: Array = ligand_concentrations.keys()
			ligand_ids.sort()
			for ligand_variant in ligand_ids:
				var metabolite_id: String = String(ligand_variant)
				if not MetaboliteCatalogScript.has_extracellular_field(metabolite_id):
					continue
				var concentration: float = maxf(0.0, float(ligand_concentrations[ligand_variant]))
				if concentration <= 0.0:
					continue
				var molecular_affinity: float = affinity(signature, metabolite_id, max_distance, distance_decay)
				if molecular_affinity <= 0.0:
					continue
				drive += molecular_affinity * concentration / binding_km
			var occupied_fraction: float = drive / (1.0 + drive)
			var bound: float = abundance * occupied_fraction
			by_signature[signature] = float(by_signature.get(signature, 0.0)) + bound
			bound_total += bound
	return {
		"by_signature": by_signature,
		"receptor_total": receptor_total,
		"bound_total": bound_total,
		"bound_fraction": bound_total / maxf(1e-12, receptor_total),
		"reference_concentration": reference_protein_count * cell_volume
	}

static func maintenance_cost(receptor_total: float, dt: float, atp_cost_per_protein_per_min: float) -> float:
	assert(receptor_total >= 0.0 and dt >= 0.0 and atp_cost_per_protein_per_min >= 0.0)
	return receptor_total * dt * atp_cost_per_protein_per_min
