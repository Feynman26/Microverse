extends RefCounted
class_name RegulatorySystem

const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")

# Stable digital ligand signatures. These are molecular identities, not semantic
# behavior labels. Any protein sequence can bind any ligand according to the same
# Hamming-affinity geometry used elsewhere. M5-B only uses local extracellular
# abundance to modulate regulatory activity; richer signalling is deferred to M11.
const LIGAND_SIGNATURES: Dictionary = {
	"G": 0x9A2B,
	"O2": 0xC4D1,
	"NH4": 0x6E35,
	"P": 0xB817,
	"X": 0xDA90
}

# Protein sequence metadata is deliberately derived from bits of the same
# 16-bit signature that also controls catalytic and binding affinities. This
# creates molecular pleiotropy without a `regulator_type` flag.
# bit 15: regulatory polarity (1 activation, 0 repression)
# bits 13-14: regulatory effect magnitude
# bit 12: environmental-gating mode
# bits 10-11: mRNA stability class
# bits 8-9: protein stability class

static func binding_affinity(protein_signature: int, promoter_motif: int) -> float:
	return CatalyticLandscapeScript.affinity(protein_signature, promoter_motif)

static func regulatory_polarity(protein_signature: int) -> float:
	return 1.0 if (protein_signature & 0x8000) != 0 else -1.0

static func regulatory_strength(protein_signature: int) -> float:
	var code: int = (protein_signature >> 13) & 0x3
	return 0.25 + 0.25 * float(code)

static func environment_gated(protein_signature: int) -> bool:
	return (protein_signature & 0x1000) != 0

static func mrna_decay_multiplier(coding_signature: int) -> float:
	var code: int = (coding_signature >> 10) & 0x3
	return [0.60, 0.85, 1.15, 1.60][code]

static func protein_decay_multiplier(coding_signature: int) -> float:
	var code: int = (coding_signature >> 8) & 0x3
	return [0.55, 0.85, 1.20, 1.70][code]

static func ligand_affinity(protein_signature: int, ligand_id: String) -> float:
	assert(LIGAND_SIGNATURES.has(ligand_id), "Unknown regulatory ligand %s" % ligand_id)
	return CatalyticLandscapeScript.affinity(protein_signature, int(LIGAND_SIGNATURES[ligand_id]))

static func environment_activation(protein_signature: int, environment: Dictionary, config) -> float:
	if not environment_gated(protein_signature):
		return 1.0
	var strongest_bound_fraction: float = 0.0
	var km: float = maxf(1e-12, float(config.regulatory_sensor_ligand_km))
	for ligand_id in LIGAND_SIGNATURES.keys():
		var concentration: float = maxf(0.0, float(environment.get(ligand_id, 0.0)))
		var concentration_fraction: float = concentration / (km + concentration)
		var bound_fraction: float = ligand_affinity(protein_signature, String(ligand_id)) * concentration_fraction
		strongest_bound_fraction = maxf(strongest_bound_fraction, bound_fraction)
	var leak: float = clampf(float(config.regulatory_sensor_leak), 0.0, 1.0)
	return leak + (1.0 - leak) * strongest_bound_fraction

static func network_snapshot(genome, expression_state: Dictionary, environment: Dictionary, config) -> Dictionary:
	var result: Dictionary = {}
	var targets: Array = genome.genes.duplicate()
	targets.sort_custom(func(a, b): return int(a.locus_id) < int(b.locus_id))
	for target_gene in targets:
		result[int(target_gene.locus_id)] = _target_regulation(genome, expression_state, target_gene, environment, config)
	return result

static func _target_regulation(genome, expression_state: Dictionary, target_gene, environment: Dictionary, config) -> Dictionary:
	var activation: float = 0.0
	var repression: float = 0.0
	var strongest_abs_contribution: float = 0.0
	var strongest_edge: Dictionary = {
		"source_locus": -1,
		"protein_signature": -1,
		"binding_affinity": 0.0,
		"occupancy": 0.0,
		"signed_contribution": 0.0
	}
	var binding_km: float = maxf(1e-12, float(config.regulatory_binding_km))
	var source_genes: Array = genome.genes.duplicate()
	source_genes.sort_custom(func(a, b): return int(a.locus_id) < int(b.locus_id))
	for source_gene in source_genes:
		var locus_id: int = int(source_gene.locus_id)
		if not expression_state.has(locus_id):
			continue
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			var affinity_value: float = binding_affinity(signature, int(target_gene.promoter_binding_motif))
			if affinity_value <= 0.0:
				continue
			var abundance: float = maxf(0.0, float(cohorts[signature_variant])) / float(config.expression_reference_protein_count)
			var active_abundance: float = abundance * environment_activation(signature, environment, config)
			var binding_drive: float = active_abundance * affinity_value
			var occupancy: float = binding_drive / (binding_km + binding_drive)
			var contribution: float = occupancy * regulatory_strength(signature) * regulatory_polarity(signature)
			if contribution >= 0.0:
				activation += contribution
			else:
				repression += -contribution
			if absf(contribution) > strongest_abs_contribution:
				strongest_abs_contribution = absf(contribution)
				strongest_edge = {
					"source_locus": locus_id,
					"protein_signature": signature,
					"binding_affinity": affinity_value,
					"occupancy": occupancy,
					"signed_contribution": contribution
				}

	var gain: float = float(config.regulatory_gain)
	var net: float = gain * (activation - repression)
	var multiplier: float = exp(clampf(net, -8.0, 8.0))
	multiplier = clampf(multiplier, float(config.regulatory_min_multiplier), float(config.regulatory_max_multiplier))
	return {
		"multiplier": multiplier,
		"activation": activation,
		"repression": repression,
		"net": activation - repression,
		"strongest_edge": strongest_edge
	}

static func active_edge_count(snapshot: Dictionary, tolerance: float = 1e-12) -> int:
	var result: int = 0
	for item in snapshot.values():
		if absf(float(item["net"])) > tolerance:
			result += 1
	return result

static func strongest_network_edge(snapshot: Dictionary) -> Dictionary:
	var best: Dictionary = {"target_locus": -1, "source_locus": -1, "signed_contribution": 0.0}
	var best_abs: float = 0.0
	var target_ids: Array = snapshot.keys()
	target_ids.sort()
	for target_variant in target_ids:
		var target_id: int = int(target_variant)
		var edge: Dictionary = snapshot[target_id]["strongest_edge"]
		var magnitude: float = absf(float(edge["signed_contribution"]))
		if magnitude > best_abs:
			best_abs = magnitude
			best = edge.duplicate(true)
			best["target_locus"] = target_id
	return best
