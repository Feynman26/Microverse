extends RefCounted
class_name ExpressionSolver

const ExpressionStateScript = preload("res://src/expression/expression_state.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")

# M5 gene expression is a causal molecular layer:
# promoter -> transcription -> mRNA -> translation -> protein -> function.
# Regulation is generic. Any expressed protein can bind a gene's inherited
# regulatory motif according to signature distance. The highest protein bit
# determines whether bound protein contributes activating or repressing effect;
# no gene is tagged as a named regulator or biological behavior.

static func initialize(genome, config):
	var mrna: Dictionary = {}
	var proteins: Dictionary = {}
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		var basal: float = float(gene.promoter_strength())
		var m_ss: float = basal * float(config.transcription_rate) / float(config.mrna_decay_rate)
		var p_ss: float = m_ss * float(config.translation_rate) / float(config.protein_decay_rate)
		mrna[locus_id] = m_ss
		proteins[locus_id] = p_ss
	var state = ExpressionStateScript.new(mrna, proteins)
	state.assert_matches_genome(genome)
	return state

static func step(state, genome, metabolites: Dictionary, dt: float, rng, config) -> Dictionary:
	assert(state != null and genome != null)
	assert(dt >= 0.0)
	state.assert_matches_genome(genome)
	if dt <= 0.0:
		return {"transcription": 0.0, "translation": 0.0, "atp_cost": 0.0, "regulation": {}}

	var old_mrna: Dictionary = state.mrna.duplicate(true)
	var old_proteins: Dictionary = state.proteins.duplicate(true)
	var transcription_requests: Dictionary = {}
	var translation_requests: Dictionary = {}
	var regulation: Dictionary = {}

	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		var regulation_factor: float = _regulation_factor(gene, genome, old_proteins, config)
		regulation[locus_id] = regulation_factor
		var transcription_noise: float = _noise_multiplier(rng, float(config.expression_noise_fraction))
		var translation_noise: float = _noise_multiplier(rng, float(config.expression_noise_fraction))
		transcription_requests[locus_id] = maxf(0.0,
			float(config.transcription_rate) * float(gene.promoter_strength()) * regulation_factor * dt * transcription_noise
		)
		translation_requests[locus_id] = maxf(0.0,
			float(config.translation_rate) * float(old_mrna.get(locus_id, 0.0)) * dt * translation_noise
		)

	# Expression has explicit energetic and structural opportunity cost. mRNA
	# synthesis consumes NUC; protein synthesis consumes AA. Degradation returns
	# those precursor units, so the represented material remains auditable.
	var total_transcription: float = _sum_values(transcription_requests)
	var total_translation: float = _sum_values(translation_requests)
	var requested_atp: float = (
		total_transcription * float(config.transcription_atp_cost_per_unit)
		+ total_translation * float(config.translation_atp_cost_per_unit)
	)
	var requested_nuc: float = total_transcription * float(config.nuc_cost_per_mrna_unit)
	var requested_aa: float = total_translation * float(config.aa_cost_per_protein_unit)
	var atp_scale: float = _availability_scale(float(metabolites.get("ATP", 0.0)), requested_atp)
	var nuc_scale: float = _availability_scale(float(metabolites.get("NUC", 0.0)), requested_nuc)
	var aa_scale: float = _availability_scale(float(metabolites.get("AA", 0.0)), requested_aa)
	var transcription_scale: float = minf(atp_scale, nuc_scale)
	var translation_scale: float = minf(atp_scale, aa_scale)

	var effective_transcription: float = 0.0
	var effective_translation: float = 0.0
	var returned_nuc: float = 0.0
	var returned_aa: float = 0.0
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		var m_old: float = float(old_mrna.get(locus_id, 0.0))
		var p_old: float = float(old_proteins.get(locus_id, 0.0))
		var m_deg: float = minf(m_old, m_old * float(config.mrna_decay_rate) * dt)
		var p_deg: float = minf(p_old, p_old * float(config.protein_decay_rate) * dt)
		var m_new: float = float(transcription_requests[locus_id]) * transcription_scale
		var p_new: float = float(translation_requests[locus_id]) * translation_scale
		state.mrna[locus_id] = maxf(0.0, m_old - m_deg + m_new)
		state.proteins[locus_id] = maxf(0.0, p_old - p_deg + p_new)
		effective_transcription += m_new
		effective_translation += p_new
		returned_nuc += m_deg * float(config.nuc_cost_per_mrna_unit)
		returned_aa += p_deg * float(config.aa_cost_per_protein_unit)

	var atp_spent: float = (
		effective_transcription * float(config.transcription_atp_cost_per_unit)
		+ effective_translation * float(config.translation_atp_cost_per_unit)
	)
	metabolites["ATP"] = maxf(0.0, float(metabolites.get("ATP", 0.0)) - atp_spent)
	metabolites["ADP"] = float(metabolites.get("ADP", 0.0)) + atp_spent
	metabolites["NUC"] = maxf(0.0,
		float(metabolites.get("NUC", 0.0))
		- effective_transcription * float(config.nuc_cost_per_mrna_unit)
		+ returned_nuc
	)
	metabolites["AA"] = maxf(0.0,
		float(metabolites.get("AA", 0.0))
		- effective_translation * float(config.aa_cost_per_protein_unit)
		+ returned_aa
	)
	state.assert_nonnegative()
	return {
		"transcription": effective_transcription,
		"translation": effective_translation,
		"atp_cost": atp_spent,
		"regulation": regulation
	}

static func partition(state, ratio: float, rng, config) -> Array:
	assert(state != null)
	assert(ratio > 0.0 and ratio < 1.0)
	var first_mrna: Dictionary = {}
	var second_mrna: Dictionary = {}
	var first_proteins: Dictionary = {}
	var second_proteins: Dictionary = {}
	var loci: Array = state.mrna.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var jitter: float = float(config.expression_partition_jitter)
		var local_ratio: float = clampf(ratio + float(rng.randf_range(-jitter, jitter)), 0.05, 0.95)
		var m: float = state.mrna_for(locus_id)
		var p: float = state.protein_for(locus_id)
		first_mrna[locus_id] = m * local_ratio
		second_mrna[locus_id] = m - float(first_mrna[locus_id])
		first_proteins[locus_id] = p * local_ratio
		second_proteins[locus_id] = p - float(first_proteins[locus_id])
	return [
		ExpressionStateScript.new(first_mrna, first_proteins),
		ExpressionStateScript.new(second_mrna, second_proteins)
	]

static func _regulation_factor(target_gene, genome, protein_snapshot: Dictionary, config) -> float:
	if not bool(config.regulation_enabled):
		return 1.0
	var activator: float = 0.0
	var repressor: float = 0.0
	for regulator_gene in genome.genes:
		var abundance: float = maxf(0.0, float(protein_snapshot.get(int(regulator_gene.locus_id), 0.0)))
		if abundance <= 0.0:
			continue
		var distance: int = CatalyticLandscapeScript.hamming_distance(
			int(regulator_gene.protein_signature), int(target_gene.regulatory_signature)
		)
		if distance > int(config.regulatory_max_distance):
			continue
		var affinity: float = exp(-float(config.regulatory_distance_decay) * float(distance))
		var occupancy: float = abundance * affinity
		if (int(regulator_gene.protein_signature) & 0x8000) == 0:
			activator += occupancy
		else:
			repressor += occupancy
	var total: float = activator + repressor
	var normalized: float = (activator - repressor) / (1.0 + total)
	return clampf(
		1.0 + float(config.regulatory_gain) * normalized,
		float(config.regulatory_min_factor),
		float(config.regulatory_max_factor)
	)

static func _noise_multiplier(rng, fraction: float) -> float:
	if fraction <= 0.0:
		return 1.0
	return maxf(0.0, 1.0 + float(rng.randf_range(-fraction, fraction)))

static func _availability_scale(available: float, requested: float) -> float:
	if requested <= 0.0 or requested <= available:
		return 1.0
	return maxf(0.0, available / requested)

static func _sum_values(values: Dictionary) -> float:
	var result: float = 0.0
	for value in values.values():
		result += float(value)
	return result
