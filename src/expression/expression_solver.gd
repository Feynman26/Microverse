extends RefCounted
class_name ExpressionSolver

const ExpressionStateScript = preload("res://src/expression/expression_state.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

static func initialize(genome, config):
	var mrna: Dictionary = {}
	var cohorts: Dictionary = {}
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		var basal: float = float(gene.promoter_strength())
		var m_ss: float = basal * float(config.transcription_rate) / float(config.mrna_decay_rate)
		var p_ss: float = m_ss * float(config.translation_rate) / float(config.protein_decay_rate)
		mrna[locus_id] = m_ss
		cohorts[locus_id] = {int(gene.protein_signature): p_ss}
	var state = ExpressionStateScript.new(mrna, cohorts)
	state.assert_matches_genome(genome)
	return state

static func step(state, genome, metabolites: Dictionary, dt: float, rng, config) -> Dictionary:
	assert(state != null and genome != null)
	assert(dt >= 0.0)
	state.assert_matches_genome(genome)
	if dt <= 0.0:
		return {"transcription": 0.0, "translation": 0.0, "atp_cost": 0.0, "nuc_committed": 0.0, "aa_committed": 0.0, "regulation": {}}

	var old_state = state.deep_copy()
	var metabolic_snapshot: Dictionary = metabolites.duplicate(true)
	var transcription_requests: Dictionary = {}
	var translation_requests: Dictionary = {}
	var regulation: Dictionary = {}
	var mrna_decay: Dictionary = {}
	var cohort_decay: Dictionary = {} # locus -> {signature: decayed_amount}
	var returned_nuc: float = 0.0
	var returned_aa: float = 0.0

	# All decay, regulation, and synthesis requests come from one immutable
	# molecular snapshot. Material released by turnover can be reused within the
	# interval but no locus receives priority from array order.
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		var m_old: float = old_state.mrna_for(locus_id)
		var m_deg: float = minf(m_old, m_old * float(config.mrna_decay_rate) * dt)
		mrna_decay[locus_id] = m_deg
		returned_nuc += m_deg * float(config.nuc_cost_per_mrna_unit)

		var decay_for_locus: Dictionary = {}
		var old_cohorts: Dictionary = old_state.protein_cohorts.get(locus_id, {})
		for signature_variant in old_cohorts.keys():
			var signature: int = int(signature_variant)
			var p_old: float = float(old_cohorts[signature])
			var p_deg: float = minf(p_old, p_old * float(config.protein_decay_rate) * dt)
			decay_for_locus[signature] = p_deg
			returned_aa += p_deg * float(config.aa_cost_per_protein_unit)
		cohort_decay[locus_id] = decay_for_locus

		var regulation_factor: float = _regulation_factor(gene, genome, old_state, metabolic_snapshot, config)
		regulation[locus_id] = regulation_factor
		var transcription_noise: float = _noise_multiplier(rng, float(config.expression_noise_fraction))
		var translation_noise: float = _noise_multiplier(rng, float(config.expression_noise_fraction))
		transcription_requests[locus_id] = maxf(0.0,
			float(config.transcription_rate) * float(gene.promoter_strength()) * regulation_factor * dt * transcription_noise
		)
		translation_requests[locus_id] = maxf(0.0,
			float(config.translation_rate) * m_old * dt * translation_noise
		)

	var total_transcription: float = _sum_values(transcription_requests)
	var total_translation: float = _sum_values(translation_requests)
	var requested_atp: float = total_transcription * float(config.transcription_atp_cost_per_unit) + total_translation * float(config.translation_atp_cost_per_unit)
	var requested_nuc: float = total_transcription * float(config.nuc_cost_per_mrna_unit)
	var requested_aa: float = total_translation * float(config.aa_cost_per_protein_unit)
	var atp_scale: float = _availability_scale(float(metabolic_snapshot.get("ATP", 0.0)), requested_atp)
	var nuc_scale: float = _availability_scale(float(metabolic_snapshot.get("NUC", 0.0)) + returned_nuc, requested_nuc)
	var aa_scale: float = _availability_scale(float(metabolic_snapshot.get("AA", 0.0)) + returned_aa, requested_aa)
	var transcription_scale: float = minf(atp_scale, nuc_scale)
	var translation_scale: float = minf(atp_scale, aa_scale)

	var effective_transcription: float = 0.0
	var effective_translation: float = 0.0
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		state.mrna[locus_id] = maxf(0.0,
			old_state.mrna_for(locus_id) - float(mrna_decay[locus_id]) + float(transcription_requests[locus_id]) * transcription_scale
		)
		effective_transcription += float(transcription_requests[locus_id]) * transcription_scale

		var updated_cohorts: Dictionary = old_state.protein_cohorts.get(locus_id, {}).duplicate(true)
		var decay_for_locus: Dictionary = cohort_decay[locus_id]
		for signature_variant in decay_for_locus.keys():
			var signature: int = int(signature_variant)
			updated_cohorts[signature] = maxf(0.0, float(updated_cohorts.get(signature, 0.0)) - float(decay_for_locus[signature]))
		var translated: float = float(translation_requests[locus_id]) * translation_scale
		var encoded_signature: int = int(gene.protein_signature) & 0xFFFF
		updated_cohorts[encoded_signature] = float(updated_cohorts.get(encoded_signature, 0.0)) + translated
		state.protein_cohorts[locus_id] = updated_cohorts
		effective_translation += translated

	var atp_spent: float = effective_transcription * float(config.transcription_atp_cost_per_unit) + effective_translation * float(config.translation_atp_cost_per_unit)
	var nuc_committed: float = effective_transcription * float(config.nuc_cost_per_mrna_unit)
	var aa_committed: float = effective_translation * float(config.aa_cost_per_protein_unit)
	metabolites["ATP"] = maxf(0.0, float(metabolic_snapshot.get("ATP", 0.0)) - atp_spent)
	metabolites["ADP"] = float(metabolic_snapshot.get("ADP", 0.0)) + atp_spent
	metabolites["NUC"] = maxf(0.0, float(metabolic_snapshot.get("NUC", 0.0)) + returned_nuc - nuc_committed)
	metabolites["AA"] = maxf(0.0, float(metabolic_snapshot.get("AA", 0.0)) + returned_aa - aa_committed)
	state.assert_nonnegative()
	return {
		"transcription": effective_transcription,
		"translation": effective_translation,
		"atp_cost": atp_spent,
		"nuc_committed": nuc_committed,
		"aa_committed": aa_committed,
		"nuc_recycled": returned_nuc,
		"aa_recycled": returned_aa,
		"regulation": regulation
	}

static func partition(state, ratio: float, rng, config) -> Array:
	assert(state != null)
	assert(ratio > 0.0 and ratio < 1.0)
	var first_mrna: Dictionary = {}
	var second_mrna: Dictionary = {}
	var first_cohorts: Dictionary = {}
	var second_cohorts: Dictionary = {}
	var loci: Array = state.mrna.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var jitter: float = float(config.expression_partition_jitter)
		var m_ratio: float = clampf(ratio + float(rng.randf_range(-jitter, jitter)), 0.05, 0.95)
		var m: float = state.mrna_for(locus_id)
		first_mrna[locus_id] = m * m_ratio
		second_mrna[locus_id] = m - float(first_mrna[locus_id])

		var first_locus: Dictionary = {}
		var second_locus: Dictionary = {}
		var cohorts: Dictionary = state.protein_cohorts.get(locus_id, {})
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			var local_ratio: float = clampf(ratio + float(rng.randf_range(-jitter, jitter)), 0.05, 0.95)
			var amount: float = float(cohorts[signature])
			first_locus[signature] = amount * local_ratio
			second_locus[signature] = amount - float(first_locus[signature])
		first_cohorts[locus_id] = first_locus
		second_cohorts[locus_id] = second_locus
	return [ExpressionStateScript.new(first_mrna, first_cohorts), ExpressionStateScript.new(second_mrna, second_cohorts)]

static func structural_totals(state, config) -> Dictionary:
	var mrna_material: float = state.total_mrna() * float(config.nuc_cost_per_mrna_unit)
	var protein_material: float = state.total_protein() * float(config.aa_cost_per_protein_unit)
	return {"C": mrna_material * 2.0 + protein_material * 2.0, "N": mrna_material + protein_material, "P": mrna_material}

static func _regulation_factor(target_gene, genome, state_snapshot, metabolites: Dictionary, config) -> float:
	if not bool(config.regulation_enabled):
		return 1.0
	var activator: float = 0.0
	var repressor: float = 0.0
	for record_variant in state_snapshot.cohort_records():
		var record: Dictionary = record_variant
		var signature: int = int(record["signature"])
		var abundance: float = maxf(0.0, float(record["amount"]))
		var distance: int = CatalyticLandscapeScript.hamming_distance(signature, int(target_gene.regulatory_signature))
		if distance > int(config.regulatory_max_distance):
			continue
		var affinity: float = exp(-float(config.regulatory_distance_decay) * float(distance))
		var effective_abundance: float = abundance * _allosteric_factor(signature, metabolites, config)
		var occupancy: float = effective_abundance * affinity
		if (signature & 0x8000) == 0:
			activator += occupancy
		else:
			repressor += occupancy
	var total: float = activator + repressor
	var normalized: float = (activator - repressor) / (1.0 + total)
	return clampf(1.0 + float(config.regulatory_gain) * normalized, float(config.regulatory_min_factor), float(config.regulatory_max_factor))

static func _allosteric_factor(protein_signature: int, metabolites: Dictionary, config) -> float:
	if not bool(config.allostery_enabled):
		return 1.0
	var strongest_occupancy: float = 0.0
	for metabolite_id in MetaboliteCatalogScript.ids():
		var amount: float = maxf(0.0, float(metabolites.get(metabolite_id, 0.0)))
		if amount <= 0.0:
			continue
		var distance: int = CatalyticLandscapeScript.hamming_distance(protein_signature, MetaboliteCatalogScript.ligand_signature(metabolite_id))
		if distance > int(config.allosteric_max_distance):
			continue
		var affinity: float = exp(-float(config.allosteric_distance_decay) * float(distance))
		var scaled_amount: float = amount * affinity
		var occupancy: float = scaled_amount / (float(config.allosteric_km_per_volume) + scaled_amount)
		strongest_occupancy = maxf(strongest_occupancy, occupancy)
	if strongest_occupancy <= 0.0:
		return 1.0
	var direction: float = -1.0 if (protein_signature & 0x4000) != 0 else 1.0
	return clampf(1.0 + direction * float(config.allosteric_gain) * strongest_occupancy, float(config.allosteric_min_factor), float(config.allosteric_max_factor))

static func _noise_multiplier(rng, fraction: float) -> float:
	if fraction <= 0.0 or rng == null:
		return 1.0
	return maxf(0.0, 1.0 + float(rng.randf_range(-fraction, fraction)))

static func _availability_scale(available: float, requested: float) -> float:
	if requested <= 0.0 or requested <= available:
		return 1.0
	return maxf(0.0, available / requested)

static func _sum_values(values: Dictionary) -> float:
	var result: float = 0.0
	for value in values.values(): result += float(value)
	return result
