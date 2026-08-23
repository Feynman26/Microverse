extends RefCounted
class_name ExpressionSystem

const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

# M5-A state remains unchanged:
# locus_id -> {"mrna": {coding_signature -> amount}, "protein": {coding_signature -> amount}}
# M5-B only changes how transcription propensity is calculated. Regulation is
# derived from physical protein cohorts + inherited promoter motifs + ordinary
# intracellular ligand concentrations; there are no named response APIs.

static func create_equilibrium_state(genome, config) -> Dictionary:
	var state: Dictionary = {}
	var tx_max: float = float(config.transcription_max_events_per_min)
	var mrna_decay: float = float(config.mrna_decay_rate_per_min)
	var translation: float = float(config.translation_events_per_mrna_per_min)
	var protein_decay: float = float(config.protein_decay_rate_per_min)
	for gene in genome.genes:
		var promoter: float = float(gene.promoter_strength())
		var signature: int = int(gene.protein_signature)
		var mrna: float = tx_max * promoter / mrna_decay
		var protein: float = translation * mrna / protein_decay
		state[int(gene.locus_id)] = {"mrna": {signature: mrna}, "protein": {signature: protein}}
	return state

static func step(state: Dictionary, genome, pools: Dictionary, dt: float, rng, config) -> Dictionary:
	assert(dt >= 0.0)
	if dt <= 0.0:
		return _empty_summary()
	_validate_state(state, genome)

	var snapshot: Dictionary = state.duplicate(true)
	var pool_snapshot: Dictionary = pools.duplicate(true)
	var tx_proposals: Dictionary = {}
	var translation_proposals: Dictionary = {}
	var mrna_decay_proposals: Dictionary = {}
	var protein_decay_proposals: Dictionary = {}
	var regulation_factors: Dictionary = {}
	var total_tx_events: float = 0.0
	var total_translation_events: float = 0.0

	var genes: Array = genome.genes.duplicate()
	genes.sort_custom(func(a, b): return int(a.locus_id) < int(b.locus_id))
	for gene in genes:
		var locus_id: int = int(gene.locus_id)
		var current_signature: int = int(gene.protein_signature)
		var prior: Dictionary = snapshot[locus_id]
		var mrna_cohorts: Dictionary = prior["mrna"]
		var protein_cohorts: Dictionary = prior["protein"]

		var regulation_factor: float = _regulation_factor(gene, snapshot, pool_snapshot, config)
		regulation_factors[locus_id] = regulation_factor
		var tx_lambda: float = (
			float(config.transcription_max_events_per_min)
			* float(gene.promoter_strength())
			* regulation_factor
			* dt
		)
		var tx_events: float = float(rng.poisson(tx_lambda))
		tx_proposals[locus_id] = {"signature": current_signature, "events": tx_events}
		total_tx_events += tx_events

		var tl_by_signature: Dictionary = {}
		var md_by_signature: Dictionary = {}
		var mrna_signatures: Array = mrna_cohorts.keys()
		mrna_signatures.sort()
		for signature_variant in mrna_signatures:
			var signature: int = int(signature_variant)
			var mrna_amount: float = maxf(0.0, float(mrna_cohorts[signature_variant]))
			var tl_lambda: float = float(config.translation_events_per_mrna_per_min) * mrna_amount * dt
			var md_lambda: float = float(config.mrna_decay_rate_per_min) * mrna_amount * dt
			var tl_events: float = float(rng.poisson(tl_lambda))
			var md_events: float = minf(mrna_amount, float(rng.poisson(md_lambda)))
			tl_by_signature[signature] = tl_events
			md_by_signature[signature] = md_events
			total_translation_events += tl_events
		translation_proposals[locus_id] = tl_by_signature
		mrna_decay_proposals[locus_id] = md_by_signature

		var pd_by_signature: Dictionary = {}
		var protein_signatures: Array = protein_cohorts.keys()
		protein_signatures.sort()
		for signature_variant in protein_signatures:
			var signature: int = int(signature_variant)
			var protein_amount: float = maxf(0.0, float(protein_cohorts[signature_variant]))
			var pd_lambda: float = float(config.protein_decay_rate_per_min) * protein_amount * dt
			pd_by_signature[signature] = minf(protein_amount, float(rng.poisson(pd_lambda)))
		protein_decay_proposals[locus_id] = pd_by_signature

	var tx_atp_demand: float = total_tx_events * float(config.transcription_atp_cost_per_event)
	var tl_atp_demand: float = total_translation_events * float(config.translation_atp_cost_per_event)
	var total_atp_demand: float = tx_atp_demand + tl_atp_demand
	var available_atp: float = maxf(0.0, float(pool_snapshot.get("ATP", 0.0)))
	var atp_scale: float = 1.0 if total_atp_demand <= available_atp or total_atp_demand <= 0.0 else available_atp / total_atp_demand
	var nuc_demand: float = total_tx_events * float(config.transcription_nuc_cost_per_event)
	var aa_demand: float = total_translation_events * float(config.translation_aa_cost_per_event)
	var available_nuc: float = maxf(0.0, float(pool_snapshot.get("NUC", 0.0)))
	var available_aa: float = maxf(0.0, float(pool_snapshot.get("AA", 0.0)))
	var nuc_scale: float = 1.0 if nuc_demand <= available_nuc or nuc_demand <= 0.0 else available_nuc / nuc_demand
	var aa_scale: float = 1.0 if aa_demand <= available_aa or aa_demand <= 0.0 else available_aa / aa_demand
	var tx_scale: float = minf(atp_scale, nuc_scale)
	var translation_scale: float = minf(atp_scale, aa_scale)

	var accepted_tx: float = 0.0
	var accepted_translation: float = 0.0
	var degraded_mrna: float = 0.0
	var degraded_protein: float = 0.0
	for gene in genes:
		var locus_id: int = int(gene.locus_id)
		var prior: Dictionary = snapshot[locus_id]
		var next_mrna: Dictionary = prior["mrna"].duplicate(true)
		var next_protein: Dictionary = prior["protein"].duplicate(true)
		for signature_variant in mrna_decay_proposals[locus_id].keys():
			var signature: int = int(signature_variant)
			var amount: float = float(mrna_decay_proposals[locus_id][signature_variant])
			next_mrna[signature] = maxf(0.0, float(next_mrna.get(signature, 0.0)) - amount)
			degraded_mrna += amount
		for signature_variant in protein_decay_proposals[locus_id].keys():
			var signature: int = int(signature_variant)
			var amount: float = float(protein_decay_proposals[locus_id][signature_variant])
			next_protein[signature] = maxf(0.0, float(next_protein.get(signature, 0.0)) - amount)
			degraded_protein += amount
		var tx_record: Dictionary = tx_proposals[locus_id]
		var tx_signature: int = int(tx_record["signature"])
		var tx_events: float = float(tx_record["events"]) * tx_scale
		next_mrna[tx_signature] = float(next_mrna.get(tx_signature, 0.0)) + tx_events
		accepted_tx += tx_events
		for signature_variant in translation_proposals[locus_id].keys():
			var signature: int = int(signature_variant)
			var events: float = float(translation_proposals[locus_id][signature_variant]) * translation_scale
			next_protein[signature] = float(next_protein.get(signature, 0.0)) + events
			accepted_translation += events
		_prune_zero_cohorts(next_mrna)
		_prune_zero_cohorts(next_protein)
		state[locus_id] = {"mrna": next_mrna, "protein": next_protein}

	var synthesis_atp: float = accepted_tx * float(config.transcription_atp_cost_per_event) + accepted_translation * float(config.translation_atp_cost_per_event)
	MetabolicSolverScript.spend_atp(pools, synthesis_atp)
	pools["NUC"] = maxf(0.0, float(pools.get("NUC", 0.0)) - accepted_tx * float(config.transcription_nuc_cost_per_event))
	pools["AA"] = maxf(0.0, float(pools.get("AA", 0.0)) - accepted_translation * float(config.translation_aa_cost_per_event))
	pools["NUC"] = float(pools["NUC"]) + degraded_mrna * float(config.transcription_nuc_cost_per_event)
	pools["AA"] = float(pools["AA"]) + degraded_protein * float(config.translation_aa_cost_per_event)

	_validate_state(state, genome)
	MetabolicSolverScript.assert_nonnegative(pools)
	return {
		"transcribed": accepted_tx,
		"translated": accepted_translation,
		"mrna_degraded": degraded_mrna,
		"protein_degraded": degraded_protein,
		"atp_spent": synthesis_atp,
		"tx_scale": tx_scale,
		"translation_scale": translation_scale,
		"regulation": regulation_factors
	}

# Generic promoter regulation. Historical protein cohorts retain their original
# sequence and therefore their own regulatory/allosteric identity after mutation.
static func _regulation_factor(target_gene, state_snapshot: Dictionary, pools: Dictionary, config) -> float:
	if not bool(config.regulation_enabled):
		return 1.0
	var activator: float = 0.0
	var repressor: float = 0.0
	var loci: Array = state_snapshot.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = state_snapshot[locus_id]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			var distance: int = _hamming_distance(signature, int(target_gene.regulatory_signature))
			if distance > int(config.regulatory_max_distance):
				continue
			var abundance: float = maxf(0.0, float(cohorts[signature_variant])) / float(config.expression_reference_protein_count)
			var affinity: float = exp(-float(config.regulatory_distance_decay) * float(distance))
			var occupancy: float = abundance * affinity * _allosteric_factor(signature, pools, config)
			if (signature & 0x8000) == 0:
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

# Any metabolite can modulate any compatible protein; molecule names never
# determine biological interpretation. The strongest compatible ligand is used
# in M5-B to keep allostery bounded and inspectable.
static func _allosteric_factor(protein_signature: int, pools: Dictionary, config) -> float:
	if not bool(config.allostery_enabled):
		return 1.0
	var strongest: float = 0.0
	for metabolite_id in MetaboliteCatalogScript.ids():
		var amount: float = maxf(0.0, float(pools.get(metabolite_id, 0.0)))
		if amount <= 0.0:
			continue
		var distance: int = _hamming_distance(protein_signature, MetaboliteCatalogScript.ligand_signature(metabolite_id))
		if distance > int(config.allosteric_max_distance):
			continue
		var affinity: float = exp(-float(config.allosteric_distance_decay) * float(distance))
		var compatible_amount: float = amount * affinity
		var occupancy: float = compatible_amount / (float(config.allosteric_km) + compatible_amount)
		strongest = maxf(strongest, occupancy)
	if strongest <= 0.0:
		return 1.0
	var direction: float = -1.0 if (protein_signature & 0x4000) != 0 else 1.0
	return clampf(
		1.0 + direction * float(config.allosteric_gain) * strongest,
		float(config.allosteric_min_factor),
		float(config.allosteric_max_factor)
	)

static func _hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & 0xFFFF
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

static func current_gene_protein(state: Dictionary, gene) -> float:
	var locus_id: int = int(gene.locus_id)
	if not state.has(locus_id): return 0.0
	return maxf(0.0, float(state[locus_id]["protein"].get(int(gene.protein_signature), 0.0)))

static func protein_cohorts_for_locus(state: Dictionary, locus_id: int) -> Dictionary:
	if not state.has(locus_id): return {}
	return state[locus_id]["protein"]

static func total_mrna(state: Dictionary) -> float:
	var result: float = 0.0
	for item in state.values():
		for amount in item["mrna"].values(): result += float(amount)
	return result

static func total_protein(state: Dictionary) -> float:
	var result: float = 0.0
	for item in state.values():
		for amount in item["protein"].values(): result += float(amount)
	return result

static func partition(state: Dictionary, ratio: float, rng, config) -> Array:
	assert(ratio > 0.0 and ratio < 1.0)
	var first: Dictionary = {}
	var second: Dictionary = {}
	var noise_scale: float = float(config.expression_partition_noise_scale)
	var loci: Array = state.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		first[locus_id] = {"mrna": {}, "protein": {}}
		second[locus_id] = {"mrna": {}, "protein": {}}
		for species_name in ["mrna", "protein"]:
			var cohorts: Dictionary = state[locus_id][species_name]
			var signatures: Array = cohorts.keys()
			signatures.sort()
			for signature_variant in signatures:
				var signature: int = int(signature_variant)
				var amount: float = maxf(0.0, float(cohorts[signature_variant]))
				var share: float = _noisy_partition_ratio(ratio, amount, noise_scale, rng)
				first[locus_id][species_name][signature] = amount * share
				second[locus_id][species_name][signature] = amount * (1.0 - share)
	return [first, second]

static func _noisy_partition_ratio(base_ratio: float, amount: float, scale: float, rng) -> float:
	if amount <= 0.0 or scale <= 0.0: return base_ratio
	var attenuation: float = 1.0 / sqrt(maxf(1.0, amount))
	return clampf(base_ratio + float(rng.randf_range(-scale, scale)) * attenuation, 0.02, 0.98)

static func structural_storage_totals(state: Dictionary, config) -> Dictionary:
	var result: Dictionary = {"C": 0.0, "N": 0.0, "P": 0.0}
	var nuc_units: Dictionary = MetaboliteCatalogScript.structural_units("NUC")
	var aa_units: Dictionary = MetaboliteCatalogScript.structural_units("AA")
	var mrna_material: float = total_mrna(state) * float(config.transcription_nuc_cost_per_event)
	var protein_material: float = total_protein(state) * float(config.translation_aa_cost_per_event)
	for element in result.keys():
		result[element] = mrna_material * float(nuc_units[element]) + protein_material * float(aa_units[element])
	return result

static func checksum(state: Dictionary) -> float:
	var loci: Array = state.keys()
	loci.sort()
	var result: float = 0.0
	for i in range(loci.size()):
		var locus_id: int = int(loci[i])
		var item: Dictionary = state[locus_id]
		result += float(locus_id) * 0.071
		for species_index in range(2):
			var species_name: String = "mrna" if species_index == 0 else "protein"
			var cohorts: Dictionary = item[species_name]
			var signatures: Array = cohorts.keys()
			signatures.sort()
			for j in range(signatures.size()):
				var signature: int = int(signatures[j])
				result += float(signature) * float(species_index + 1) * 0.0000007
				result += float(cohorts[signature]) * float((i + 3) * (j + 5) * (species_index + 1)) * 0.00031
	return result

static func _validate_state(state: Dictionary, genome) -> void:
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		assert(state.has(locus_id), "Missing expression state for locus %d" % locus_id)
		for species_name in ["mrna", "protein"]:
			assert(state[locus_id].has(species_name))
			for signature_variant in state[locus_id][species_name].keys():
				var signature: int = int(signature_variant)
				assert(signature >= 0 and signature <= 0xFFFF)
				assert(float(state[locus_id][species_name][signature_variant]) >= -1e-10)

static func _prune_zero_cohorts(cohorts: Dictionary) -> void:
	var to_remove: Array = []
	for signature_variant in cohorts.keys():
		if float(cohorts[signature_variant]) <= 1e-12: to_remove.append(signature_variant)
	for signature_variant in to_remove: cohorts.erase(signature_variant)

static func _empty_summary() -> Dictionary:
	return {
		"transcribed": 0.0,
		"translated": 0.0,
		"mrna_degraded": 0.0,
		"protein_degraded": 0.0,
		"atp_spent": 0.0,
		"tx_scale": 1.0,
		"translation_scale": 1.0,
		"regulation": {}
	}
