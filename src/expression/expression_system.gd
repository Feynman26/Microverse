extends RefCounted
class_name ExpressionSystem

const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

# M5-A state is keyed by stable locus ID. mRNA/protein are digital molecular
# amounts. Birth/death proposals are Poisson-distributed using the universe RNG;
# accepted synthesis is then allocated simultaneously across loci under shared
# ATP/NUC/AA scarcity, preventing gene-array order from becoming expression
# priority.

static func create_equilibrium_state(genome, config) -> Dictionary:
	var state: Dictionary = {}
	var tx_max: float = float(config.transcription_max_events_per_min)
	var mrna_decay: float = float(config.mrna_decay_rate_per_min)
	var translation: float = float(config.translation_events_per_mrna_per_min)
	var protein_decay: float = float(config.protein_decay_rate_per_min)
	for gene in genome.genes:
		var promoter: float = float(gene.promoter_strength())
		var mrna: float = tx_max * promoter / mrna_decay
		var protein: float = translation * mrna / protein_decay
		state[int(gene.locus_id)] = {"mrna": mrna, "protein": protein}
	return state

static func step(state: Dictionary, genome, pools: Dictionary, dt: float, rng, config) -> Dictionary:
	assert(dt >= 0.0)
	if dt <= 0.0:
		return _empty_summary()
	_validate_state(state, genome)

	var snapshot: Dictionary = state.duplicate(true)
	var proposals: Dictionary = {}
	var total_tx_events: float = 0.0
	var total_translation_events: float = 0.0

	# Propose stochastic molecular events from one immutable expression snapshot.
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		var prior: Dictionary = snapshot[locus_id]
		var mrna: float = maxf(0.0, float(prior["mrna"]))
		var protein: float = maxf(0.0, float(prior["protein"]))
		var tx_lambda: float = float(config.transcription_max_events_per_min) * float(gene.promoter_strength()) * dt
		var tl_lambda: float = float(config.translation_events_per_mrna_per_min) * mrna * dt
		var mrna_decay_lambda: float = float(config.mrna_decay_rate_per_min) * mrna * dt
		var protein_decay_lambda: float = float(config.protein_decay_rate_per_min) * protein * dt

		var tx_events: float = float(rng.poisson(tx_lambda))
		var translation_events: float = float(rng.poisson(tl_lambda))
		var mrna_decay_events: float = minf(mrna, float(rng.poisson(mrna_decay_lambda)))
		var protein_decay_events: float = minf(protein, float(rng.poisson(protein_decay_lambda)))
		proposals[locus_id] = {
			"tx": tx_events,
			"translation": translation_events,
			"mrna_decay": mrna_decay_events,
			"protein_decay": protein_decay_events
		}
		total_tx_events += tx_events
		total_translation_events += translation_events

	# Shared synthesis resources are allocated with one scale derived from the
	# pre-synthesis pool snapshot. This is the same fairness principle used for
	# cellular nutrient competition and metabolic substrate allocation.
	var tx_atp_demand: float = total_tx_events * float(config.transcription_atp_cost_per_event)
	var tl_atp_demand: float = total_translation_events * float(config.translation_atp_cost_per_event)
	var total_atp_demand: float = tx_atp_demand + tl_atp_demand
	var available_atp: float = maxf(0.0, float(pools.get("ATP", 0.0)))
	var atp_scale: float = 1.0 if total_atp_demand <= available_atp or total_atp_demand <= 0.0 else available_atp / total_atp_demand

	var nuc_demand: float = total_tx_events * float(config.transcription_nuc_cost_per_event)
	var aa_demand: float = total_translation_events * float(config.translation_aa_cost_per_event)
	var available_nuc: float = maxf(0.0, float(pools.get("NUC", 0.0)))
	var available_aa: float = maxf(0.0, float(pools.get("AA", 0.0)))
	var nuc_scale: float = 1.0 if nuc_demand <= available_nuc or nuc_demand <= 0.0 else available_nuc / nuc_demand
	var aa_scale: float = 1.0 if aa_demand <= available_aa or aa_demand <= 0.0 else available_aa / aa_demand
	var tx_scale: float = minf(atp_scale, nuc_scale)
	var translation_scale: float = minf(atp_scale, aa_scale)

	var accepted_tx: float = 0.0
	var accepted_translation: float = 0.0
	var degraded_mrna: float = 0.0
	var degraded_protein: float = 0.0
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		var prior: Dictionary = snapshot[locus_id]
		var proposal: Dictionary = proposals[locus_id]
		var tx_events: float = float(proposal["tx"]) * tx_scale
		var translation_events: float = float(proposal["translation"]) * translation_scale
		var mrna_decay_events: float = float(proposal["mrna_decay"])
		var protein_decay_events: float = float(proposal["protein_decay"])
		state[locus_id] = {
			"mrna": maxf(0.0, float(prior["mrna"]) + tx_events - mrna_decay_events),
			"protein": maxf(0.0, float(prior["protein"]) + translation_events - protein_decay_events)
		}
		accepted_tx += tx_events
		accepted_translation += translation_events
		degraded_mrna += mrna_decay_events
		degraded_protein += protein_decay_events

	# Synthesis transfers material from precursor pools into expression storage;
	# degradation returns it (compressed complete recycling). Energy expenditure
	# is ATP -> ADP and therefore remains inside the adenylate currency pool.
	var synthesis_atp: float = (
		accepted_tx * float(config.transcription_atp_cost_per_event)
		+ accepted_translation * float(config.translation_atp_cost_per_event)
	)
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
		"translation_scale": translation_scale
	}

static func normalized_protein(state: Dictionary, locus_id: int, config) -> float:
	if not state.has(locus_id):
		return 0.0
	return maxf(0.0, float(state[locus_id]["protein"])) / float(config.expression_reference_protein_count)

static func total_mrna(state: Dictionary) -> float:
	var result: float = 0.0
	for item in state.values():
		result += float(item["mrna"])
	return result

static func total_protein(state: Dictionary) -> float:
	var result: float = 0.0
	for item in state.values():
		result += float(item["protein"])
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
		var item: Dictionary = state[locus_id]
		var mrna: float = maxf(0.0, float(item["mrna"]))
		var protein: float = maxf(0.0, float(item["protein"]))
		var mrna_ratio: float = _noisy_partition_ratio(ratio, mrna, noise_scale, rng)
		var protein_ratio: float = _noisy_partition_ratio(ratio, protein, noise_scale, rng)
		first[locus_id] = {"mrna": mrna * mrna_ratio, "protein": protein * protein_ratio}
		second[locus_id] = {"mrna": mrna * (1.0 - mrna_ratio), "protein": protein * (1.0 - protein_ratio)}
	return [first, second]

static func _noisy_partition_ratio(base_ratio: float, amount: float, scale: float, rng) -> float:
	if amount <= 0.0 or scale <= 0.0:
		return base_ratio
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
		result += float(item["mrna"]) * float(i + 3) * 0.013
		result += float(item["protein"]) * float(i + 5) * 0.0017
	return result

static func _validate_state(state: Dictionary, genome) -> void:
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		assert(state.has(locus_id), "Missing expression state for locus %d" % locus_id)
		assert(float(state[locus_id]["mrna"]) >= -1e-10)
		assert(float(state[locus_id]["protein"]) >= -1e-10)

static func _empty_summary() -> Dictionary:
	return {
		"transcribed": 0.0,
		"translated": 0.0,
		"mrna_degraded": 0.0,
		"protein_degraded": 0.0,
		"atp_spent": 0.0,
		"tx_scale": 1.0,
		"translation_scale": 1.0
	}
