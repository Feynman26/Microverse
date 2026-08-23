extends RefCounted
class_name DNAReplication

const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")

# Deliberately outside the ancestral active radius. The closest ancestral
# protein is locus 12 (0xC136) at Hamming distance 5; one ordinary coding bit
# flip can therefore enter the configured distance-4 repair landscape without
# gifting the ancestor a repair phenotype.
const REPAIR_TARGET_SIGNATURE: int = 0xD15A

static func step(cell, dt: float, config) -> Dictionary:
	assert(dt >= 0.0)
	assert(cell.genome != null)
	if not bool(config.evolvable_replication_enabled):
		return _summary(cell, 0.0, 0.0, 0.0, 0.0, true, config)
	if dt <= 0.0:
		return _summary(cell, 0.0, 0.0, 0.0, 0.0, replication_complete(cell), config)
	if replication_complete(cell):
		return _summary(cell, 0.0, 0.0, 0.0, 0.0, true, config)
	if float(cell.volume) < float(config.genome_replication_initiation_volume):
		return _summary(cell, 0.0, 0.0, 0.0, 0.0, false, config)

	var gene_count: float = float(maxi(1, cell.genome.gene_count()))
	var remaining_gene_equivalents: float = maxf(0.0, gene_count - float(cell.replication_gene_equivalents_copied))
	if remaining_gene_equivalents <= 1e-12:
		cell.replication_gene_equivalents_copied = gene_count
		cell.replication_progress = 1.0
		return _summary(cell, 0.0, 0.0, 0.0, 0.0, true, config)

	var repair: float = repair_activity(cell.expression_state, float(cell.volume), config)
	var requested_gene_equivalents: float = minf(
		remaining_gene_equivalents,
		float(config.genome_replication_gene_copy_rate_per_min) * dt
	)
	var atp_cost_per_gene: float = (
		float(config.genome_replication_atp_cost_per_gene)
		+ repair * float(config.dna_repair_atp_cost_per_gene_activity)
	)
	var nuc_cost_per_gene: float = float(config.genome_replication_nuc_cost_per_gene)
	var requested_atp: float = requested_gene_equivalents * atp_cost_per_gene
	var requested_nuc: float = requested_gene_equivalents * nuc_cost_per_gene
	var available_atp: float = maxf(0.0, cell.pool("ATP"))
	var available_nuc: float = maxf(0.0, cell.pool("NUC"))
	var atp_scale: float = 1.0 if requested_atp <= available_atp or requested_atp <= 0.0 else available_atp / requested_atp
	var nuc_scale: float = 1.0 if requested_nuc <= available_nuc or requested_nuc <= 0.0 else available_nuc / requested_nuc
	var resource_scale: float = minf(atp_scale, nuc_scale)
	var copied: float = requested_gene_equivalents * resource_scale
	var atp_spent: float = copied * atp_cost_per_gene
	var nuc_spent: float = copied * nuc_cost_per_gene

	var actual_atp: float = MetabolicSolverScript.spend_atp(cell.metabolites, atp_spent)
	assert(absf(actual_atp - atp_spent) <= 1e-9, "Replication ATP pre-scaling failed")
	cell.metabolites["NUC"] = maxf(0.0, cell.pool("NUC") - nuc_spent)
	cell.replication_gene_equivalents_copied += copied
	cell.replication_atp_spent += atp_spent
	cell.replication_nuc_spent += nuc_spent
	cell.replication_repair_activity_integral += copied * repair
	cell.replication_progress = clampf(cell.replication_gene_equivalents_copied / gene_count, 0.0, 1.0)
	if gene_count - cell.replication_gene_equivalents_copied <= 1e-10:
		cell.replication_gene_equivalents_copied = gene_count
		cell.replication_progress = 1.0

	return _summary(cell, copied, atp_spent, nuc_spent, repair, replication_complete(cell), config)

static func repair_activity(expression_state: Dictionary, cell_volume: float, config) -> float:
	assert(cell_volume > 0.0)
	var reference_abundance: float = (
		float(config.expression_reference_protein_count)
		* cell_volume
		/ maxf(1e-12, float(config.ancestor_volume))
	)
	var activity: float = 0.0
	var loci: Array = expression_state.keys()
	loci.sort()
	for locus_variant in loci:
		var cohorts: Dictionary = expression_state[locus_variant]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			var distance: int = hamming_distance(signature, REPAIR_TARGET_SIGNATURE)
			if distance > int(config.dna_repair_max_distance):
				continue
			var abundance: float = maxf(0.0, float(cohorts[signature_variant])) / maxf(1e-12, reference_abundance)
			activity += abundance * exp(-float(config.dna_repair_distance_decay) * float(distance))
	return maxf(0.0, activity)

static func mutation_profile(cell, config) -> Dictionary:
	var mean_repair: float = mean_repair_activity(cell)
	var suppression: float = 1.0 + float(config.dna_repair_fidelity_gain) * mean_repair
	var point_probability: float = clampf(
		float(config.baseline_point_error_rate_per_gene) / suppression,
		float(config.minimum_point_error_rate_per_gene),
		float(config.baseline_point_error_rate_per_gene)
	)
	var structural_probability: float = clampf(
		float(config.baseline_structural_error_rate_per_genome) / suppression,
		float(config.minimum_structural_error_rate_per_genome),
		float(config.baseline_structural_error_rate_per_genome)
	)
	return {
		"mean_repair_activity": mean_repair,
		"point_error_rate_per_gene": point_probability,
		"structural_error_rate_per_genome": structural_probability,
		"replication_gene_equivalents": float(cell.replication_gene_equivalents_copied),
		"replication_atp_spent": float(cell.replication_atp_spent),
		"replication_nuc_spent": float(cell.replication_nuc_spent)
	}

static func mean_repair_activity(cell) -> float:
	var copied: float = maxf(0.0, float(cell.replication_gene_equivalents_copied))
	if copied <= 1e-12:
		return 0.0
	return maxf(0.0, float(cell.replication_repair_activity_integral) / copied)

static func replication_complete(cell) -> bool:
	return float(cell.replication_progress) >= 1.0 - 1e-10

static func hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & 0xFFFF
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

static func _summary(cell, copied: float, atp_spent: float, nuc_spent: float, repair: float, complete: bool, config) -> Dictionary:
	var gene_count: int = maxi(1, cell.genome.gene_count())
	return {
		"enabled": bool(config.evolvable_replication_enabled),
		"started": float(cell.volume) >= float(config.genome_replication_initiation_volume),
		"complete": complete,
		"gene_count": gene_count,
		"progress": float(cell.replication_progress),
		"gene_equivalents_copied": float(cell.replication_gene_equivalents_copied),
		"copied_this_tick": copied,
		"atp_spent_this_tick": atp_spent,
		"nuc_spent_this_tick": nuc_spent,
		"repair_activity_this_tick": repair,
		"mean_repair_activity": mean_repair_activity(cell),
		"cumulative_atp_spent": float(cell.replication_atp_spent),
		"cumulative_nuc_spent": float(cell.replication_nuc_spent)
	}