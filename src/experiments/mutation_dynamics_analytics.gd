extends RefCounted
class_name MutationDynamicsAnalytics

const GenomeScript = preload("res://src/genetics/genome.gd")
const DNAReplicationScript = preload("res://src/genetics/dna_replication.gd")

# Pure observational M10 analytics. These functions read simulation/event state
# and return measurements only; they never consume RNG, alter genomes/resources,
# assign fitness, or classify a mutation as beneficial/deleterious.

static func sample_population(sim) -> Dictionary:
	var lineages: Dictionary = {}
	var functional_genotypes: Dictionary = {}
	var point_rates: Array = []
	var structural_rates: Array = []
	var repair_values: Array = []
	var gene_counts: Array = []
	var replication_units: Array = []
	var departure_values: Array = []

	for cell in sim.cells:
		if not cell.alive or cell.genome == null:
			continue
		var repair: float = DNAReplicationScript.repair_activity(cell.expression_state, float(cell.volume), sim.config)
		var rates: Dictionary = _rates_from_repair(repair, sim.config)
		var key: String = str(cell.genome.fingerprint())
		if not lineages.has(key):
			lineages[key] = {
				"count": 0,
				"point_rate_sum": 0.0,
				"structural_rate_sum": 0.0,
				"repair_sum": 0.0,
				"gene_count_sum": 0.0,
				"replication_units_sum": 0.0,
				"ancestral_departure_sum": 0.0
			}
		var item: Dictionary = lineages[key]
		var departure: Dictionary = genome_departure_from_ancestor(cell.genome)
		item["count"] = int(item["count"]) + 1
		item["point_rate_sum"] = float(item["point_rate_sum"]) + float(rates["point_error_rate_per_gene"])
		item["structural_rate_sum"] = float(item["structural_rate_sum"]) + float(rates["structural_error_rate_per_genome"])
		item["repair_sum"] = float(item["repair_sum"]) + repair
		item["gene_count_sum"] = float(item["gene_count_sum"]) + float(cell.genome.gene_count())
		item["replication_units_sum"] = float(item["replication_units_sum"]) + float(cell.genome.replication_unit_count())
		item["ancestral_departure_sum"] = float(item["ancestral_departure_sum"]) + float(departure["departure_count"])
		lineages[key] = item
		var functional: String = functional_key(cell.genome)
		functional_genotypes[functional] = int(functional_genotypes.get(functional, 0)) + 1
		point_rates.append(float(rates["point_error_rate_per_gene"]))
		structural_rates.append(float(rates["structural_error_rate_per_genome"]))
		repair_values.append(repair)
		gene_counts.append(float(cell.genome.gene_count()))
		replication_units.append(float(cell.genome.replication_unit_count()))
		departure_values.append(float(departure["departure_count"]))

	for key_variant in lineages.keys():
		var key: String = String(key_variant)
		var item: Dictionary = lineages[key]
		var count: float = maxf(1.0, float(item["count"]))
		lineages[key] = {
			"count": int(item["count"]),
			"mean_point_error_rate_per_gene": float(item["point_rate_sum"]) / count,
			"mean_structural_error_rate_per_genome": float(item["structural_rate_sum"]) / count,
			"mean_repair_activity": float(item["repair_sum"]) / count,
			"mean_gene_count": float(item["gene_count_sum"]) / count,
			"mean_replication_units": float(item["replication_units_sum"]) / count,
			"mean_ancestral_departure": float(item["ancestral_departure_sum"]) / count
		}

	return {
		"population": int(sim.population_size()),
		"tick": int(sim.tick_index),
		"lineages": lineages,
		"genetic_genotype_count": lineages.size(),
		"functional_genotype_count": functional_genotypes.size(),
		"functional_genotype_frequencies": functional_genotypes,
		"point_error_rate_per_gene": _distribution(point_rates),
		"structural_error_rate_per_genome": _distribution(structural_rates),
		"repair_activity": _distribution(repair_values),
		"gene_count": _distribution(gene_counts),
		"replication_units": _distribution(replication_units),
		"ancestral_departure": _distribution(departure_values)
	}

# Functional identity is deliberately an analytical projection, never an
# authoritative genome identity. Stable locus IDs, neutral markers and current
# gene order are ancestry/history information. They are excluded because the
# present execution model does not use them to determine expression, binding or
# catalysis. Copy number is retained by keeping one tuple per physical locus.
static func functional_key(genome) -> String:
	assert(genome != null)
	var parts: Array[String] = []
	for gene in genome.genes:
		parts.append("%d:%d:%d:PC%d:RC%d" % [
			int(gene.promoter_code),
			int(gene.protein_signature),
			int(gene.regulatory_signature),
			int(gene.promoter_copy_number),
			int(gene.regulatory_copy_number)
		])
	parts.sort()
	return "|".join(parts)

static func summarize_event_log(event_log: Array) -> Dictionary:
	var birth_ids: Dictionary = {}
	var daughter_profiles: Dictionary = {}
	var mutations_by_cell: Dictionary = {}
	var attempts_by_cell: Dictionary = {}
	var mutation_types: Dictionary = {}
	var attempted_mutation_types: Dictionary = {}
	var blocked_mutation_types: Dictionary = {}
	var genome_size_deltas: Array = []
	var point_rates: Array = []
	var structural_rates: Array = []
	var repair_values: Array = []

	for event_variant in event_log:
		var event: Dictionary = event_variant
		var kind: String = String(event.get("kind", ""))
		if kind == "division":
			var profile: Dictionary = event.get("replication_profile", {})
			for daughter_variant in event.get("daughter_ids", []):
				var daughter_id: int = int(daughter_variant)
				daughter_profiles[daughter_id] = profile.duplicate(true)
				if profile.has("point_error_rate_per_gene"):
					point_rates.append(float(profile["point_error_rate_per_gene"]))
				if profile.has("structural_error_rate_per_genome"):
					structural_rates.append(float(profile["structural_error_rate_per_genome"]))
				if profile.has("mean_repair_activity"):
					repair_values.append(float(profile["mean_repair_activity"]))
		elif kind == "birth" and int(event.get("parent_id", -1)) >= 0:
			birth_ids[int(event.get("cell_id", -1))] = true
		elif kind == "mutation":
			var cell_id: int = int(event.get("cell_id", -1))
			var mutation_type: String = String(event.get("mutation_type", "unknown"))
			attempts_by_cell[cell_id] = int(attempts_by_cell.get(cell_id, 0)) + 1
			attempted_mutation_types[mutation_type] = int(attempted_mutation_types.get(mutation_type, 0)) + 1
			if mutation_type.ends_with("_blocked"):
				blocked_mutation_types[mutation_type] = int(blocked_mutation_types.get(mutation_type, 0)) + 1
				continue
			mutations_by_cell[cell_id] = int(mutations_by_cell.get(cell_id, 0)) + 1
			mutation_types[mutation_type] = int(mutation_types.get(mutation_type, 0)) + 1
			if event.has("parent_genome_size") and event.has("resulting_genome_size"):
				genome_size_deltas.append(float(event["resulting_genome_size"]) - float(event["parent_genome_size"]))

	var realized_per_birth: Array = []
	var attempted_per_birth: Array = []
	for birth_id_variant in birth_ids.keys():
		var birth_id: int = int(birth_id_variant)
		realized_per_birth.append(float(mutations_by_cell.get(birth_id, 0)))
		attempted_per_birth.append(float(attempts_by_cell.get(birth_id, 0)))

	return {
		"births": birth_ids.size(),
		"mutations": _sum_counts(mutation_types),
		"mutation_attempts": _sum_counts(attempted_mutation_types),
		"blocked_mutation_attempts": _sum_counts(blocked_mutation_types),
		"mutation_types": mutation_types,
		"attempted_mutation_types": attempted_mutation_types,
		"blocked_mutation_types": blocked_mutation_types,
		"realized_mutations_per_birth": _distribution(realized_per_birth),
		"attempted_mutations_per_birth": _distribution(attempted_per_birth),
		"expected_point_error_rate_per_gene": _distribution(point_rates),
		"expected_structural_error_rate_per_genome": _distribution(structural_rates),
		"replication_repair_activity": _distribution(repair_values),
		"genome_size_delta_per_mutation_event": _distribution(genome_size_deltas),
		"births_with_replication_profile": daughter_profiles.size()
	}

static func genome_departure_from_ancestor(genome) -> Dictionary:
	assert(genome != null)
	var ancestor = GenomeScript.create_ancestor()
	var missing_loci: Array = []
	var extra_loci: Array = []
	var coding_changed: Array = []
	var promoter_changed: Array = []
	var regulatory_changed: Array = []
	var cis_changed: Array = []

	for ancestral_gene in ancestor.genes:
		var current = genome.get_gene_by_locus(int(ancestral_gene.locus_id))
		if current == null:
			missing_loci.append(int(ancestral_gene.locus_id))
			continue
		if int(current.protein_signature) != int(ancestral_gene.protein_signature):
			coding_changed.append(int(ancestral_gene.locus_id))
		if int(current.promoter_code) != int(ancestral_gene.promoter_code):
			promoter_changed.append(int(ancestral_gene.locus_id))
		if int(current.regulatory_signature) != int(ancestral_gene.regulatory_signature):
			regulatory_changed.append(int(ancestral_gene.locus_id))
		if int(current.promoter_copy_number) != 1 or int(current.regulatory_copy_number) != 1:
			cis_changed.append(int(ancestral_gene.locus_id))

	for current_gene in genome.genes:
		if ancestor.get_gene_by_locus(int(current_gene.locus_id)) == null:
			extra_loci.append(int(current_gene.locus_id))

	var departure_count: int = (
		missing_loci.size() + extra_loci.size() + coding_changed.size()
		+ promoter_changed.size() + regulatory_changed.size() + cis_changed.size()
	)
	return {
		"departure_count": departure_count,
		"missing_ancestral_loci": missing_loci,
		"extra_loci": extra_loci,
		"coding_signature_changes": coding_changed,
		"promoter_changes": promoter_changed,
		"regulatory_signature_changes": regulatory_changed,
		"cis_copy_changes": cis_changed,
		"note": "molecular departure proxy only; no beneficial/deleterious or fitness interpretation"
	}

# Descriptive association around a known environmental shift. This measures
# whether higher-copy-error lineages gained or lost abundance; it does not claim
# the environmental change caused that shift without a paired causal fork.
static func environment_shift_association(before: Dictionary, after: Dictionary) -> Dictionary:
	var before_mean: float = _weighted_lineage_rate(before.get("lineages", {}))
	var after_mean: float = _weighted_lineage_rate(after.get("lineages", {}))
	var before_high: float = _high_error_abundance_fraction(before.get("lineages", {}), before_mean)
	var after_high: float = _high_error_abundance_fraction(after.get("lineages", {}), before_mean)
	return {
		"before_tick": int(before.get("tick", -1)),
		"after_tick": int(after.get("tick", -1)),
		"weighted_point_error_rate_before": before_mean,
		"weighted_point_error_rate_after": after_mean,
		"weighted_point_error_rate_change": after_mean - before_mean,
		"above_baseline_error_abundance_before": before_high,
		"above_baseline_error_abundance_after": after_high,
		"above_baseline_error_abundance_change": after_high - before_high,
		"interpretation": "descriptive association; use M9 paired fork/revert for causal attribution"
	}

static func _rates_from_repair(repair: float, config) -> Dictionary:
	var suppression: float = 1.0 + float(config.dna_repair_fidelity_gain) * maxf(0.0, repair)
	return {
		"point_error_rate_per_gene": clampf(
			float(config.baseline_point_error_rate_per_gene) / suppression,
			float(config.minimum_point_error_rate_per_gene),
			float(config.baseline_point_error_rate_per_gene)
		),
		"structural_error_rate_per_genome": clampf(
			float(config.baseline_structural_error_rate_per_genome) / suppression,
			float(config.minimum_structural_error_rate_per_genome),
			float(config.baseline_structural_error_rate_per_genome)
		)
	}

static func _weighted_lineage_rate(lineages: Dictionary) -> float:
	var weighted: float = 0.0
	var population: float = 0.0
	for item_variant in lineages.values():
		var item: Dictionary = item_variant
		var count: float = float(item.get("count", 0))
		weighted += count * float(item.get("mean_point_error_rate_per_gene", 0.0))
		population += count
	return weighted / maxf(1.0, population)

static func _high_error_abundance_fraction(lineages: Dictionary, threshold: float) -> float:
	var high: float = 0.0
	var population: float = 0.0
	for item_variant in lineages.values():
		var item: Dictionary = item_variant
		var count: float = float(item.get("count", 0))
		population += count
		if float(item.get("mean_point_error_rate_per_gene", 0.0)) > threshold + 1e-15:
			high += count
	return high / maxf(1.0, population)

static func _distribution(values: Array) -> Dictionary:
	if values.is_empty():
		return {"count": 0, "min": 0.0, "mean": 0.0, "max": 0.0}
	var minimum: float = INF
	var maximum: float = -INF
	var total: float = 0.0
	for value_variant in values:
		var value: float = float(value_variant)
		minimum = minf(minimum, value)
		maximum = maxf(maximum, value)
		total += value
	return {"count": values.size(), "min": minimum, "mean": total / float(values.size()), "max": maximum}

static func _sum_counts(counts: Dictionary) -> int:
	var total: int = 0
	for value_variant in counts.values():
		total += int(value_variant)
	return total
