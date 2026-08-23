extends RefCounted
class_name MutationEngine

const GeneScript = preload("res://src/genetics/gene.gd")

const STRUCTURAL_OPERATIONS: Array[String] = [
	"gene_duplication",
	"gene_deletion",
	"segment_inversion",
	"local_rearrangement"
]

# Historical M3 mutation harness. It remains intentionally unchanged so the
# original Bernoulli-rate tests continue to characterize that milestone. The
# production M10 simulation uses mutate_replicated_copy() instead.
func mutate_copy(parent_genome, rng, config) -> Dictionary:
	var child_genome = parent_genome.deep_copy()
	var events: Array = []
	if not bool(config.mutation_enabled):
		return {"genome": child_genome, "events": events}

	for gene in child_genome.genes:
		if float(rng.randf()) < float(config.promoter_mutation_rate_per_gene):
			events.append(_mutate_promoter(gene, rng, config))
		if float(rng.randf()) < float(config.signature_mutation_rate_per_gene):
			events.append(_mutate_signature(gene, rng, config))
		if float(rng.randf()) < float(config.regulatory_signature_mutation_rate_per_gene):
			events.append(_mutate_regulatory_signature(gene, rng, config))
		if float(rng.randf()) < float(config.neutral_marker_mutation_rate_per_gene):
			events.append(_mutate_neutral_marker(gene, rng))

	child_genome.validate()
	return {"genome": child_genome, "events": events}

# M10 production heredity. Point-error probability is supplied by the completed
# physical replication cycle and therefore already incorporates sequence-derived
# repair investment and its resource cost. One error draw is made per copied
# gene; if an error occurs, its molecular target is chosen without a beneficial
# or deleterious classification. Structural copying errors are genome-level and
# use the same replication-derived profile.
func mutate_replicated_copy(parent_genome, rng, config, replication_profile: Dictionary) -> Dictionary:
	var child_genome = parent_genome.deep_copy()
	var events: Array = []
	if not bool(config.mutation_enabled):
		return {"genome": child_genome, "events": events}

	var point_probability: float = clampf(float(replication_profile.get("point_error_rate_per_gene", 0.0)), 0.0, 1.0)
	var structural_probability: float = clampf(float(replication_profile.get("structural_error_rate_per_genome", 0.0)), 0.0, 1.0)
	var copied_genes: Array = child_genome.genes.duplicate()
	for gene in copied_genes:
		if float(rng.randf()) >= point_probability:
			continue
		var channel: int = int(rng.randi_range(0, 3))
		var event: Dictionary
		match channel:
			0: event = _mutate_promoter(gene, rng, config)
			1: event = _mutate_signature(gene, rng, config)
			2: event = _mutate_regulatory_signature(gene, rng, config)
			_: event = _mutate_neutral_marker(gene, rng)
		event["replication_derived"] = true
		event["error_probability"] = point_probability
		events.append(event)

	if float(rng.randf()) < structural_probability:
		var operation: String = STRUCTURAL_OPERATIONS[int(rng.randi_range(0, STRUCTURAL_OPERATIONS.size() - 1))]
		var structural_event: Dictionary = apply_structural_mutation(child_genome, operation, rng)
		structural_event["replication_derived"] = true
		structural_event["error_probability"] = structural_probability
		events.append(structural_event)

	child_genome.validate()
	return {"genome": child_genome, "events": events}

# Public controlled operator used by M10 validation. It edits only the supplied
# genome and never assigns biological value to the change.
func apply_structural_mutation(genome, operation: String, rng) -> Dictionary:
	assert(genome != null)
	assert(STRUCTURAL_OPERATIONS.has(operation), "Unknown structural mutation operation: %s" % operation)
	match operation:
		"gene_duplication": return _duplicate_gene(genome, rng)
		"gene_deletion": return _delete_gene(genome, rng)
		"segment_inversion": return _invert_segment(genome, rng)
		"local_rearrangement": return _rearrange_local(genome, rng)
	assert(false)
	return {}

func _duplicate_gene(genome, rng) -> Dictionary:
	assert(genome.genes.size() >= 1)
	var source_index: int = int(rng.randi_range(0, genome.genes.size() - 1))
	var source = genome.genes[source_index]
	var duplicate = source.deep_copy()
	duplicate.locus_id = _draw_unique_locus_id(genome, rng)
	genome.genes.insert(source_index + 1, duplicate)
	return {
		"mutation_type": "gene_duplication",
		"source_locus_id": int(source.locus_id),
		"new_locus_id": int(duplicate.locus_id),
		"insert_index": source_index + 1,
		"copied_gene": duplicate.canonical_key()
	}

func _delete_gene(genome, rng) -> Dictionary:
	if genome.genes.size() <= 1:
		# A zero-gene genome is outside the current executable chemistry model.
		# Preserve a valid genome while still returning an explicit no-op record.
		return {
			"mutation_type": "gene_deletion_blocked",
			"reason": "minimum_one_gene",
			"locus_id": int(genome.genes[0].locus_id)
		}
	var index: int = int(rng.randi_range(0, genome.genes.size() - 1))
	var removed = genome.genes[index]
	var removed_key: String = removed.canonical_key()
	var removed_locus: int = int(removed.locus_id)
	genome.genes.remove_at(index)
	return {
		"mutation_type": "gene_deletion",
		"locus_id": removed_locus,
		"former_index": index,
		"deleted_gene": removed_key
	}

func _invert_segment(genome, rng) -> Dictionary:
	if genome.genes.size() <= 1:
		return {"mutation_type": "segment_inversion_blocked", "reason": "single_gene_genome"}
	var first: int = int(rng.randi_range(0, genome.genes.size() - 1))
	var second: int = int(rng.randi_range(0, genome.genes.size() - 2))
	if second >= first:
		second += 1
	var start_index: int = mini(first, second)
	var end_index: int = maxi(first, second)
	var before_loci: Array = []
	for i in range(start_index, end_index + 1):
		before_loci.append(int(genome.genes[i].locus_id))
	var reversed: Array = genome.genes.slice(start_index, end_index + 1)
	reversed.reverse()
	for offset in range(reversed.size()):
		genome.genes[start_index + offset] = reversed[offset]
	var after_loci: Array = []
	for i in range(start_index, end_index + 1):
		after_loci.append(int(genome.genes[i].locus_id))
	return {
		"mutation_type": "segment_inversion",
		"start_index": start_index,
		"end_index": end_index,
		"before_loci": before_loci,
		"after_loci": after_loci
	}

func _rearrange_local(genome, rng) -> Dictionary:
	if genome.genes.size() <= 1:
		return {"mutation_type": "local_rearrangement_blocked", "reason": "single_gene_genome"}
	var first_index: int = int(rng.randi_range(0, genome.genes.size() - 2))
	var second_index: int = first_index + 1
	var first_locus: int = int(genome.genes[first_index].locus_id)
	var second_locus: int = int(genome.genes[second_index].locus_id)
	var temp = genome.genes[first_index]
	genome.genes[first_index] = genome.genes[second_index]
	genome.genes[second_index] = temp
	return {
		"mutation_type": "local_rearrangement",
		"first_index": first_index,
		"second_index": second_index,
		"before_loci": [first_locus, second_locus],
		"after_loci": [second_locus, first_locus]
	}

func _draw_unique_locus_id(genome, rng) -> int:
	var used: Dictionary = {}
	for gene in genome.genes:
		used[int(gene.locus_id)] = true
	var candidate: int = int(rng.randi_range(1, 2147483646))
	while used.has(candidate):
		candidate = int(rng.randi_range(1, 2147483646))
	return candidate

func _mutate_promoter(gene, rng, config) -> Dictionary:
	var old_value: int = int(gene.promoter_code)
	var max_step: int = int(config.promoter_mutation_step_max)
	var delta: int = int(rng.randi_range(-max_step, max_step))
	if delta == 0: delta = 1
	var new_value: int = clampi(old_value + delta, GeneScript.PROMOTER_CODE_MIN, GeneScript.PROMOTER_CODE_MAX)
	if new_value == old_value:
		new_value = old_value - 1 if old_value >= GeneScript.PROMOTER_CODE_MAX else old_value + 1
	gene.promoter_code = new_value
	return {"mutation_type": "promoter_code", "locus_id": int(gene.locus_id), "old_value": old_value, "new_value": new_value}

func _mutate_signature(gene, rng, config) -> Dictionary:
	var result: Dictionary = _flip_signature(int(gene.protein_signature), rng, config)
	gene.protein_signature = int(result["new_value"])
	result["mutation_type"] = "protein_signature_bit_flip"
	result["locus_id"] = int(gene.locus_id)
	return result

func _mutate_regulatory_signature(gene, rng, config) -> Dictionary:
	var result: Dictionary = _flip_signature(int(gene.regulatory_signature), rng, config)
	gene.regulatory_signature = int(result["new_value"])
	result["mutation_type"] = "regulatory_signature_bit_flip"
	result["locus_id"] = int(gene.locus_id)
	return result

func _flip_signature(old_value: int, rng, config) -> Dictionary:
	var bit_index: int = int(rng.randi_range(0, int(config.protein_signature_bits) - 1))
	var new_value: int = (old_value ^ (1 << bit_index)) & GeneScript.SIGNATURE_MASK
	return {"bit_index": bit_index, "old_value": old_value, "new_value": new_value}

func _mutate_neutral_marker(gene, rng) -> Dictionary:
	var old_value: int = int(gene.neutral_marker)
	var new_value: int = int(rng.randi_range(0, GeneScript.NEUTRAL_MARKER_MASK))
	if new_value == old_value: new_value = (old_value + 1) & GeneScript.NEUTRAL_MARKER_MASK
	gene.neutral_marker = new_value
	return {"mutation_type": "neutral_marker", "locus_id": int(gene.locus_id), "old_value": old_value, "new_value": new_value}
