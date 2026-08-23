extends RefCounted
class_name MutationEngine

const GeneScript = preload("res://src/genetics/gene.gd")

# M3 mutations operate only on inherited molecular information. They do not
# alter basal physiology directly. M4/M5 give promoter/signature mutations
# functional consequences through proteins, reactions and regulation.

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
		if float(rng.randf()) < float(config.neutral_marker_mutation_rate_per_gene):
			events.append(_mutate_neutral_marker(gene, rng))

	child_genome.validate()
	return {"genome": child_genome, "events": events}

func _mutate_promoter(gene, rng, config) -> Dictionary:
	var old_value: int = int(gene.promoter_code)
	var max_step: int = int(config.promoter_mutation_step_max)
	var delta: int = int(rng.randi_range(-max_step, max_step))
	if delta == 0:
		delta = 1
	var new_value: int = clampi(old_value + delta, GeneScript.PROMOTER_CODE_MIN, GeneScript.PROMOTER_CODE_MAX)
	if new_value == old_value:
		new_value = old_value - 1 if old_value >= GeneScript.PROMOTER_CODE_MAX else old_value + 1
	gene.promoter_code = new_value
	return {
		"mutation_type": "promoter_code",
		"locus_id": int(gene.locus_id),
		"old_value": old_value,
		"new_value": new_value
	}

func _mutate_signature(gene, rng, config) -> Dictionary:
	var old_value: int = int(gene.protein_signature)
	var bit_count: int = int(config.protein_signature_bits)
	var bit_index: int = int(rng.randi_range(0, bit_count - 1))
	var new_value: int = (old_value ^ (1 << bit_index)) & GeneScript.SIGNATURE_MASK
	gene.protein_signature = new_value
	return {
		"mutation_type": "protein_signature_bit_flip",
		"locus_id": int(gene.locus_id),
		"bit_index": bit_index,
		"old_value": old_value,
		"new_value": new_value
	}

func _mutate_neutral_marker(gene, rng) -> Dictionary:
	var old_value: int = int(gene.neutral_marker)
	var new_value: int = int(rng.randi_range(0, GeneScript.NEUTRAL_MARKER_MASK))
	if new_value == old_value:
		new_value = (old_value + 1) & GeneScript.NEUTRAL_MARKER_MASK
	gene.neutral_marker = new_value
	return {
		"mutation_type": "neutral_marker",
		"locus_id": int(gene.locus_id),
		"old_value": old_value,
		"new_value": new_value
	}
