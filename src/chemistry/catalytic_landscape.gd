extends RefCounted
class_name CatalyticLandscape

const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

const SIGNATURE_BITS: int = 16
const ACTIVE_MAX_DISTANCE: int = 4
const DISTANCE_DECAY: float = 0.70

static func hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & 0xFFFF
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

static func affinity(protein_signature: int, reaction_signature: int) -> float:
	var distance: int = hamming_distance(protein_signature, reaction_signature)
	if distance > ACTIVE_MAX_DISTANCE:
		return 0.0
	return exp(-DISTANCE_DECAY * float(distance))

static func is_active(protein_signature: int, reaction_signature: int) -> bool:
	return hamming_distance(protein_signature, reaction_signature) <= ACTIVE_MAX_DISTANCE

# M4 proxy retained only for landscape characterization and controlled backward
# comparisons. Production M5 physiology must use proteome_activity().
static func gene_activity(gene, reaction) -> float:
	return float(gene.promoter_strength()) * affinity(int(gene.protein_signature), int(reaction.signature)) * float(reaction.catalytic_ceiling)

static func genome_activity(genome, reaction) -> float:
	var result: float = 0.0
	for gene in genome.genes:
		result += gene_activity(gene, reaction)
	return result

# M5 catalytic capacity is a consequence of molecules that actually exist in
# this cell. Promoter code no longer enters metabolism directly: it affects
# transcription upstream, and protein abundance then enters catalysis here.
static func proteome_activity(genome, expression_state: Dictionary, reaction, config) -> float:
	var result: float = 0.0
	for gene in genome.genes:
		var abundance: float = ExpressionSystemScript.normalized_protein(expression_state, int(gene.locus_id), config)
		result += abundance * affinity(int(gene.protein_signature), int(reaction.signature)) * float(reaction.catalytic_ceiling)
	return result

static func strongest_protein_contribution(genome, expression_state: Dictionary, reaction, config) -> Dictionary:
	var best_locus: int = -1
	var best_activity: float = 0.0
	var best_distance: int = SIGNATURE_BITS + 1
	for gene in genome.genes:
		var distance: int = hamming_distance(int(gene.protein_signature), int(reaction.signature))
		var abundance: float = ExpressionSystemScript.normalized_protein(expression_state, int(gene.locus_id), config)
		var activity: float = abundance * affinity(int(gene.protein_signature), int(reaction.signature)) * float(reaction.catalytic_ceiling)
		if activity > best_activity:
			best_activity = activity
			best_locus = int(gene.locus_id)
			best_distance = distance
	return {"locus_id": best_locus, "activity": best_activity, "distance": best_distance}

static func strongest_gene_activity(genome, reaction) -> Dictionary:
	var best_locus: int = -1
	var best_activity: float = 0.0
	var best_distance: int = SIGNATURE_BITS + 1
	for gene in genome.genes:
		var distance: int = hamming_distance(int(gene.protein_signature), int(reaction.signature))
		var activity: float = gene_activity(gene, reaction)
		if activity > best_activity:
			best_activity = activity
			best_locus = int(gene.locus_id)
			best_distance = distance
	return {"locus_id": best_locus, "activity": best_activity, "distance": best_distance}

static func active_reaction_count_for_signature(protein_signature: int, reactions: Array) -> int:
	var result: int = 0
	for reaction in reactions:
		if is_active(protein_signature, int(reaction.signature)):
			result += 1
	return result
