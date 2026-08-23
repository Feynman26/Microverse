extends RefCounted
class_name CatalyticLandscape

const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

const SIGNATURE_BITS: int = 16
const ACTIVE_MAX_DISTANCE: int = 4
const DISTANCE_DECAY: float = 0.70

# Sequence affinity is a pure function of two immutable 16-bit signatures.
# Cache only that molecular lookup, never abundance or cellular state. This
# removes repeated XOR/popcount/exp work while returning the exact first-call
# value on every later recognition event.
static var _affinity_cache: Dictionary = {}

static func hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & 0xFFFF
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

static func affinity(protein_signature: int, reaction_signature: int) -> float:
	var protein: int = protein_signature & 0xFFFF
	var reaction: int = reaction_signature & 0xFFFF
	var key: int = (protein << 16) | reaction
	if _affinity_cache.has(key):
		return float(_affinity_cache[key])
	var distance: int = hamming_distance(protein, reaction)
	var result: float = 0.0 if distance > ACTIVE_MAX_DISTANCE else exp(-DISTANCE_DECAY * float(distance))
	_affinity_cache[key] = result
	return result

static func is_active(protein_signature: int, reaction_signature: int) -> bool:
	return hamming_distance(protein_signature, reaction_signature) <= ACTIVE_MAX_DISTANCE

# M4 proxy retained only for landscape characterization and historical controls.
static func gene_activity(gene, reaction) -> float:
	return float(gene.promoter_strength()) * affinity(int(gene.protein_signature), int(reaction.signature)) * float(reaction.catalytic_ceiling)

static func genome_activity(genome, reaction) -> float:
	var result: float = 0.0
	for gene in genome.genes:
		result += gene_activity(gene, reaction)
	return result

# Production M5 catalysis reads the proteins that physically exist, including
# old-sequence cohorts inherited across a coding mutation. DNA is used for new
# expression upstream, not as a shortcut for current catalytic abundance.
static func proteome_activity(genome, expression_state: Dictionary, reaction, config) -> float:
	var result: float = 0.0
	for gene in genome.genes:
		var cohorts: Dictionary = ExpressionSystemScript.protein_cohorts_for_locus(expression_state, int(gene.locus_id))
		for signature_variant in cohorts.keys():
			var protein_signature: int = int(signature_variant)
			var abundance: float = maxf(0.0, float(cohorts[signature_variant])) / float(config.expression_reference_protein_count)
			result += abundance * affinity(protein_signature, int(reaction.signature)) * float(reaction.catalytic_ceiling)
	return result

static func strongest_protein_contribution(genome, expression_state: Dictionary, reaction, config) -> Dictionary:
	var best_locus: int = -1
	var best_signature: int = -1
	var best_activity: float = 0.0
	var best_distance: int = SIGNATURE_BITS + 1
	for gene in genome.genes:
		var cohorts: Dictionary = ExpressionSystemScript.protein_cohorts_for_locus(expression_state, int(gene.locus_id))
		for signature_variant in cohorts.keys():
			var protein_signature: int = int(signature_variant)
			var distance: int = hamming_distance(protein_signature, int(reaction.signature))
			var abundance: float = maxf(0.0, float(cohorts[signature_variant])) / float(config.expression_reference_protein_count)
			var activity: float = abundance * affinity(protein_signature, int(reaction.signature)) * float(reaction.catalytic_ceiling)
			if activity > best_activity:
				best_activity = activity
				best_locus = int(gene.locus_id)
				best_signature = protein_signature
				best_distance = distance
	return {"locus_id": best_locus, "protein_signature": best_signature, "activity": best_activity, "distance": best_distance}

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
