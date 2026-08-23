extends RefCounted
class_name ExpressionState

# Authoritative M5 expression state. Keys are stable locus IDs; values are
# non-negative normalized molecular abundance units. The state is intentionally
# separate from Genome: genotype is inherited information, while mRNA/protein
# abundance is a transient phenotype that can differ between clonal cells.

var mrna: Dictionary = {}
var proteins: Dictionary = {}

func _init(p_mrna: Dictionary = {}, p_proteins: Dictionary = {}) -> void:
	mrna = p_mrna.duplicate(true)
	proteins = p_proteins.duplicate(true)
	assert_nonnegative()

func deep_copy():
	return ExpressionState.new(mrna, proteins)

func mrna_for(locus_id: int) -> float:
	return float(mrna.get(locus_id, 0.0))

func protein_for(locus_id: int) -> float:
	return float(proteins.get(locus_id, 0.0))

func total_mrna() -> float:
	var result: float = 0.0
	for value in mrna.values():
		result += float(value)
	return result

func total_protein() -> float:
	var result: float = 0.0
	for value in proteins.values():
		result += float(value)
	return result

func assert_matches_genome(genome) -> void:
	assert(genome != null)
	for gene in genome.genes:
		assert(mrna.has(int(gene.locus_id)), "Missing mRNA state for locus %d" % int(gene.locus_id))
		assert(proteins.has(int(gene.locus_id)), "Missing protein state for locus %d" % int(gene.locus_id))

func assert_nonnegative() -> void:
	for value in mrna.values():
		assert(float(value) >= -1e-10, "Negative mRNA abundance")
	for value in proteins.values():
		assert(float(value) >= -1e-10, "Negative protein abundance")

func checksum() -> float:
	var loci: Array = mrna.keys()
	loci.sort()
	var result: float = 0.0
	for i in range(loci.size()):
		var locus_id: int = int(loci[i])
		result += mrna_for(locus_id) * float((i + 3) * 29)
		result += protein_for(locus_id) * float((i + 5) * 43)
	return result
