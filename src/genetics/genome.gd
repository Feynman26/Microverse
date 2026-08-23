extends RefCounted
class_name Genome

const GeneScript = preload("res://src/genetics/gene.gd")
const FINGERPRINT_MODULUS: int = 2147483647

var genes: Array = []

func _init(p_genes: Array = []) -> void:
	for gene in p_genes:
		genes.append(gene.deep_copy())
	validate()

static func create_ancestor():
	# Each ancestral promoter motif is outside the active regulatory-binding
	# radius of the complete ancestral proteome. Its own locus protein is exactly
	# distance 5 while regulation activates at <=4. Thus the ancestor begins with
	# no hidden regulation, yet a single appropriate coding or motif bit mutation
	# can create a new regulatory edge.
	var ancestor_genes: Array = [
		GeneScript.new(1, 6200, 0x1357, 101, 0xD359),
		GeneScript.new(2, 5400, 0x2468, 102, 0xF768),
		GeneScript.new(3, 7100, 0x369C, 103, 0xF7D8),
		GeneScript.new(4, 4800, 0x48AD, 104, 0xC5A9),
		GeneScript.new(5, 3900, 0x5ACE, 105, 0xD3C8),
		GeneScript.new(6, 6600, 0x6BDF, 106, 0xE7CB),
		GeneScript.new(7, 5700, 0x7CE1, 107, 0xF7E9),
		GeneScript.new(8, 4500, 0x8DF2, 108, 0xF5FA),
		GeneScript.new(9, 5200, 0x9E03, 109, 0xF70B),
		GeneScript.new(10, 6000, 0xAF14, 110, 0xE719),
		GeneScript.new(11, 4300, 0xB025, 111, 0xF329),
		GeneScript.new(12, 5000, 0xC136, 112, 0xF33A)
	]
	return Genome.new(ancestor_genes)

func validate() -> void:
	var seen: Dictionary = {}
	for gene in genes:
		gene.validate()
		assert(not seen.has(gene.locus_id), "Duplicate locus_id in genome")
		seen[gene.locus_id] = true

func deep_copy():
	return Genome.new(genes)

func gene_count() -> int:
	return genes.size()

func get_gene_by_locus(locus_id: int):
	for gene in genes:
		if int(gene.locus_id) == locus_id:
			return gene
	return null

func canonical_key() -> String:
	var ordered: Array = genes.duplicate()
	ordered.sort_custom(func(a, b): return int(a.locus_id) < int(b.locus_id))
	var parts: PackedStringArray = PackedStringArray()
	for gene in ordered:
		parts.append(gene.canonical_key())
	return "|".join(parts)

# Deterministic rolling fingerprint. Canonical key remains collision-resolving.
func fingerprint() -> int:
	var ordered: Array = genes.duplicate()
	ordered.sort_custom(func(a, b): return int(a.locus_id) < int(b.locus_id))
	var value: int = 104729
	for gene in ordered:
		value = int((value * 131 + int(gene.locus_id) * 17 + int(gene.promoter_code)) % FINGERPRINT_MODULUS)
		value = int((value * 131 + int(gene.protein_signature)) % FINGERPRINT_MODULUS)
		value = int((value * 131 + int(gene.promoter_binding_motif)) % FINGERPRINT_MODULUS)
		value = int((value * 131 + int(gene.neutral_marker % 1000003)) % FINGERPRINT_MODULUS)
	return value

func exact_equals(other) -> bool:
	if other == null or genes.size() != other.genes.size():
		return false
	return canonical_key() == other.canonical_key()

func checksum() -> float:
	var ordered: Array = genes.duplicate()
	ordered.sort_custom(func(a, b): return int(a.locus_id) < int(b.locus_id))
	var result: float = float(ordered.size()) * 0.17
	for i in range(ordered.size()):
		result += float(ordered[i].checksum()) * float(i + 1)
	return result
