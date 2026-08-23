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
	# M5 assigns each promoter a regulatory motif. The motifs form a sparse ring
	# through existing protein signatures rather than named biological pathways:
	# each target can be influenced by a protein from another locus if that
	# protein is present. Because protein-signature high bit also determines
	# regulatory sign, the same sequence can carry catalytic and regulatory
	# pleiotropy and both properties remain mutable.
	var ancestor_genes: Array = [
		GeneScript.new(1, 6200, 0x1357, 101, 0xC136),
		GeneScript.new(2, 5400, 0x2468, 102, 0x1357),
		GeneScript.new(3, 7100, 0x369C, 103, 0x2468),
		GeneScript.new(4, 4800, 0x48AD, 104, 0x369C),
		GeneScript.new(5, 3900, 0x5ACE, 105, 0x48AD),
		GeneScript.new(6, 6600, 0x6BDF, 106, 0x5ACE),
		GeneScript.new(7, 5700, 0x7CE1, 107, 0x6BDF),
		GeneScript.new(8, 4500, 0x8DF2, 108, 0x7CE1),
		GeneScript.new(9, 5200, 0x9E03, 109, 0x8DF2),
		GeneScript.new(10, 6000, 0xAF14, 110, 0x9E03),
		GeneScript.new(11, 4300, 0xB025, 111, 0xAF14),
		GeneScript.new(12, 5000, 0xC136, 112, 0xB025)
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
	var parts: PackedStringArray = PackedStringArray()
	for gene in genes:
		parts.append(gene.canonical_key())
	return "|".join(parts)

func fingerprint() -> int:
	var value: int = 104729
	for gene in genes:
		value = int((value * 131 + int(gene.locus_id) * 17 + int(gene.promoter_code)) % FINGERPRINT_MODULUS)
		value = int((value * 131 + int(gene.protein_signature)) % FINGERPRINT_MODULUS)
		value = int((value * 131 + int(gene.regulatory_signature)) % FINGERPRINT_MODULUS)
		value = int((value * 131 + int(gene.neutral_marker % 1000003)) % FINGERPRINT_MODULUS)
	return value

func exact_equals(other) -> bool:
	if other == null or genes.size() != other.genes.size():
		return false
	return canonical_key() == other.canonical_key()

func checksum() -> float:
	var result: float = float(genes.size()) * 0.17
	for i in range(genes.size()):
		result += float(genes[i].checksum()) * float(i + 1)
	return result
