extends RefCounted
class_name Genome

const GeneScript = preload("res://src/genetics/gene.gd")
const FINGERPRINT_MODULUS: int = 2147483647
# Coarse DNA-copying units. The coding region is one unit; each promoter or
# regulatory site is one eighth of that unit. Absolute units are model units,
# but relative genome expansion/reduction now has an explicit copying burden.
const CIS_REGION_REPLICATION_UNITS: float = 0.125

var genes: Array = []

func _init(p_genes: Array = []) -> void:
	for gene in p_genes:
		genes.append(gene.deep_copy())
	validate()

static func create_ancestor():
	# M5-B regulatory motifs are deliberately dormant: each motif is distance 4
	# from its designated ancestral protein and at least distance 4 from every
	# ancestral protein. The active regulatory radius is 3, so the ancestor is
	# not gifted a regulatory circuit, while one ordinary motif mutation can
	# make an edge accessible.
	var ancestor_genes: Array = [
		GeneScript.new(1, 6200, 0x1357, 101, 0x1358),
		GeneScript.new(2, 5400, 0x2468, 102, 0x2467),
		GeneScript.new(3, 7100, 0x369C, 103, 0x3693),
		GeneScript.new(4, 4800, 0x48AD, 104, 0x48A2),
		GeneScript.new(5, 3900, 0x5ACE, 105, 0x5AC1),
		GeneScript.new(6, 6600, 0x6BDF, 106, 0x6BD0),
		GeneScript.new(7, 5700, 0x7CE1, 107, 0x7CEE),
		GeneScript.new(8, 4500, 0x8DF2, 108, 0x8DFD),
		GeneScript.new(9, 5200, 0x9E03, 109, 0x9E0C),
		GeneScript.new(10, 6000, 0xAF14, 110, 0xAF1B),
		GeneScript.new(11, 4300, 0xB025, 111, 0xB02A),
		GeneScript.new(12, 5000, 0xC136, 112, 0xC139)
	]
	return Genome.new(ancestor_genes)

func validate() -> void:
	assert(not genes.is_empty(), "Genome must retain at least one coding locus")
	var seen: Dictionary = {}
	for gene in genes:
		gene.validate()
		assert(not seen.has(gene.locus_id), "Duplicate locus_id in genome")
		seen[gene.locus_id] = true

func deep_copy():
	return Genome.new(genes)

func gene_count() -> int:
	return genes.size()

func replication_unit_count() -> float:
	var units: float = float(genes.size())
	for gene in genes:
		units += float(gene.promoter_copy_number) * CIS_REGION_REPLICATION_UNITS
		units += float(gene.regulatory_copy_number) * CIS_REGION_REPLICATION_UNITS
	return maxf(CIS_REGION_REPLICATION_UNITS, units)

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
		# Preserve historical fingerprints for single-copy cis architecture.
		if int(gene.promoter_copy_number) != 1 or int(gene.regulatory_copy_number) != 1:
			value = int((value * 131 + int(gene.promoter_copy_number) * 19 + int(gene.regulatory_copy_number) * 23) % FINGERPRINT_MODULUS)
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
