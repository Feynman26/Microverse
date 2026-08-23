extends RefCounted
class_name ExpressionState

# mRNA is tracked by locus because a transcript belongs to a gene locus.
# Protein is tracked by locus *and the signature that was encoded when that
# cohort was translated*. This matters after mutation: inherited old proteins
# must not magically acquire the daughter's new DNA sequence.

var mrna: Dictionary = {}
var protein_cohorts: Dictionary = {} # locus_id -> {protein_signature: amount}

func _init(p_mrna: Dictionary = {}, p_protein_cohorts: Dictionary = {}) -> void:
	mrna = p_mrna.duplicate(true)
	protein_cohorts = p_protein_cohorts.duplicate(true)
	assert_nonnegative()

func deep_copy():
	return ExpressionState.new(mrna, protein_cohorts)

func mrna_for(locus_id: int) -> float:
	return float(mrna.get(locus_id, 0.0))

func protein_for(locus_id: int) -> float:
	var result: float = 0.0
	var cohorts: Dictionary = protein_cohorts.get(locus_id, {})
	for amount in cohorts.values():
		result += float(amount)
	return result

func protein_for_signature(locus_id: int, signature: int) -> float:
	var cohorts: Dictionary = protein_cohorts.get(locus_id, {})
	return float(cohorts.get(signature & 0xFFFF, 0.0))

func set_single_protein_cohort(locus_id: int, signature: int, amount: float) -> void:
	assert(amount >= 0.0)
	protein_cohorts[locus_id] = {signature & 0xFFFF: amount}

func add_protein(locus_id: int, signature: int, amount: float) -> void:
	assert(amount >= 0.0)
	var cohorts: Dictionary = protein_cohorts.get(locus_id, {}).duplicate(true)
	var normalized_signature: int = signature & 0xFFFF
	cohorts[normalized_signature] = float(cohorts.get(normalized_signature, 0.0)) + amount
	protein_cohorts[locus_id] = cohorts

func total_mrna() -> float:
	var result: float = 0.0
	for value in mrna.values(): result += float(value)
	return result

func total_protein() -> float:
	var result: float = 0.0
	for locus_variant in protein_cohorts.keys():
		result += protein_for(int(locus_variant))
	return result

func cohort_records() -> Array:
	var records: Array = []
	var loci: Array = protein_cohorts.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = protein_cohorts[locus_id]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			var amount: float = float(cohorts[signature])
			if amount > 0.0:
				records.append({"locus_id": locus_id, "signature": signature, "amount": amount})
	return records

func assert_matches_genome(genome) -> void:
	assert(genome != null)
	for gene in genome.genes:
		var locus_id: int = int(gene.locus_id)
		assert(mrna.has(locus_id), "Missing mRNA state for locus %d" % locus_id)
		assert(protein_cohorts.has(locus_id), "Missing protein cohort ledger for locus %d" % locus_id)

func assert_nonnegative() -> void:
	for value in mrna.values():
		assert(float(value) >= -1e-10, "Negative mRNA abundance")
	for cohorts_variant in protein_cohorts.values():
		var cohorts: Dictionary = cohorts_variant
		for amount in cohorts.values():
			assert(float(amount) >= -1e-10, "Negative protein cohort abundance")

func checksum() -> float:
	var loci: Array = mrna.keys()
	loci.sort()
	var result: float = 0.0
	for i in range(loci.size()):
		var locus_id: int = int(loci[i])
		result += mrna_for(locus_id) * float((i + 3) * 29)
		var cohorts: Dictionary = protein_cohorts.get(locus_id, {})
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for j in range(signatures.size()):
			var signature: int = int(signatures[j])
			result += float(cohorts[signature]) * float((i + 5) * 43 + (j + 1) * 7)
			result += float(signature) * 0.00000031 * float(i + 1)
	return result
