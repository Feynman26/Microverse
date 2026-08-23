extends RefCounted
class_name GeneticAssay

const GenomeScript = preload("res://src/genetics/genome.gd")

static func clone_genome(genome) -> Dictionary:
	assert(genome != null)
	return {
		"genome": genome.deep_copy(),
		"intervention": {
			"kind": "genetic_intervention",
			"operation": "clone_genotype",
			"source_genotype": genome.canonical_key()
		}
	}

static func revert_mutation(genome, mutation_event: Dictionary) -> Dictionary:
	assert(genome != null)
	assert(mutation_event.has("mutation_type") and mutation_event.has("locus_id"))
	assert(mutation_event.has("old_value") and mutation_event.has("new_value"))
	var result = genome.deep_copy()
	var locus_id: int = int(mutation_event["locus_id"])
	var gene = result.get_gene_by_locus(locus_id)
	assert(gene != null, "Cannot revert mutation at missing locus")
	var mutation_type: String = String(mutation_event["mutation_type"])
	var current_value: int = _mutation_value(gene, mutation_type)
	assert(current_value == int(mutation_event["new_value"]), "Reversion target no longer carries the mutation's derived allele")
	_set_mutation_value(gene, mutation_type, int(mutation_event["old_value"]))
	result.validate()
	return {
		"genome": result,
		"intervention": _intervention("revert_mutation", locus_id, mutation_type, current_value, int(mutation_event["old_value"]))
	}

static func introduce_mutation(genome, mutation_event: Dictionary) -> Dictionary:
	assert(genome != null)
	assert(mutation_event.has("mutation_type") and mutation_event.has("locus_id") and mutation_event.has("new_value"))
	var result = genome.deep_copy()
	var locus_id: int = int(mutation_event["locus_id"])
	var mutation_type: String = String(mutation_event["mutation_type"])
	var gene = result.get_gene_by_locus(locus_id)
	assert(gene != null, "Cannot introduce mutation at missing locus")
	var old_value: int = _mutation_value(gene, mutation_type)
	if mutation_event.has("old_value"):
		assert(old_value == int(mutation_event["old_value"]), "Mutation background does not carry the declared ancestral allele")
	_set_mutation_value(gene, mutation_type, int(mutation_event["new_value"]))
	result.validate()
	return {
		"genome": result,
		"intervention": _intervention("introduce_mutation", locus_id, mutation_type, old_value, int(mutation_event["new_value"]))
	}

static func knock_out_gene(genome, locus_id: int) -> Dictionary:
	assert(genome != null and locus_id > 0)
	var result = genome.deep_copy()
	var gene = result.get_gene_by_locus(locus_id)
	assert(gene != null, "Cannot knock out missing locus")
	var old_value: int = int(gene.promoter_code)
	gene.promoter_code = 0
	result.validate()
	return {
		"genome": result,
		"intervention": _intervention("knockout_promoter", locus_id, "promoter_code", old_value, 0)
	}

static func restore_ancestral_allele(genome, locus_id: int) -> Dictionary:
	assert(genome != null and locus_id > 0)
	var result = genome.deep_copy()
	var ancestor = GenomeScript.create_ancestor()
	var target = result.get_gene_by_locus(locus_id)
	var ancestral = ancestor.get_gene_by_locus(locus_id)
	assert(target != null and ancestral != null, "Ancestral restoration requires a canonical ancestral locus")
	var before: String = target.canonical_key()
	target.promoter_code = int(ancestral.promoter_code)
	target.protein_signature = int(ancestral.protein_signature)
	target.neutral_marker = int(ancestral.neutral_marker)
	target.regulatory_signature = int(ancestral.regulatory_signature)
	result.validate()
	return {
		"genome": result,
		"intervention": {
			"kind": "genetic_intervention",
			"operation": "restore_ancestral_allele",
			"locus_id": locus_id,
			"before": before,
			"after": target.canonical_key()
		}
	}

static func differences(reference, candidate) -> Array:
	assert(reference != null and candidate != null)
	var differences_found: Array = []
	var locus_ids: Dictionary = {}
	for gene in reference.genes:
		locus_ids[int(gene.locus_id)] = true
	for gene in candidate.genes:
		locus_ids[int(gene.locus_id)] = true
	var loci: Array = locus_ids.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var a = reference.get_gene_by_locus(locus_id)
		var b = candidate.get_gene_by_locus(locus_id)
		if a == null or b == null:
			differences_found.append({"locus_id": locus_id, "field": "presence", "reference": a != null, "candidate": b != null})
			continue
		for field_name in ["promoter_code", "protein_signature", "regulatory_signature", "neutral_marker"]:
			var av: int = int(a.get(field_name))
			var bv: int = int(b.get(field_name))
			if av != bv:
				differences_found.append({"locus_id": locus_id, "field": field_name, "reference": av, "candidate": bv})
	return differences_found

static func single_reversion_candidates(reference, candidate) -> Array:
	var result: Array = []
	for difference_variant in differences(reference, candidate):
		var difference: Dictionary = difference_variant
		if String(difference["field"]) == "presence":
			continue
		var mutation_type: String = _field_to_mutation_type(String(difference["field"]))
		var event: Dictionary = {
			"mutation_type": mutation_type,
			"locus_id": int(difference["locus_id"]),
			"old_value": int(difference["reference"]),
			"new_value": int(difference["candidate"])
		}
		var reverted: Dictionary = revert_mutation(candidate, event)
		result.append({"difference": difference, "genome": reverted["genome"], "intervention": reverted["intervention"]})
	return result

static func _mutation_value(gene, mutation_type: String) -> int:
	match mutation_type:
		"promoter_code": return int(gene.promoter_code)
		"protein_signature_bit_flip": return int(gene.protein_signature)
		"regulatory_signature_bit_flip": return int(gene.regulatory_signature)
		"neutral_marker": return int(gene.neutral_marker)
		_:
			assert(false, "Unsupported experimental mutation type: %s" % mutation_type)
	return 0

static func _set_mutation_value(gene, mutation_type: String, value: int) -> void:
	match mutation_type:
		"promoter_code": gene.promoter_code = value
		"protein_signature_bit_flip": gene.protein_signature = value & 0xFFFF
		"regulatory_signature_bit_flip": gene.regulatory_signature = value & 0xFFFF
		"neutral_marker": gene.neutral_marker = value & 0x7FFFFFFF
		_:
			assert(false, "Unsupported experimental mutation type: %s" % mutation_type)

static func _field_to_mutation_type(field_name: String) -> String:
	match field_name:
		"promoter_code": return "promoter_code"
		"protein_signature": return "protein_signature_bit_flip"
		"regulatory_signature": return "regulatory_signature_bit_flip"
		"neutral_marker": return "neutral_marker"
		_:
			assert(false, "No mutation event mapping for field: %s" % field_name)
	return ""

static func _intervention(operation: String, locus_id: int, mutation_type: String, old_value: int, new_value: int) -> Dictionary:
	return {
		"kind": "genetic_intervention",
		"operation": operation,
		"locus_id": locus_id,
		"mutation_type": mutation_type,
		"old_value": old_value,
		"new_value": new_value
	}
