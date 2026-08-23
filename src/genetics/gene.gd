extends RefCounted
class_name Gene

# Genotype is discrete so inheritance, mutation, hashing and replay do not depend
# on float serialization. M5-B adds a promoter-binding motif. Protein regulatory
# polarity/strength, sensor mode and stability are deliberately derived from the
# protein signature itself, so existing protein cohorts preserve those molecular
# properties after a DNA mutation without needing hidden gene-type flags.

const PROMOTER_CODE_MIN: int = 0
const PROMOTER_CODE_MAX: int = 10000
const SIGNATURE_MASK: int = 0xFFFF
const NEUTRAL_MARKER_MASK: int = 0x7FFFFFFF

var locus_id: int
var promoter_code: int
var protein_signature: int
var neutral_marker: int
var promoter_binding_motif: int

func _init(
	p_locus_id: int = 0,
	p_promoter_code: int = 5000,
	p_protein_signature: int = 0,
	p_neutral_marker: int = 0,
	p_promoter_binding_motif: int = 0xFFFF
) -> void:
	locus_id = p_locus_id
	promoter_code = clampi(p_promoter_code, PROMOTER_CODE_MIN, PROMOTER_CODE_MAX)
	protein_signature = p_protein_signature & SIGNATURE_MASK
	neutral_marker = p_neutral_marker & NEUTRAL_MARKER_MASK
	promoter_binding_motif = p_promoter_binding_motif & SIGNATURE_MASK
	validate()

func validate() -> void:
	assert(locus_id > 0, "Gene locus_id must be positive")
	assert(promoter_code >= PROMOTER_CODE_MIN and promoter_code <= PROMOTER_CODE_MAX)
	assert(protein_signature >= 0 and protein_signature <= SIGNATURE_MASK)
	assert(neutral_marker >= 0 and neutral_marker <= NEUTRAL_MARKER_MASK)
	assert(promoter_binding_motif >= 0 and promoter_binding_motif <= SIGNATURE_MASK)

func promoter_strength() -> float:
	return float(promoter_code) / float(PROMOTER_CODE_MAX)

func deep_copy():
	return Gene.new(locus_id, promoter_code, protein_signature, neutral_marker, promoter_binding_motif)

func canonical_key() -> String:
	return "%d:%d:%d:%d:%d" % [locus_id, promoter_code, protein_signature, neutral_marker, promoter_binding_motif]

func checksum() -> float:
	return (
		float(locus_id) * 3.0
		+ float(promoter_code) * 0.007
		+ float(protein_signature) * 0.00011
		+ float(neutral_marker % 1000003) * 0.000001
		+ float(promoter_binding_motif) * 0.000037
	)
