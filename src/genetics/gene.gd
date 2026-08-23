extends RefCounted
class_name Gene

const PROMOTER_CODE_MIN: int = 0
const PROMOTER_CODE_MAX: int = 10000
const SIGNATURE_MASK: int = 0xFFFF
const NEUTRAL_MARKER_MASK: int = 0x7FFFFFFF
const REGION_COPY_MAX: int = 8

var locus_id: int
var promoter_code: int
var protein_signature: int
var neutral_marker: int
# M5-B: inherited molecular motif recognized by regulatory protein cohorts.
var regulatory_signature: int
# M10 cis-regulatory architecture. A value of one reproduces every historical
# genome exactly. Duplication/deletion changes copy number without creating a
# special phenotype flag; promoter copies change transcription opportunity and
# regulatory copies change the number of identical binding sites.
var promoter_copy_number: int
var regulatory_copy_number: int

func _init(
	p_locus_id: int = 0,
	p_promoter_code: int = 5000,
	p_protein_signature: int = 0,
	p_neutral_marker: int = 0,
	p_regulatory_signature: int = 0,
	p_promoter_copy_number: int = 1,
	p_regulatory_copy_number: int = 1
) -> void:
	locus_id = p_locus_id
	promoter_code = clampi(p_promoter_code, PROMOTER_CODE_MIN, PROMOTER_CODE_MAX)
	protein_signature = p_protein_signature & SIGNATURE_MASK
	neutral_marker = p_neutral_marker & NEUTRAL_MARKER_MASK
	regulatory_signature = p_regulatory_signature & SIGNATURE_MASK
	promoter_copy_number = p_promoter_copy_number
	regulatory_copy_number = p_regulatory_copy_number
	validate()

func validate() -> void:
	assert(locus_id > 0, "Gene locus_id must be positive")
	assert(promoter_code >= PROMOTER_CODE_MIN and promoter_code <= PROMOTER_CODE_MAX)
	assert(protein_signature >= 0 and protein_signature <= SIGNATURE_MASK)
	assert(neutral_marker >= 0 and neutral_marker <= NEUTRAL_MARKER_MASK)
	assert(regulatory_signature >= 0 and regulatory_signature <= SIGNATURE_MASK)
	assert(promoter_copy_number >= 0 and promoter_copy_number <= REGION_COPY_MAX)
	assert(regulatory_copy_number >= 0 and regulatory_copy_number <= REGION_COPY_MAX)

func promoter_strength() -> float:
	# Multiple identical promoter copies create multiple transcription-initiation
	# opportunities. Zero copies is a true cis-promoter deletion.
	return (float(promoter_code) / float(PROMOTER_CODE_MAX)) * float(promoter_copy_number)

func deep_copy():
	return Gene.new(
		locus_id, promoter_code, protein_signature, neutral_marker, regulatory_signature,
		promoter_copy_number, regulatory_copy_number
	)

func canonical_key() -> String:
	var base: String = "%d:%d:%d:%d:%d" % [locus_id, promoter_code, protein_signature, neutral_marker, regulatory_signature]
	# Preserve the exact historical canonical representation for the overwhelmingly
	# common single-copy architecture. Only an actual M10 cis-architecture change
	# extends the key, so existing M3-M9 freezer fixtures remain stable.
	if promoter_copy_number == 1 and regulatory_copy_number == 1:
		return base
	return "%s:PC%d:RC%d" % [base, promoter_copy_number, regulatory_copy_number]

func checksum() -> float:
	var result: float = (
		float(locus_id) * 3.0
		+ float(promoter_code) * 0.007
		+ float(protein_signature) * 0.00011
		+ float(neutral_marker % 1000003) * 0.000001
		+ float(regulatory_signature) * 0.000013
	)
	# Single-copy historical genes retain their exact old checksum. Region-copy
	# changes become authoritative only when architecture actually diverges.
	result += float(promoter_copy_number - 1) * 0.000017
	result += float(regulatory_copy_number - 1) * 0.000019
	return result
