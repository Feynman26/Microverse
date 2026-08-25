extends RefCounted
class_name ExtracellularProteinTurnover

# Generic first-order environmental protein hydrolysis. The rate is supplied by
# the experiment/configuration layer; this module assigns no biological meaning
# to a sequence and gives no protein special stability. Lost protein material is
# returned locally to the extracellular amino-acid pool using the same modeled
# AA material per protein unit as intracellular translation/degradation.
static func step(
	world,
	dt: float,
	decay_rate_per_min: float,
	aa_material_per_protein: float,
	aa_field_name: String = "amino_acids"
) -> Dictionary:
	assert(dt >= 0.0)
	assert(decay_rate_per_min >= 0.0)
	assert(aa_material_per_protein >= 0.0)
	assert(world.has_field(aa_field_name), "Extracellular protein turnover requires an AA field")
	if dt <= 0.0 or decay_rate_per_min <= 0.0 or world.protein_fields.is_empty():
		return {"degraded_protein": 0.0, "recycled_aa": 0.0, "removed_signatures": []}

	var decay_fraction: float = 1.0 - exp(-decay_rate_per_min * dt)
	var degraded_total: float = 0.0
	var recycled_total: float = 0.0
	var removed_signatures: Array = []
	var signatures: Array = world.protein_signatures()
	var aa_field = world.get_field(aa_field_name)

	for signature_variant in signatures:
		var signature: int = int(signature_variant)
		var field = world.get_protein_field(signature)
		var degraded_signature: float = 0.0
		for y in range(world.height):
			for x in range(world.width):
				var current: float = maxf(0.0, float(field.get_value(x, y)))
				if current <= 0.0:
					continue
				var degraded: float = minf(current, current * decay_fraction)
				if degraded <= 0.0:
					continue
				field.set_value(x, y, current - degraded)
				var recycled: float = degraded * aa_material_per_protein
				if recycled > 0.0:
					aa_field.add_amount(x, y, recycled)
				degraded_signature += degraded
				recycled_total += recycled
		degraded_total += degraded_signature

		# First-order decay is asymptotic. Once only floating-point dust remains,
		# recycle that residue before removing the dynamic sequence field so long
		# evolutionary runs do not accumulate dead protein signatures forever.
		var residual: float = float(field.total_amount())
		if residual <= 1e-12:
			if residual > 0.0:
				for y in range(world.height):
					for x in range(world.width):
						var amount: float = maxf(0.0, float(field.get_value(x, y)))
						if amount > 0.0:
							aa_field.add_amount(x, y, amount * aa_material_per_protein)
							degraded_total += amount
							recycled_total += amount * aa_material_per_protein
			world.protein_fields.erase(signature)
			removed_signatures.append(signature)

	return {
		"degraded_protein": degraded_total,
		"recycled_aa": recycled_total,
		"removed_signatures": removed_signatures,
		"decay_fraction": decay_fraction
	}
