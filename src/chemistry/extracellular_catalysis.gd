extends RefCounted
class_name ExtracellularCatalysis

const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")

# The upper four sequence bits are a compressed secretion/localization motif.
# This is a molecular property, not a behavior flag. No ancestral protein has
# the Dxxx motif. An ordinary one-bit C136 -> D136 coding mutation can therefore
# create a secreted cohort while preserving all lower sequence information.
const SECRETION_SIGNAL_MASK: int = 0xF000
const SECRETION_SIGNAL_VALUE: int = 0xD000

static func has_secretion_signal(protein_signature: int) -> bool:
	return (protein_signature & SECRETION_SIGNAL_MASK) == SECRETION_SIGNAL_VALUE

# Requests are derived only from protein molecules that physically exist.
# DNA with a secretion-compatible coding signature does nothing until expression
# has produced the corresponding protein cohort.
static func secretion_proposals(expression_state: Dictionary, dt: float, config) -> Dictionary:
	var result: Dictionary = {}
	if dt <= 0.0:
		return result
	var loci: Array = expression_state.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			if not has_secretion_signal(signature):
				continue
			var abundance: float = maxf(0.0, float(cohorts[signature_variant]))
			var requested: float = minf(
				abundance,
				abundance * float(config.extracellular_protein_secretion_fraction_per_min) * dt
			)
			if requested > 0.0:
				result[signature] = float(result.get(signature, 0.0)) + requested
	return result

static func secretion_cost(proposals: Dictionary, config) -> float:
	var total: float = 0.0
	for amount in proposals.values():
		total += maxf(0.0, float(amount))
	return total * float(config.extracellular_protein_secretion_atp_cost_per_unit)

static func secretion_energy_scale(proposals: Dictionary, available_atp: float, config) -> float:
	var requested_cost: float = secretion_cost(proposals, config)
	if requested_cost <= 0.0:
		return 1.0
	return minf(1.0, maxf(0.0, available_atp) / requested_cost)

# Remove a signature proportionally from every locus that currently contains
# that cohort. This avoids first-locus privilege if duplicated genes later make
# the same physical protein sequence.
static func remove_protein_signature(expression_state: Dictionary, protein_signature: int, requested: float) -> float:
	var amount_requested: float = maxf(0.0, requested)
	if amount_requested <= 0.0:
		return 0.0
	var loci_with_signature: Array[int] = []
	var total_available: float = 0.0
	var loci: Array = expression_state.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		if cohorts.has(protein_signature):
			var abundance: float = maxf(0.0, float(cohorts[protein_signature]))
			if abundance > 0.0:
				loci_with_signature.append(locus_id)
				total_available += abundance
	if total_available <= 0.0:
		return 0.0

	var removed_total: float = minf(amount_requested, total_available)
	var fraction: float = removed_total / total_available
	for locus_id in loci_with_signature:
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		var old_amount: float = maxf(0.0, float(cohorts[protein_signature]))
		var new_amount: float = maxf(0.0, old_amount * (1.0 - fraction))
		if new_amount <= 1e-12:
			cohorts.erase(protein_signature)
		else:
			cohorts[protein_signature] = new_amount
	return removed_total

# Generic simultaneous extracellular catalytic solver. Secreted protein fields
# provide catalytic activity through the same Hamming-distance chemistry used
# intracellularly. Reactions read one local substrate snapshot, compete
# proportionally for shared substrates, and apply bounded nonnegative deltas.
static func step(world, reactions: Array, dt: float, config) -> Dictionary:
	var cumulative_fluxes: Dictionary = {}
	for reaction in reactions:
		cumulative_fluxes[reaction.reaction_id] = 0.0
	if dt <= 0.0 or reactions.is_empty():
		return {"fluxes": cumulative_fluxes, "total_flux": 0.0}

	var protein_signatures: Array = world.protein_signatures()
	if protein_signatures.is_empty():
		return {"fluxes": cumulative_fluxes, "total_flux": 0.0}

	# Precompute sequence contribution coefficients once per tick. A coefficient
	# of zero means that extracellular protein cannot catalyse that reaction.
	var coefficients: Dictionary = {}
	var any_active_pair: bool = false
	for reaction in reactions:
		var by_signature: Dictionary = {}
		for signature_variant in protein_signatures:
			var signature: int = int(signature_variant)
			var coefficient: float = (
				CatalyticLandscapeScript.affinity(signature, int(reaction.signature))
				* float(reaction.catalytic_ceiling)
				/ float(config.expression_reference_protein_count)
			)
			by_signature[signature] = coefficient
			if coefficient > 0.0:
				any_active_pair = true
		coefficients[reaction.reaction_id] = by_signature
	if not any_active_pair:
		return {"fluxes": cumulative_fluxes, "total_flux": 0.0}

	var total_flux: float = 0.0
	for y in range(world.height):
		for x in range(world.width):
			var potential_flux: Dictionary = {}
			var substrate_demand: Dictionary = {}

			for reaction in reactions:
				var activity: float = 0.0
				var by_signature: Dictionary = coefficients[reaction.reaction_id]
				for signature_variant in protein_signatures:
					var signature: int = int(signature_variant)
					var coefficient: float = float(by_signature[signature])
					if coefficient <= 0.0:
						continue
					var enzyme_amount: float = float(world.get_protein_field(signature).get_value(x, y))
					activity += enzyme_amount * coefficient
				if activity <= 0.0:
					potential_flux[reaction.reaction_id] = 0.0
					continue

				var saturation: float = _limiting_saturation(
					world,
					x,
					y,
					reaction.substrates,
					float(config.extracellular_catalysis_km)
				)
				var requested_flux: float = maxf(
					0.0,
					activity * float(config.extracellular_catalysis_rate_scale) * saturation * dt
				)
				potential_flux[reaction.reaction_id] = requested_flux
				for field_variant in reaction.substrates.keys():
					var field_name: String = String(field_variant)
					var demand: float = requested_flux * float(reaction.substrates[field_variant])
					substrate_demand[field_name] = float(substrate_demand.get(field_name, 0.0)) + demand

			if substrate_demand.is_empty():
				continue

			var substrate_scale: Dictionary = {}
			for field_variant in substrate_demand.keys():
				var field_name: String = String(field_variant)
				var available: float = maxf(0.0, float(world.get_field(field_name).get_value(x, y)))
				var demand: float = float(substrate_demand[field_variant])
				substrate_scale[field_name] = 1.0 if demand <= available or demand <= 0.0 else available / demand

			var deltas: Dictionary = {}
			for reaction in reactions:
				var requested_flux: float = float(potential_flux.get(reaction.reaction_id, 0.0))
				var effective_flux: float = requested_flux
				for field_variant in reaction.substrates.keys():
					var field_name: String = String(field_variant)
					effective_flux = minf(
						effective_flux,
						requested_flux * float(substrate_scale.get(field_name, 1.0))
					)
				if effective_flux <= 0.0:
					continue
				cumulative_fluxes[reaction.reaction_id] = float(cumulative_fluxes[reaction.reaction_id]) + effective_flux
				total_flux += effective_flux
				for field_variant in reaction.substrates.keys():
					var field_name: String = String(field_variant)
					deltas[field_name] = float(deltas.get(field_name, 0.0)) - effective_flux * float(reaction.substrates[field_variant])
				for field_variant in reaction.products.keys():
					var field_name: String = String(field_variant)
					deltas[field_name] = float(deltas.get(field_name, 0.0)) + effective_flux * float(reaction.products[field_variant])

			for field_variant in deltas.keys():
				var field_name: String = String(field_variant)
				var field = world.get_field(field_name)
				var current: float = float(field.get_value(x, y))
				var next_value: float = current + float(deltas[field_variant])
				assert(next_value >= -1e-10, "Extracellular catalysis produced negative field %s" % field_name)
				field.set_value(x, y, maxf(0.0, next_value))

	return {"fluxes": cumulative_fluxes, "total_flux": total_flux}

static func _limiting_saturation(world, x: int, y: int, substrates: Dictionary, km: float) -> float:
	var result: float = 1.0
	var effective_km: float = maxf(1e-12, km)
	for field_variant in substrates.keys():
		var field_name: String = String(field_variant)
		var coefficient: float = float(substrates[field_variant])
		var available: float = maxf(0.0, float(world.get_field(field_name).get_value(x, y)))
		var available_equivalents: float = available / coefficient
		var saturation: float = available_equivalents / (effective_km + available_equivalents)
		result = minf(result, saturation)
	return result
