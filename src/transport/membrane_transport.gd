extends RefCounted
class_name MembraneTransport

const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

# Transport recognition uses the same compact 16-bit protein sequence space as
# catalysis/regulation, but in a distinct molecular-recognition domain. The XOR
# transform was characterized against the M6 ancestor: every M7 secondary
# transport target is exactly Hamming distance 5 from its nearest ancestral
# protein. With ACTIVE_MAX_DISTANCE=4 the ancestor therefore receives no
# secondary transport for free, while an ordinary one-bit coding mutation can
# make each capability reachable.
const TRANSPORT_DOMAIN_MASK: int = 0x4190
const ACTIVE_MAX_DISTANCE: int = 4
const DISTANCE_DECAY: float = 0.70

# Pure molecular recognition cache. A key contains only protein and target
# signatures; no concentration, resource, RNG or phenotype state is cached.
static var _affinity_cache: Dictionary = {}

static func hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & 0xFFFF
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

static func target_signature(metabolite_id: String) -> int:
	assert(MetaboliteCatalogScript.has_extracellular_field(metabolite_id))
	return (MetaboliteCatalogScript.ligand_signature(metabolite_id) ^ TRANSPORT_DOMAIN_MASK) & 0xFFFF

static func affinity(protein_signature: int, metabolite_id: String) -> float:
	var protein: int = protein_signature & 0xFFFF
	var target: int = target_signature(metabolite_id)
	return affinity_to_target(protein, target)

static func affinity_to_target(protein_signature: int, target: int) -> float:
	var protein: int = protein_signature & 0xFFFF
	var key: int = (protein << 16) | target
	if _affinity_cache.has(key):
		return float(_affinity_cache[key])
	var distance: int = hamming_distance(protein, target)
	var result: float = 0.0 if distance > ACTIVE_MAX_DISTANCE else exp(-DISTANCE_DECAY * float(distance))
	_affinity_cache[key] = result
	return result

# Transport is a property of proteins that physically exist, not of DNA alone.
# A coding mutation therefore changes transport only after expression builds a
# compatible protein cohort, preserving the M5 genotype->proteome causal chain.
static func proteome_activity(expression_state: Dictionary, metabolite_id: String, config) -> float:
	var result: float = 0.0
	var loci: Array = expression_state.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var protein_signature: int = int(signature_variant)
			var abundance: float = maxf(0.0, float(cohorts[signature_variant])) / float(config.expression_reference_protein_count)
			result += abundance * affinity(protein_signature, metabolite_id)
	return result

# P3-B exact dense recognition. Every metabolite accumulator sees protein
# cohorts in precisely the same sorted locus/signature order as the historical
# scalar function above. The accumulations are merely interleaved, allowing one
# proteome traversal and one pair of sorts per cell instead of one per target.
static func proteome_activities(
	expression_state: Dictionary,
	metabolite_ids: Array[String],
	config
) -> PackedFloat64Array:
	return proteome_activities_for_targets(
		expression_state,
		transport_target_signatures(metabolite_ids),
		config,
		{}
	)

static func transport_target_signatures(metabolite_ids: Array[String]) -> PackedInt32Array:
	var targets := PackedInt32Array()
	targets.resize(metabolite_ids.size())
	for target_index in range(metabolite_ids.size()):
		targets[target_index] = target_signature(metabolite_ids[target_index])
	return targets

# affinity_rows belongs to one compiled target vector in one simulation engine.
# It is a derived cache keyed only by physical 16-bit protein signature.
static func proteome_activities_for_targets(
	expression_state: Dictionary,
	targets: PackedInt32Array,
	config,
	affinity_rows: Dictionary
) -> PackedFloat64Array:
	var target_count: int = targets.size()
	var result := PackedFloat64Array()
	result.resize(target_count)

	var loci: Array = expression_state.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var protein_signature: int = int(signature_variant)
			var abundance: float = maxf(0.0, float(cohorts[signature_variant])) / float(config.expression_reference_protein_count)
			var affinity_row: PackedFloat64Array
			if affinity_rows.has(protein_signature):
				affinity_row = affinity_rows[protein_signature]
			else:
				affinity_row = PackedFloat64Array()
				affinity_row.resize(target_count)
				for target_index in range(target_count):
					affinity_row[target_index] = affinity_to_target(protein_signature, targets[target_index])
				affinity_rows[protein_signature] = affinity_row
			for target_index in range(target_count):
				result[target_index] += abundance * affinity_row[target_index]
	return result

# Signed desired exchange: positive means extracellular -> intracellular,
# negative means intracellular -> extracellular. Direction is determined only
# by the concentration gradient; there is no import/export behavioral flag.
static func desired_exchange(
	internal_amount: float,
	volume: float,
	external_amount: float,
	activity: float,
	dt: float,
	config
) -> float:
	assert(internal_amount >= 0.0)
	assert(volume > 0.0)
	assert(external_amount >= 0.0)
	assert(activity >= 0.0)
	assert(dt >= 0.0)
	if activity <= 0.0 or dt <= 0.0:
		return 0.0

	var internal_concentration: float = internal_amount / volume
	var gradient: float = external_amount - internal_concentration
	var magnitude: float = absf(gradient)
	if magnitude <= 1e-12:
		return 0.0
	var saturation: float = magnitude / (float(config.secondary_transport_gradient_km) + magnitude)
	var capacity: float = activity * float(config.secondary_transport_vmax_per_reference_protein) * volume * dt
	var proposed: float = capacity * saturation
	if gradient > 0.0:
		var intracellular_capacity: float = maxf(
			0.0,
			float(config.intracellular_pool_capacity_per_volume) * volume - internal_amount
		)
		return minf(proposed, intracellular_capacity)
	return -minf(proposed, internal_amount)

static func total_movement(exchanges: Dictionary) -> float:
	var result: float = 0.0
	for value in exchanges.values():
		result += absf(float(value))
	return result

static func movement_cost(moved_units: float, config) -> float:
	return maxf(0.0, moved_units) * float(config.secondary_transport_atp_cost_per_unit)

# All simultaneous secondary transport proposals from one cell share its ATP
# pool proportionally. No metabolite is privileged by iteration order.
static func energy_scale(exchanges: Dictionary, available_atp: float, config) -> float:
	var demand: float = movement_cost(total_movement(exchanges), config)
	if demand <= 0.0:
		return 1.0
	return minf(1.0, maxf(0.0, available_atp) / demand)

# Array equivalent used by the dense allocator. Iteration order is the
# canonical metabolite order, matching Dictionary.values() in the legacy path.
static func energy_scale_dense(exchanges: PackedFloat64Array, available_atp: float, config) -> float:
	var moved: float = 0.0
	for signed_amount in exchanges:
		moved += absf(signed_amount)
	var demand: float = movement_cost(moved, config)
	if demand <= 0.0:
		return 1.0
	return minf(1.0, maxf(0.0, available_atp) / demand)
