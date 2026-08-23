extends RefCounted
class_name MetabolicSolver

const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

static func create_initial_pools(volume: float, config) -> Dictionary:
	var pools: Dictionary = {}
	for metabolite_id in MetaboliteCatalogScript.ids(): pools[metabolite_id] = 0.0
	pools["BIO"] = volume * float(config.biomass_units_per_volume)
	pools["ATP"] = float(config.initial_atp_per_volume) * volume
	pools["ADP"] = float(config.initial_adp_per_volume) * volume
	pools["NAD"] = float(config.initial_nad_per_volume) * volume
	pools["NADH"] = float(config.initial_nadh_per_volume) * volume
	return pools

# Compatibility hierarchy:
# - protein_cohorts: authoritative M5 molecular proteome preserving old protein identity;
# - protein_abundance: M4/M5 controlled assays using current locus signature;
# - neither: legacy M4 promoter-as-abundance characterization only.
static func step(
	pools: Dictionary,
	genome,
	reactions: Array,
	dt: float,
	volume: float,
	config,
	protein_abundance: Dictionary = {},
	protein_cohorts: Dictionary = {}
) -> Dictionary:
	assert(genome != null)
	assert(dt >= 0.0)
	assert(volume > 0.0)
	var cumulative_fluxes: Dictionary = {}
	for reaction in reactions: cumulative_fluxes[reaction.reaction_id] = 0.0
	if dt <= 0.0: return cumulative_fluxes
	var substeps: int = maxi(1, int(config.metabolic_substeps_per_tick))
	var sub_dt: float = dt / float(substeps)
	for _substep in range(substeps):
		var fluxes: Dictionary = _solve_substep(pools, genome, reactions, sub_dt, volume, config, protein_abundance, protein_cohorts)
		for reaction_id in fluxes.keys(): cumulative_fluxes[reaction_id] = float(cumulative_fluxes[reaction_id]) + float(fluxes[reaction_id])
		assert_nonnegative(pools)
	return cumulative_fluxes

static func _solve_substep(
	pools: Dictionary,
	genome,
	reactions: Array,
	dt: float,
	volume: float,
	config,
	protein_abundance: Dictionary,
	protein_cohorts: Dictionary
) -> Dictionary:
	var snapshot: Dictionary = pools.duplicate(true)
	var potential_flux: Dictionary = {}
	var substrate_demand: Dictionary = {}
	for reaction in reactions:
		var activity: float
		if not protein_cohorts.is_empty():
			activity = CatalyticLandscapeScript.cohort_activity(protein_cohorts, reaction)
		elif not protein_abundance.is_empty():
			activity = CatalyticLandscapeScript.proteome_activity(genome, protein_abundance, reaction)
		else:
			activity = CatalyticLandscapeScript.genome_activity(genome, reaction)
		var saturation: float = _limiting_saturation(snapshot, reaction.substrates, volume, float(config.metabolic_km_per_volume))
		var capacity: float = activity * float(config.metabolic_rate_scale) * volume * dt
		var requested_flux: float = maxf(0.0, capacity * saturation)
		potential_flux[reaction.reaction_id] = requested_flux
		for metabolite_id in reaction.substrates.keys():
			var demand: float = requested_flux * float(reaction.substrates[metabolite_id])
			substrate_demand[metabolite_id] = float(substrate_demand.get(metabolite_id, 0.0)) + demand

	var substrate_scale: Dictionary = {}
	for metabolite_id in substrate_demand.keys():
		var available: float = maxf(0.0, float(snapshot.get(metabolite_id, 0.0)))
		var demand: float = float(substrate_demand[metabolite_id])
		substrate_scale[metabolite_id] = 1.0 if demand <= available or demand <= 0.0 else available / demand

	var effective_flux: Dictionary = {}
	var deltas: Dictionary = {}
	for reaction in reactions:
		var flux: float = float(potential_flux[reaction.reaction_id])
		for metabolite_id in reaction.substrates.keys():
			flux = minf(flux, float(potential_flux[reaction.reaction_id]) * float(substrate_scale.get(metabolite_id, 1.0)))
		effective_flux[reaction.reaction_id] = flux
		if flux <= 0.0: continue
		for metabolite_id in reaction.substrates.keys():
			deltas[metabolite_id] = float(deltas.get(metabolite_id, 0.0)) - flux * float(reaction.substrates[metabolite_id])
		for metabolite_id in reaction.products.keys():
			deltas[metabolite_id] = float(deltas.get(metabolite_id, 0.0)) + flux * float(reaction.products[metabolite_id])
	for metabolite_id in deltas.keys():
		pools[metabolite_id] = float(pools.get(metabolite_id, 0.0)) + float(deltas[metabolite_id])
		if float(pools[metabolite_id]) < 0.0 and float(pools[metabolite_id]) > -1e-10: pools[metabolite_id] = 0.0
	return effective_flux

static func _limiting_saturation(snapshot: Dictionary, substrates: Dictionary, volume: float, km_per_volume: float) -> float:
	var result: float = 1.0
	var km: float = maxf(1e-12, km_per_volume * volume)
	for metabolite_id in substrates.keys():
		var coefficient: float = float(substrates[metabolite_id])
		var available_equivalents: float = maxf(0.0, float(snapshot.get(metabolite_id, 0.0))) / coefficient
		var saturation: float = available_equivalents / (km + available_equivalents)
		result = minf(result, saturation)
	return result

static func spend_atp(pools: Dictionary, requested_amount: float) -> float:
	var requested: float = maxf(0.0, requested_amount)
	var available: float = maxf(0.0, float(pools.get("ATP", 0.0)))
	var spent: float = minf(requested, available)
	pools["ATP"] = available - spent
	pools["ADP"] = float(pools.get("ADP", 0.0)) + spent
	return spent

static func add_pool(pools: Dictionary, metabolite_id: String, amount: float) -> void:
	assert(MetaboliteCatalogScript.has(metabolite_id))
	assert(amount >= 0.0)
	pools[metabolite_id] = float(pools.get(metabolite_id, 0.0)) + amount

static func partition(source: Dictionary, ratio: float) -> Array:
	assert(ratio > 0.0 and ratio < 1.0)
	var first: Dictionary = {}
	var second: Dictionary = {}
	for metabolite_id in MetaboliteCatalogScript.ids():
		var amount: float = float(source.get(metabolite_id, 0.0))
		first[metabolite_id] = amount * ratio
		second[metabolite_id] = amount * (1.0 - ratio)
	return [first, second]

static func structural_totals(pools: Dictionary) -> Dictionary:
	var result: Dictionary = {"C": 0.0, "N": 0.0, "P": 0.0}
	for metabolite_id in MetaboliteCatalogScript.ids():
		var amount: float = float(pools.get(metabolite_id, 0.0))
		var units: Dictionary = MetaboliteCatalogScript.structural_units(metabolite_id)
		for element in result.keys(): result[element] = float(result[element]) + amount * float(units[element])
	return result

static func assert_nonnegative(pools: Dictionary) -> void:
	for metabolite_id in MetaboliteCatalogScript.ids():
		assert(float(pools.get(metabolite_id, 0.0)) >= -1e-10, "Negative metabolite pool %s=%s" % [metabolite_id, pools.get(metabolite_id, 0.0)])

static func checksum(pools: Dictionary) -> float:
	var ids: Array[String] = MetaboliteCatalogScript.ids()
	var result: float = 0.0
	for i in range(ids.size()): result += float(pools.get(ids[i], 0.0)) * float((i + 3) * 17)
	return result
