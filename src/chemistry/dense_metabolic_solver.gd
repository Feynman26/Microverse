extends RefCounted
class_name DenseMetabolicSolver

const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

# Exact P3 replacement for the dictionary-heavy M10 substep. Authoritative
# pools remain externally compatible dictionaries in this subgate, while all
# temporary solve state uses reusable dense Float64 buffers.

static func step(
	pools: Dictionary,
	genome,
	expression_state: Dictionary,
	reactions: Array,
	compiled_network,
	dt: float,
	volume: float,
	config,
	workspace: Dictionary
) -> Dictionary:
	assert(compiled_network != null)
	if dt <= 0.0:
		var zero_fluxes: Dictionary = {}
		for reaction_id in compiled_network.reaction_ids:
			zero_fluxes[reaction_id] = 0.0
		return zero_fluxes
	var metabolite_ids: Array[String] = MetaboliteCatalogScript.ORDERED_IDS
	var metabolite_count: int = metabolite_ids.size()
	var reaction_count: int = reactions.size()
	var pool_values: PackedFloat64Array = _take_buffer(workspace, "pool_values", metabolite_count)
	var snapshot_values: PackedFloat64Array = _take_buffer(workspace, "snapshot_values", metabolite_count)
	var substrate_demand: PackedFloat64Array = _take_buffer(workspace, "substrate_demand", metabolite_count)
	var substrate_scale: PackedFloat64Array = _take_buffer(workspace, "substrate_scale", metabolite_count)
	var deltas: PackedFloat64Array = _take_buffer(workspace, "deltas", metabolite_count)
	var catalytic_activities: PackedFloat64Array = _take_buffer(workspace, "catalytic_activities", reaction_count)
	var potential_flux: PackedFloat64Array = _take_buffer(workspace, "potential_flux", reaction_count)
	var cumulative_flux: PackedFloat64Array = _take_buffer(workspace, "cumulative_flux", reaction_count, true)

	for metabolite_index in range(metabolite_count):
		pool_values[metabolite_index] = float(pools.get(metabolite_ids[metabolite_index], 0.0))
	for reaction_index in range(reaction_count):
		catalytic_activities[reaction_index] = CatalyticLandscapeScript.proteome_activity(
			genome,
			expression_state,
			reactions[reaction_index],
			config
		)

	# Keep the numerical loop in the function that exclusively owns these
	# PackedFloat64Arrays. Passing them to a helper would introduce an extra
	# reference and trigger copy-on-write on the first mutation.
	var substeps: int = maxi(1, int(config.metabolic_substeps_per_tick))
	var sub_dt: float = dt / float(substeps)
	var km: float = maxf(1e-12, float(config.metabolic_km_per_volume) * volume)
	for _substep in range(substeps):
		for metabolite_index in range(metabolite_count):
			snapshot_values[metabolite_index] = pool_values[metabolite_index]
			substrate_demand[metabolite_index] = 0.0
			substrate_scale[metabolite_index] = 1.0
			deltas[metabolite_index] = 0.0

		for reaction_index in range(reaction_count):
			var saturation: float = 1.0
			var indices: PackedInt32Array = compiled_network.substrate_indices[reaction_index]
			var coefficients: PackedFloat64Array = compiled_network.substrate_coefficients[reaction_index]
			for item_index in range(indices.size()):
				var metabolite_index: int = indices[item_index]
				var available_equivalents: float = maxf(0.0, snapshot_values[metabolite_index]) / coefficients[item_index]
				var local_saturation: float = available_equivalents / (km + available_equivalents)
				saturation = minf(saturation, local_saturation)
			var capacity: float = catalytic_activities[reaction_index] * float(config.metabolic_rate_scale) * volume * sub_dt
			var requested_flux: float = maxf(0.0, capacity * saturation)
			potential_flux[reaction_index] = requested_flux
			for item_index in range(indices.size()):
				substrate_demand[indices[item_index]] += requested_flux * coefficients[item_index]

		for metabolite_index in range(metabolite_count):
			var demand: float = substrate_demand[metabolite_index]
			if demand > 0.0:
				var available: float = maxf(0.0, snapshot_values[metabolite_index])
				substrate_scale[metabolite_index] = 1.0 if demand <= available else available / demand

		for reaction_index in range(reaction_count):
			var requested_flux: float = potential_flux[reaction_index]
			var flux: float = requested_flux
			var substrate_indices: PackedInt32Array = compiled_network.substrate_indices[reaction_index]
			var substrate_coefficients: PackedFloat64Array = compiled_network.substrate_coefficients[reaction_index]
			for item_index in range(substrate_indices.size()):
				flux = minf(flux, requested_flux * substrate_scale[substrate_indices[item_index]])
			cumulative_flux[reaction_index] += flux
			if flux <= 0.0:
				continue
			for item_index in range(substrate_indices.size()):
				deltas[substrate_indices[item_index]] -= flux * substrate_coefficients[item_index]
			var product_indices: PackedInt32Array = compiled_network.product_indices[reaction_index]
			var product_coefficients: PackedFloat64Array = compiled_network.product_coefficients[reaction_index]
			for item_index in range(product_indices.size()):
				deltas[product_indices[item_index]] += flux * product_coefficients[item_index]

		for metabolite_index in range(metabolite_count):
			pool_values[metabolite_index] += deltas[metabolite_index]
			if pool_values[metabolite_index] < 0.0 and pool_values[metabolite_index] > -1e-10:
				pool_values[metabolite_index] = 0.0
			assert(pool_values[metabolite_index] >= -1e-10)

	for metabolite_index in range(metabolite_count):
		pools[metabolite_ids[metabolite_index]] = pool_values[metabolite_index]
	var result: Dictionary = {}
	for reaction_index in range(reaction_count):
		result[compiled_network.reaction_ids[reaction_index]] = cumulative_flux[reaction_index]

	_store_buffer(workspace, "pool_values", pool_values)
	_store_buffer(workspace, "snapshot_values", snapshot_values)
	_store_buffer(workspace, "substrate_demand", substrate_demand)
	_store_buffer(workspace, "substrate_scale", substrate_scale)
	_store_buffer(workspace, "deltas", deltas)
	_store_buffer(workspace, "catalytic_activities", catalytic_activities)
	_store_buffer(workspace, "potential_flux", potential_flux)
	_store_buffer(workspace, "cumulative_flux", cumulative_flux)
	return result

static func _take_buffer(workspace: Dictionary, key: String, size: int, clear: bool = false) -> PackedFloat64Array:
	var buffer: PackedFloat64Array = workspace.get(key, PackedFloat64Array())
	workspace.erase(key)
	if buffer.size() != size:
		buffer.resize(size)
	if clear:
		buffer.fill(0.0)
	return buffer

static func _store_buffer(workspace: Dictionary, key: String, buffer: PackedFloat64Array) -> void:
	workspace[key] = buffer
