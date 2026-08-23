extends RefCounted
class_name SimulationEngine

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const ExtracellularReactionCatalogScript = preload("res://src/chemistry/extracellular_reaction_catalog.gd")
const ExtracellularCatalysisScript = preload("res://src/chemistry/extracellular_catalysis.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const CellMechanicsScript = preload("res://src/physics/cell_mechanics.gd")

# Historical basal uptake is preserved separately while M7 adds sequence-derived
# exchange for the secondary extracellular chemistry.
const TRANSPORTED_RESOURCES: Array[String] = ["glucose", "oxygen", "nitrogen", "phosphorus"]

var config
var rng
var world
var mutation_engine
var reactions: Array = []
var extracellular_reactions: Array = []
var cells: Array = []
var tick_index: int = 0
var simulation_time_min: float = 0.0
var next_cell_id: int = 1
var next_mutation_id: int = 1
var event_log: Array = []
var last_mechanics_summary: Dictionary = {}
var last_secondary_transport_summary: Dictionary = {}
var last_protein_secretion_summary: Dictionary = {}
var last_extracellular_catalysis_summary: Dictionary = {}

func _init(p_config = null) -> void:
	config = p_config if p_config != null else SimConfigScript.new()
	config.validate()
	rng = DeterministicRngScript.new(config.seed)
	mutation_engine = MutationEngineScript.new()
	reactions = ReactionCatalogScript.create_m4_candidate()
	ReactionCatalogScript.validate_unique(reactions)
	extracellular_reactions = ExtracellularReactionCatalogScript.create_m7_candidate()
	ExtracellularReactionCatalogScript.validate_unique(extracellular_reactions)
	world = WorldStateScript.new(config.world_width, config.world_height, config.grid_cell_size_um)
	_register_extracellular_fields()

func _register_extracellular_fields() -> void:
	# Preserve the historical field order for the four basal resources so old
	# diagnostics remain easy to compare, then append every new field in sorted
	# molecular-ID order for deterministic checksums/replay.
	world.register_field("glucose", config.glucose_diffusion, config.initial_glucose)
	world.register_field("oxygen", config.oxygen_diffusion, config.initial_oxygen)
	world.register_field("nitrogen", config.nitrogen_diffusion, config.initial_nitrogen)
	world.register_field("phosphorus", config.phosphorus_diffusion, config.initial_phosphorus)

	for metabolite_id in MetaboliteCatalogScript.extracellular_ids():
		var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
		if world.has_field(field_name):
			continue
		world.register_field(
			field_name,
			config.extracellular_diffusion_coefficient(metabolite_id),
			config.extracellular_initial_amount(metabolite_id)
		)

func seed_ancestor(position: Vector2 = Vector2(-1.0, -1.0)):
	if position.x < 0.0 or position.y < 0.0:
		position = Vector2(float(config.world_width - 1) * 0.5, float(config.world_height - 1) * 0.5)
	var cell = CellStateScript.new(_allocate_cell_id(), -1, 0, tick_index, world.clamp_position(position), config.ancestor_volume)
	cell.genome = GenomeScript.create_ancestor()
	cell.initialize_molecular_state(config)
	cell.position = CellMechanicsScript.clamp_position(cell.position, CellMechanicsScript.radius_for_cell(cell, config), world)
	cells.append(cell)
	_record_event("birth", {
		"cell_id": cell.id,
		"parent_id": -1,
		"generation": 0,
		"genotype_fingerprint": cell.genome.fingerprint()
	})
	return cell

func step(tick_count: int = 1) -> void:
	assert(tick_count >= 0)
	for _unused in range(tick_count):
		_step_once()

func _step_once() -> void:
	var dt: float = float(config.tick_dt_min)
	world.diffuse(dt)
	_allocate_membrane_transport(dt)
	last_secondary_transport_summary = _allocate_secondary_membrane_transport(dt)
	last_protein_secretion_summary = _secrete_extracellular_proteins(dt)
	last_extracellular_catalysis_summary = ExtracellularCatalysisScript.step(world, extracellular_reactions, dt, config)

	# Extracellular products are created after membrane allocation, so another
	# cell cannot consume a just-created public molecule at zero elapsed time.
	# It becomes available to transport on the next simulation tick after the
	# ordinary diffusion phase has had a chance to establish spatial structure.
	for cell in cells:
		if cell.alive:
			cell.step_intracellular(dt, config, reactions, rng)

	_process_deaths()
	_process_divisions()
	last_mechanics_summary = CellMechanicsScript.relax(cells, world, config, bool(config.mechanical_use_spatial_index))
	world.assert_nonnegative()
	tick_index += 1
	simulation_time_min += dt

func _allocate_membrane_transport(dt: float) -> void:
	var records: Array = []
	var totals: Dictionary = {}
	for resource in TRANSPORTED_RESOURCES:
		totals[resource] = {}

	for cell in cells:
		if not cell.alive:
			continue
		var request: Dictionary = cell.transport_requests(dt, world, config)
		var key: Vector2i = _grid_key(cell.position)
		records.append({"cell": cell, "key": key, "request": request})
		for resource in TRANSPORTED_RESOURCES:
			var resource_totals: Dictionary = totals[resource]
			resource_totals[key] = float(resource_totals.get(key, 0.0)) + float(request.get(resource, 0.0))
			totals[resource] = resource_totals

	var scales: Dictionary = {}
	for resource in TRANSPORTED_RESOURCES:
		var resource_scales: Dictionary = {}
		var resource_totals: Dictionary = totals[resource]
		for key_variant in resource_totals.keys():
			var key: Vector2i = key_variant
			var available: float = float(world.get_field(resource).get_value(key.x, key.y))
			var total_request: float = float(resource_totals[key])
			resource_scales[key] = 1.0 if total_request <= available or total_request <= 0.0 else available / total_request
		scales[resource] = resource_scales

	for record in records:
		var cell = record["cell"]
		var key: Vector2i = record["key"]
		var request: Dictionary = record["request"]
		var uptake: Dictionary = {}
		for resource in TRANSPORTED_RESOURCES:
			var resource_scales: Dictionary = scales[resource]
			var allocated: float = float(request.get(resource, 0.0)) * float(resource_scales.get(key, 1.0))
			uptake[resource] = float(world.get_field(resource).remove_amount(key.x, key.y, allocated))
		cell.apply_uptake(uptake)

# M7-B generic secondary membrane transport. Every proposal is computed from a
# common pre-exchange world/cell state. Positive exchange is import and negative
# exchange is export. ATP scarcity scales all proposals from a given cell by the
# same factor, while extracellular scarcity scales all importers at one lattice
# site proportionally. Exports become available only after that allocation, so
# a producer cannot feed a consumer instantaneously within the same transport
# phase.
func _allocate_secondary_membrane_transport(dt: float) -> Dictionary:
	var records: Array = []
	var import_totals: Dictionary = {}
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		import_totals[metabolite_id] = {}

	for cell in cells:
		if not cell.alive:
			continue
		var key: Vector2i = _grid_key(cell.position)
		var desired: Dictionary = {}
		var activities: Dictionary = {}
		for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
			var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
			var external_amount: float = maxf(0.0, float(world.get_field(field_name).get_value(key.x, key.y)))
			var activity: float = MembraneTransportScript.proteome_activity(cell.expression_state, metabolite_id, config)
			activities[metabolite_id] = activity
			desired[metabolite_id] = MembraneTransportScript.desired_exchange(
				cell.pool(metabolite_id),
				cell.volume,
				external_amount,
				activity,
				dt,
				config
			)

		var energy_scale: float = MembraneTransportScript.energy_scale(desired, cell.pool("ATP"), config)
		var proposed: Dictionary = {}
		for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
			var signed_amount: float = float(desired[metabolite_id]) * energy_scale
			proposed[metabolite_id] = signed_amount
			if signed_amount > 0.0:
				var totals_by_key: Dictionary = import_totals[metabolite_id]
				totals_by_key[key] = float(totals_by_key.get(key, 0.0)) + signed_amount
				import_totals[metabolite_id] = totals_by_key
		records.append({
			"cell": cell,
			"key": key,
			"activities": activities,
			"energy_scale": energy_scale,
			"proposed": proposed
		})

	var import_scales: Dictionary = {}
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
		var scales_by_key: Dictionary = {}
		var totals_by_key: Dictionary = import_totals[metabolite_id]
		for key_variant in totals_by_key.keys():
			var key: Vector2i = key_variant
			var requested: float = float(totals_by_key[key])
			var available: float = maxf(0.0, float(world.get_field(field_name).get_value(key.x, key.y)))
			scales_by_key[key] = 1.0 if requested <= available or requested <= 0.0 else available / requested
		import_scales[metabolite_id] = scales_by_key

	var world_imports: Dictionary = {}
	var world_exports: Dictionary = {}
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		world_imports[metabolite_id] = {}
		world_exports[metabolite_id] = {}

	for record in records:
		var key: Vector2i = record["key"]
		var proposed: Dictionary = record["proposed"]
		var actual: Dictionary = {}
		for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
			var signed_amount: float = float(proposed[metabolite_id])
			if signed_amount > 0.0:
				var scale: float = float(import_scales[metabolite_id].get(key, 1.0))
				var imported: float = signed_amount * scale
				actual[metabolite_id] = imported
				var imports_by_key: Dictionary = world_imports[metabolite_id]
				imports_by_key[key] = float(imports_by_key.get(key, 0.0)) + imported
				world_imports[metabolite_id] = imports_by_key
			elif signed_amount < 0.0:
				actual[metabolite_id] = signed_amount
				var exports_by_key: Dictionary = world_exports[metabolite_id]
				exports_by_key[key] = float(exports_by_key.get(key, 0.0)) + absf(signed_amount)
				world_exports[metabolite_id] = exports_by_key
			else:
				actual[metabolite_id] = 0.0
		record["actual"] = actual

	# World deltas are applied in aggregate per field/lattice site, eliminating
	# sequential first-come access from the authoritative resource update.
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
		var field = world.get_field(field_name)
		var imports_by_key: Dictionary = world_imports[metabolite_id]
		for key_variant in imports_by_key.keys():
			var key: Vector2i = key_variant
			var requested_remove: float = float(imports_by_key[key])
			var removed: float = float(field.remove_amount(key.x, key.y, requested_remove))
			assert(absf(removed - requested_remove) <= 1e-9, "Secondary import allocation exceeded snapshot availability")
		var exports_by_key: Dictionary = world_exports[metabolite_id]
		for key_variant in exports_by_key.keys():
			var key: Vector2i = key_variant
			field.add_amount(key.x, key.y, float(exports_by_key[key]))

	var by_cell: Dictionary = {}
	var total_moved: float = 0.0
	var total_atp_spent: float = 0.0
	for record in records:
		var cell = record["cell"]
		var actual: Dictionary = record["actual"]
		var moved: float = 0.0
		var signed_actual: Dictionary = {}
		for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
			var signed_amount: float = float(actual[metabolite_id])
			signed_actual[metabolite_id] = signed_amount
			if signed_amount > 0.0:
				MetabolicSolverScript.add_pool(cell.metabolites, metabolite_id, signed_amount)
			elif signed_amount < 0.0:
				cell.set_pool(metabolite_id, maxf(0.0, cell.pool(metabolite_id) + signed_amount))
			moved += absf(signed_amount)
		var cost: float = MembraneTransportScript.movement_cost(moved, config)
		var spent: float = MetabolicSolverScript.spend_atp(cell.metabolites, cost)
		assert(absf(spent - cost) <= 1e-9, "ATP pre-scaling failed to fund secondary transport")
		total_moved += moved
		total_atp_spent += spent
		by_cell[int(cell.id)] = {
			"exchange": signed_actual,
			"activities": record["activities"],
			"energy_scale": record["energy_scale"],
			"moved": moved,
			"atp_spent": spent
		}

	return {
		"by_cell": by_cell,
		"total_moved": total_moved,
		"total_atp_spent": total_atp_spent
	}

# M7-E protein secretion is sequence-derived and acts on realized protein
# cohorts. Every actual secreted molecule is removed from the intracellular
# proteome, deposited into an extracellular sequence-specific field, and paid
# for with ATP. ATP scarcity scales all simultaneously secreted signatures from
# one cell proportionally.
func _secrete_extracellular_proteins(dt: float) -> Dictionary:
	var by_cell: Dictionary = {}
	var total_secreted: float = 0.0
	var total_atp_spent: float = 0.0

	for cell in cells:
		if not cell.alive:
			continue
		var proposals: Dictionary = ExtracellularCatalysisScript.secretion_proposals(cell.expression_state, dt, config)
		var energy_scale: float = ExtracellularCatalysisScript.secretion_energy_scale(proposals, cell.pool("ATP"), config)
		var actual_by_signature: Dictionary = {}
		var secreted: float = 0.0
		var signatures: Array = proposals.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			var requested: float = maxf(0.0, float(proposals[signature_variant]) * energy_scale)
			var removed: float = ExtracellularCatalysisScript.remove_protein_signature(cell.expression_state, signature, requested)
			assert(absf(removed - requested) <= 1e-9, "Protein secretion proposal exceeded realized cohort")
			if removed > 0.0:
				world.release_protein(signature, cell.position, removed, float(config.extracellular_protein_diffusion))
			actual_by_signature[signature] = removed
			secreted += removed

		var cost: float = secreted * float(config.extracellular_protein_secretion_atp_cost_per_unit)
		var spent: float = MetabolicSolverScript.spend_atp(cell.metabolites, cost)
		assert(absf(spent - cost) <= 1e-9, "ATP pre-scaling failed to fund protein secretion")
		total_secreted += secreted
		total_atp_spent += spent
		by_cell[int(cell.id)] = {
			"secreted": actual_by_signature,
			"total_secreted": secreted,
			"energy_scale": energy_scale,
			"atp_spent": spent
		}

	return {
		"by_cell": by_cell,
		"total_secreted": total_secreted,
		"total_atp_spent": total_atp_spent
	}

func _process_deaths() -> void:
	var survivors: Array = []
	for cell in cells:
		if cell.alive:
			survivors.append(cell)
			continue
		var pools: Dictionary = cell.releasable_pools(config)
		var released: Dictionary = {}
		var field_names: Array = pools.keys()
		field_names.sort()
		for field_variant in field_names:
			var field_name := String(field_variant)
			assert(world.has_field(field_name), "Lysis target has no extracellular field: %s" % field_name)
			var amount: float = maxf(0.0, float(pools[field_name]))
			if amount <= 0.0:
				continue
			world.release(field_name, cell.position, amount)
			released[field_name] = amount
		_record_event("death", {
			"cell_id": cell.id,
			"generation": cell.generation,
			"reason": cell.death_reason,
			"genotype_fingerprint": cell.genome.fingerprint() if cell.genome != null else -1,
			"released_pools": released
		})
	cells = survivors

func _process_divisions() -> void:
	var next_population: Array = []
	var projected_population: int = cells.size()
	for cell in cells:
		if not cell.alive:
			continue
		if cell.ready_to_divide(config) and projected_population < int(config.max_cells):
			var daughters: Array = cell.create_daughters(_allocate_cell_id(), _allocate_cell_id(), tick_index, rng, world, config)
			projected_population += 1
			_record_event("division", {
				"parent_id": cell.id,
				"parent_genotype_fingerprint": cell.genome.fingerprint() if cell.genome != null else -1,
				"daughter_ids": [daughters[0].id, daughters[1].id],
				"generation": cell.generation + 1
			})
			for daughter in daughters:
				var parent_fingerprint: int = int(cell.genome.fingerprint())
				var mutation_result: Dictionary = mutation_engine.mutate_copy(daughter.genome, rng, config)
				daughter.genome = mutation_result["genome"]
				var daughter_fingerprint: int = int(daughter.genome.fingerprint())
				_record_event("birth", {
					"cell_id": daughter.id,
					"parent_id": cell.id,
					"generation": daughter.generation,
					"genotype_fingerprint": daughter_fingerprint
				})
				for raw_event in mutation_result["events"]:
					var mutation_payload: Dictionary = raw_event.duplicate(true)
					mutation_payload["mutation_id"] = _allocate_mutation_id()
					mutation_payload["cell_id"] = daughter.id
					mutation_payload["parent_id"] = cell.id
					mutation_payload["generation"] = daughter.generation
					mutation_payload["parent_genotype_fingerprint"] = parent_fingerprint
					mutation_payload["resulting_genotype_fingerprint"] = daughter_fingerprint
					_record_event("mutation", mutation_payload)
			next_population.append_array(daughters)
		else:
			next_population.append(cell)
	cells = next_population

func relax_mechanics() -> Dictionary:
	last_mechanics_summary = CellMechanicsScript.relax(cells, world, config, bool(config.mechanical_use_spatial_index))
	return last_mechanics_summary

func maximum_overlap() -> float:
	return CellMechanicsScript.max_overlap(cells, config)

func _grid_key(position: Vector2) -> Vector2i:
	return Vector2i(
		clampi(roundi(position.x), 0, int(config.world_width) - 1),
		clampi(roundi(position.y), 0, int(config.world_height) - 1)
	)

func _allocate_cell_id() -> int:
	var result: int = next_cell_id
	next_cell_id += 1
	return result

func _allocate_mutation_id() -> int:
	var result: int = next_mutation_id
	next_mutation_id += 1
	return result

func _record_event(kind: String, payload: Dictionary) -> void:
	var event: Dictionary = payload.duplicate(true)
	event["kind"] = kind
	event["tick"] = tick_index
	event["time_min"] = simulation_time_min
	event_log.append(event)

func population_size() -> int:
	return cells.size()

func maximum_generation() -> int:
	var result: int = 0
	for cell in cells:
		result = maxi(result, int(cell.generation))
	return result

func total_cell_volume() -> float:
	var result: float = 0.0
	for cell in cells:
		result += float(cell.volume)
	return result

func mutation_event_count() -> int:
	var result: int = 0
	for event in event_log:
		if event["kind"] == "mutation":
			result += 1
	return result

func genotype_count() -> int:
	var genotypes: Dictionary = {}
	for cell in cells:
		if cell.genome != null:
			genotypes[cell.genome.canonical_key()] = true
	return genotypes.size()

func checksum() -> float:
	var result: float = float(world.checksum()) + float(tick_index) * 37.0 + simulation_time_min * 41.0
	for cell in cells:
		result += float(cell.checksum())
	result += float(rng.get_state() % 1000003) * 1e-6
	result += float(next_cell_id) * 0.00017 + float(next_mutation_id) * 0.00019
	return result
