extends RefCounted
class_name CellState

const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

var id: int
var parent_id: int
var generation: int
var birth_tick: int
var position: Vector2
var genome = null

var alive: bool = true
var death_reason: String = ""

# Authoritative molecular state.
var metabolites: Dictionary = {}
var expression_state: Dictionary = {}
var last_fluxes: Dictionary = {}
var last_expression_summary: Dictionary = {}
var volume: float = 1.0
var damage: float = 0.0
var energy_debt: float = 0.0

func _init(
	p_id: int = 0,
	p_parent_id: int = -1,
	p_generation: int = 0,
	p_birth_tick: int = 0,
	p_position: Vector2 = Vector2.ZERO,
	p_volume: float = 1.0
) -> void:
	id = p_id
	parent_id = p_parent_id
	generation = p_generation
	birth_tick = p_birth_tick
	position = p_position
	volume = p_volume

func initialize_molecular_state(config) -> void:
	assert(genome != null, "Genome must exist before molecular initialization")
	metabolites = MetabolicSolverScript.create_initial_pools(volume, config)
	expression_state = ExpressionSystemScript.create_equilibrium_state(genome, config)
	# The equilibrium constructor predates the finite M5-C proteome budget and
	# expresses every locus independently. Normalize the initial condition into
	# the physically admissible proteome without recycling: initialization is a
	# state construction, not a biochemical degradation event.
	_enforce_proteome_budget(config, false)
	last_fluxes = {}
	last_expression_summary = {}
	_sync_volume_from_biomass(config)

# Compatibility alias for setup code introduced in M4.
func initialize_metabolism(config) -> void:
	initialize_molecular_state(config)

func pool(metabolite_id: String) -> float:
	return float(metabolites.get(metabolite_id, 0.0))

func set_pool(metabolite_id: String, amount: float) -> void:
	assert(amount >= 0.0)
	metabolites[metabolite_id] = amount

func total_mrna() -> float:
	return ExpressionSystemScript.total_mrna(expression_state)

func total_protein() -> float:
	return ExpressionSystemScript.total_protein(expression_state)

func proteome_capacity(config) -> float:
	return float(config.proteome_capacity_reference_units) * float(config.expression_reference_protein_count)

func proteome_utilization(config) -> float:
	return total_protein() / maxf(1e-12, proteome_capacity(config))

func transport_requests(dt: float, world, config) -> Dictionary:
	if not alive:
		return {"glucose": 0.0, "oxygen": 0.0, "nitrogen": 0.0, "phosphorus": 0.0}
	assert(not metabolites.is_empty(), "Cell metabolism must be initialized before transport")
	return {
		"glucose": _transport_request("glucose", "G", float(config.glucose_transport_vmax), float(config.glucose_transport_km), dt, world, config),
		"oxygen": _transport_request("oxygen", "O2", float(config.oxygen_transport_vmax), float(config.oxygen_transport_km), dt, world, config),
		"nitrogen": _transport_request("nitrogen", "NH4", float(config.nitrogen_transport_vmax), float(config.nitrogen_transport_km), dt, world, config),
		"phosphorus": _transport_request("phosphorus", "P", float(config.phosphorus_transport_vmax), float(config.phosphorus_transport_km), dt, world, config)
	}

func _transport_request(world_field: String, internal_id: String, vmax: float, km: float, dt: float, world, config) -> float:
	var local_amount: float = maxf(0.0, float(world.sample(world_field, position)))
	var rate: float = vmax * local_amount / (km + local_amount)
	var capacity: float = maxf(0.0, float(config.intracellular_pool_capacity_per_volume) * volume - pool(internal_id))
	return minf(rate * dt, capacity)

func apply_uptake(uptake: Dictionary) -> void:
	assert(not metabolites.is_empty())
	var mapping: Dictionary = {
		"glucose": "G",
		"oxygen": "O2",
		"nitrogen": "NH4",
		"phosphorus": "P"
	}
	for world_field in mapping.keys():
		var amount: float = maxf(0.0, float(uptake.get(world_field, 0.0)))
		MetabolicSolverScript.add_pool(metabolites, String(mapping[world_field]), amount)

func step_intracellular(dt: float, config, reactions: Array, rng) -> void:
	if not alive:
		return
	assert(genome != null, "M5 physiology requires a genome")
	assert(not metabolites.is_empty(), "M5 physiology requires initialized metabolite pools")
	assert(not expression_state.is_empty(), "M5 physiology requires explicit expression state")

	# Expression happens before metabolism for this tick. It consumes current
	# ATP/material and changes the proteome; metabolism then reads that realized
	# proteome. No promoter value enters catalytic flux directly.
	last_expression_summary = ExpressionSystemScript.step(expression_state, genome, metabolites, dt, rng, config)
	var proteome_summary: Dictionary = _enforce_proteome_budget(config, true)
	for key in proteome_summary.keys():
		last_expression_summary[key] = proteome_summary[key]
	last_fluxes = MetabolicSolverScript.step(metabolites, genome, expression_state, reactions, dt, volume, config)
	_sync_volume_from_biomass(config)
	_pay_maintenance(dt, config)
	_update_damage_and_repair(dt, config)
	_check_viability(config)
	_assert_state(config)

# Generic finite proteome constraint. All realized protein cohorts compete for
# the same physical budget; no locus, sequence or biological function is given
# priority. If stochastic expression overshoots the budget, every cohort is
# scaled by the same factor. Removed protein returns only its modeled amino-acid
# material. ATP spent to synthesize it remains spent, so unnecessary expression
# carries an energetic opportunity cost instead of being free.
func _enforce_proteome_budget(config, recycle_excess: bool) -> Dictionary:
	var capacity: float = proteome_capacity(config)
	var before: float = total_protein()
	var scale: float = 1.0
	var removed: float = 0.0
	if before > capacity and before > 0.0:
		scale = capacity / before
		var loci: Array = expression_state.keys()
		loci.sort()
		for locus_variant in loci:
			var locus_id: int = int(locus_variant)
			var cohorts: Dictionary = expression_state[locus_id]["protein"]
			var signatures: Array = cohorts.keys()
			signatures.sort()
			for signature_variant in signatures:
				var old_amount: float = maxf(0.0, float(cohorts[signature_variant]))
				var new_amount: float = old_amount * scale
				cohorts[signature_variant] = new_amount
				removed += old_amount - new_amount
	if recycle_excess and removed > 0.0:
		metabolites["AA"] = float(metabolites.get("AA", 0.0)) + removed * float(config.translation_aa_cost_per_event)
	return {
		"proteome_capacity": capacity,
		"proteome_before_constraint": before,
		"proteome_removed": removed,
		"proteome_scale": scale,
		"proteome_utilization": total_protein() / maxf(1e-12, capacity)
	}

func _pay_maintenance(dt: float, config) -> void:
	var required: float = float(config.maintenance_atp_rate_per_volume) * volume * dt
	var paid: float = MetabolicSolverScript.spend_atp(metabolites, required)
	var unmet: float = required - paid
	if unmet > 0.0:
		energy_debt += unmet
	else:
		energy_debt = maxf(0.0, energy_debt - required * 0.25)

func _update_damage_and_repair(dt: float, config) -> void:
	var current_ros: float = pool("ROS")
	var spontaneous_decay: float = minf(current_ros, current_ros * float(config.spontaneous_ros_decay_rate) * dt)
	metabolites["ROS"] = current_ros - spontaneous_decay
	damage += float(config.ros_damage_rate) * pool("ROS") * dt

	var possible_repair: float = minf(damage, float(config.basal_repair_rate) * dt)
	var requested_cost: float = possible_repair * float(config.repair_atp_cost)
	var paid: float = MetabolicSolverScript.spend_atp(metabolites, requested_cost)
	if requested_cost > 0.0:
		possible_repair *= paid / requested_cost
	damage = maxf(0.0, damage - possible_repair)

func _sync_volume_from_biomass(config) -> void:
	volume = pool("BIO") / float(config.biomass_units_per_volume)
	assert(volume > 0.0, "Living cell cannot have zero structural biomass")

func _check_viability(config) -> void:
	if damage >= float(config.lethal_damage):
		alive = false
		death_reason = "damage"
	elif energy_debt >= float(config.lethal_energy_debt):
		alive = false
		death_reason = "energy_failure"

func ready_to_divide(config) -> bool:
	return alive and volume >= float(config.division_volume) and pool("ATP") >= float(config.division_atp_cost)

func create_daughters(first_id: int, second_id: int, tick: int, rng, world, config) -> Array:
	assert(ready_to_divide(config))
	MetabolicSolverScript.spend_atp(metabolites, float(config.division_atp_cost))

	var partition_jitter: float = float(config.partition_jitter)
	var ratio: float = 0.5 + float(rng.randf_range(-partition_jitter, partition_jitter))
	var offset_scale: float = float(config.daughter_offset_grid)
	var first_offset: Vector2 = Vector2(float(rng.randf_range(-1.0, 1.0)), float(rng.randf_range(-1.0, 1.0))).normalized() * offset_scale
	if first_offset == Vector2.ZERO:
		first_offset = Vector2(offset_scale, 0.0)
	var second_offset: Vector2 = -first_offset

	var metabolite_partitions: Array = MetabolicSolverScript.partition(metabolites, ratio)
	var expression_partitions: Array = ExpressionSystemScript.partition(expression_state, ratio, rng, config)
	var first = CellState.new(first_id, id, generation + 1, tick, world.clamp_position(position + first_offset), volume * ratio)
	var second = CellState.new(second_id, id, generation + 1, tick, world.clamp_position(position + second_offset), volume * (1.0 - ratio))
	first.metabolites = metabolite_partitions[0]
	second.metabolites = metabolite_partitions[1]
	first.expression_state = expression_partitions[0]
	second.expression_state = expression_partitions[1]
	first.damage = damage * ratio
	second.damage = damage * (1.0 - ratio)
	first.energy_debt = energy_debt * ratio
	second.energy_debt = energy_debt * (1.0 - ratio)
	first.genome = genome.deep_copy() if genome != null else null
	second.genome = genome.deep_copy() if genome != null else null
	first._sync_volume_from_biomass(config)
	second._sync_volume_from_biomass(config)
	alive = false
	death_reason = "division"
	return [first, second]

func releasable_pools() -> Dictionary:
	return {
		"glucose": maxf(0.0, pool("G")),
		"oxygen": maxf(0.0, pool("O2")),
		"nitrogen": maxf(0.0, pool("NH4")),
		"phosphorus": maxf(0.0, pool("P"))
	}

func total_adenylate() -> float:
	return pool("ATP") + pool("ADP")

func total_redox_currency() -> float:
	return pool("NAD") + pool("NADH")

func _assert_state(config) -> void:
	assert(volume > 0.0)
	assert(damage >= -1e-10)
	assert(energy_debt >= -1e-10)
	MetabolicSolverScript.assert_nonnegative(metabolites)
	assert(absf(volume - pool("BIO") / float(config.biomass_units_per_volume)) <= 1e-10)
	if genome != null:
		genome.validate()
	assert(total_mrna() >= -1e-10 and total_protein() >= -1e-10)
	assert(total_protein() <= proteome_capacity(config) + 1e-8)

func checksum() -> float:
	var result: float = (
		float(id) * 0.001
		+ position.x * 3.0
		+ position.y * 5.0
		+ volume * 7.0
		+ damage * 29.0
		+ energy_debt * 31.0
		+ MetabolicSolverScript.checksum(metabolites) * 0.001
		+ ExpressionSystemScript.checksum(expression_state) * 0.017
	)
	if genome != null:
		result += float(genome.checksum()) * 0.013
	return result