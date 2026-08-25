extends RefCounted
class_name CellState

const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const DNAReplicationScript = preload("res://src/genetics/dna_replication.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const SignallingSystemScript = preload("res://src/signalling/signalling_system.gd")

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
var last_replication_summary: Dictionary = {}
# M11 reversible activation molecules and physical orientation. Neither is a
# behavioral label: signalling_state is keyed only by molecular signature and
# motor_heading is ordinary geometry required by any active force generator.
var signalling_state: Dictionary = {}
var motor_heading: Vector2 = Vector2.RIGHT
var last_receptor_summary: Dictionary = {}
var last_signalling_summary: Dictionary = {}
var last_motility_summary: Dictionary = {}
var volume: float = 1.0
var damage: float = 0.0
var energy_debt: float = 0.0

# M10 authoritative DNA-copying state. Progress is physical work completed on
# the current cell cycle, not a fitness score. Daughters begin a new cycle at 0.
var replication_progress: float = 0.0
var replication_gene_equivalents_copied: float = 0.0
var replication_atp_spent: float = 0.0
var replication_nuc_spent: float = 0.0
var replication_repair_activity_integral: float = 0.0

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
	expression_state = ExpressionSystemScript.create_equilibrium_state(genome, config, volume)
	# The equilibrium constructor predates the finite M5-C proteome budget and
	# expresses every locus independently. Normalize the initial condition into
	# the physically admissible proteome without recycling: initialization is a
	# state construction, not a biochemical degradation event.
	_enforce_proteome_budget(config, false)
	last_fluxes = {}
	last_expression_summary = {}
	last_replication_summary = {}
	signalling_state = {}
	motor_heading = Vector2.RIGHT
	last_receptor_summary = {}
	last_signalling_summary = {}
	last_motility_summary = {}
	_sync_volume_from_biomass(config)
	_reset_replication_cycle()
	# Controlled fixtures historically construct already division-sized cells.
	# Such a constructed state represents a mother whose genome is already
	# copied; ordinary simulations seed at ancestor_volume and must replicate.
	if bool(config.evolvable_replication_enabled) and volume >= float(config.division_volume):
		replication_gene_equivalents_copied = float(genome.gene_count())
		replication_progress = 1.0

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

# Proteome is an extensive physical inventory. The historical M5 reference cap
# applies at ancestor volume; a cell with twice that biomass can physically hold
# twice as much protein without changing its proteome concentration.
func proteome_capacity(config) -> float:
	var reference_capacity: float = float(config.proteome_capacity_reference_units) * float(config.expression_reference_protein_count)
	var volume_scale: float = volume / maxf(1e-12, float(config.ancestor_volume))
	return reference_capacity * maxf(1e-12, volume_scale)

func proteome_utilization(config) -> float:
	return total_protein() / maxf(1e-12, proteome_capacity(config))

# M11 basal membrane migration. The historical chemical Vmax/Km values remain
# reference kinetic scales, but actual uptake capacity now comes from proteins
# that physically exist. DNA mutation alone therefore cannot change uptake until
# expression/turnover changes the realized membrane-compatible cohort.
func transport_requests(dt: float, world, config) -> Dictionary:
	if not alive:
		return {"glucose": 0.0, "oxygen": 0.0, "nitrogen": 0.0, "phosphorus": 0.0}
	assert(not metabolites.is_empty(), "Cell metabolism must be initialized before transport")
	return {
		"glucose": _molecular_transport_request("glucose", "G", float(config.glucose_transport_vmax), float(config.glucose_transport_km), dt, world, config),
		"oxygen": _molecular_transport_request("oxygen", "O2", float(config.oxygen_transport_vmax), float(config.oxygen_transport_km), dt, world, config),
		"nitrogen": _molecular_transport_request("nitrogen", "NH4", float(config.nitrogen_transport_vmax), float(config.nitrogen_transport_km), dt, world, config),
		"phosphorus": _molecular_transport_request("phosphorus", "P", float(config.phosphorus_transport_vmax), float(config.phosphorus_transport_km), dt, world, config)
	}

func _molecular_transport_request(world_field: String, internal_id: String, reference_vmax: float, km: float, dt: float, world, config) -> float:
	var local_amount: float = maxf(0.0, float(world.sample(world_field, position)))
	var activity: float = MembraneTransportScript.basal_proteome_activity(expression_state, internal_id, config)
	return MembraneTransportScript.basal_import_request(
		pool(internal_id), volume, local_amount, activity, internal_id,
		reference_vmax, km, dt, config
	)

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
	last_expression_summary = ExpressionSystemScript.step(expression_state, genome, metabolites, dt, rng, config, volume)
	var proteome_summary: Dictionary = _enforce_proteome_budget(config, true)
	for key in proteome_summary.keys():
		last_expression_summary[key] = proteome_summary[key]
	last_fluxes = MetabolicSolverScript.step(metabolites, genome, expression_state, reactions, dt, volume, config)
	_sync_volume_from_biomass(config)
	# DNA copying occurs after this tick's metabolism, so the replication fork
	# pays from explicit ATP/NUC actually available in the cell. This cost then
	# competes with maintenance and repair rather than being an abstract delay.
	last_replication_summary = DNAReplicationScript.step(self, dt, config)
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

# ATP/ADP and NAD/NADH are explicit energetic/redox carrier currencies but are
# deliberately outside the model's structural C/N/P bookkeeping. Their charge
# state is physical (metabolism must convert ADP->ATP and NAD->NADH), while the
# carrier scaffold is implicit cellular material. New biomass therefore adds
# only discharged/oxidized carrier capacity. The previous fixed founder-wide
# carrier inventory was an artificial population ceiling: after enough divisions
# no cell could hold the ATP required for another division even with nutrients.
func _sync_volume_from_biomass(config) -> void:
	var previous_volume: float = volume
	var next_volume: float = pool("BIO") / float(config.biomass_units_per_volume)
	assert(next_volume > 0.0, "Living cell cannot have zero structural biomass")
	var added_volume: float = maxf(0.0, next_volume - previous_volume)
	if added_volume > 0.0:
		var adenylate_capacity_per_volume: float = (
			float(config.initial_atp_per_volume) + float(config.initial_adp_per_volume)
		)
		var redox_capacity_per_volume: float = (
			float(config.initial_nad_per_volume) + float(config.initial_nadh_per_volume)
		)
		# Capacity is born uncharged: no ATP or NADH is gifted by cell growth.
		metabolites["ADP"] = float(metabolites.get("ADP", 0.0)) + added_volume * adenylate_capacity_per_volume
		metabolites["NAD"] = float(metabolites.get("NAD", 0.0)) + added_volume * redox_capacity_per_volume
	volume = next_volume

func _check_viability(config) -> void:
	if damage >= float(config.lethal_damage):
		alive = false
		death_reason = "damage"
	elif energy_debt >= float(config.lethal_energy_debt):
		alive = false
		death_reason = "energy_failure"

func ready_to_divide(config) -> bool:
	var replication_ready: bool = (
		not bool(config.evolvable_replication_enabled)
		or DNAReplicationScript.replication_complete(self)
	)
	return alive and replication_ready and volume >= float(config.division_volume) and pool("ATP") >= float(config.division_atp_cost)

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
	var signalling_partitions: Array = SignallingSystemScript.partition(signalling_state, ratio)
	var first = CellState.new(first_id, id, generation + 1, tick, world.clamp_position(position + first_offset), volume * ratio)
	var second = CellState.new(second_id, id, generation + 1, tick, world.clamp_position(position + second_offset), volume * (1.0 - ratio))
	first.metabolites = metabolite_partitions[0]
	second.metabolites = metabolite_partitions[1]
	first.expression_state = expression_partitions[0]
	second.expression_state = expression_partitions[1]
	first.signalling_state = signalling_partitions[0]
	second.signalling_state = signalling_partitions[1]
	first.motor_heading = motor_heading
	second.motor_heading = motor_heading
	first.damage = damage * ratio
	second.damage = damage * (1.0 - ratio)
	first.energy_debt = energy_debt * ratio
	second.energy_debt = energy_debt * (1.0 - ratio)
	first.genome = genome.deep_copy() if genome != null else null
	second.genome = genome.deep_copy() if genome != null else null
	first._sync_volume_from_biomass(config)
	second._sync_volume_from_biomass(config)
	# Daughter constructors already initialize a fresh replication cycle at zero;
	# copied parental DNA has been partitioned into one inherited genome per cell.
	alive = false
	death_reason = "division"
	return [first, second]

func _reset_replication_cycle() -> void:
	replication_progress = 0.0
	replication_gene_equivalents_copied = 0.0
	replication_atp_spent = 0.0
	replication_nuc_spent = 0.0
	replication_repair_activity_integral = 0.0
	last_replication_summary = {}

# Lysis exposes every metabolite that has a physical extracellular field. BIO
# is not discarded: it is hydrolyzed into the exact precursor stoichiometry of
# the reverse structural-biomass reaction (2 AA + 1 LIP + 2 NUC per BIO). The
# material physically stored in mRNA/protein and M10 DNA is also returned to
# NUC/AA. Resident chromosome material is derived from current architecture;
# any partially or fully synthesized second copy contributes the NUC actually
# spent on replication. Thus cell death cannot silently erase modeled C/N/P.
# ATP/ADP and NAD/NADH are carrier-state variables with zero modeled C/N/P;
# stored chemical energy/redox state dissipates at lysis while their implicit
# scaffold material is represented by recycled cellular biomass.
func releasable_pools(config) -> Dictionary:
	var result: Dictionary = {}
	for metabolite_id in MetaboliteCatalogScript.extracellular_ids():
		var field_name: String = MetaboliteCatalogScript.extracellular_field(metabolite_id)
		result[field_name] = maxf(0.0, pool(metabolite_id))

	var structural_biomass: float = maxf(0.0, pool("BIO"))
	var aa_field: String = MetaboliteCatalogScript.extracellular_field("AA")
	var lip_field: String = MetaboliteCatalogScript.extracellular_field("LIP")
	var nuc_field: String = MetaboliteCatalogScript.extracellular_field("NUC")
	var dna_nuc_material: float = 0.0
	if genome != null:
		dna_nuc_material = DNAReplicationScript.total_cell_dna_nuc_material(self, config)
	result[aa_field] = (
		float(result.get(aa_field, 0.0))
		+ 2.0 * structural_biomass
		+ total_protein() * float(config.translation_aa_cost_per_event)
	)
	result[lip_field] = float(result.get(lip_field, 0.0)) + structural_biomass
	result[nuc_field] = (
		float(result.get(nuc_field, 0.0))
		+ 2.0 * structural_biomass
		+ total_mrna() * float(config.transcription_nuc_cost_per_event)
		+ dna_nuc_material
	)
	return result

func total_adenylate() -> float:
	return pool("ATP") + pool("ADP")

func total_redox_currency() -> float:
	return pool("NAD") + pool("NADH")

func _assert_state(config) -> void:
	assert(volume > 0.0)
	assert(damage >= -1e-10)
	assert(energy_debt >= -1e-10)
	assert(replication_progress >= -1e-10 and replication_progress <= 1.0 + 1e-10)
	assert(replication_gene_equivalents_copied >= -1e-10)
	assert(replication_atp_spent >= -1e-10 and replication_nuc_spent >= -1e-10)
	assert(replication_repair_activity_integral >= -1e-10)
	for amount_variant in signalling_state.values():
		assert(float(amount_variant) >= -1e-10)
	assert(motor_heading.length_squared() > 1e-18)
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
		+ replication_progress * 43.0
		+ replication_gene_equivalents_copied * 0.047
		+ replication_atp_spent * 0.053
		+ replication_nuc_spent * 0.059
		+ replication_repair_activity_integral * 0.061
		+ motor_heading.x * 0.067
		+ motor_heading.y * 0.071
		+ SignallingSystemScript.checksum(signalling_state) * 0.073
		+ MetabolicSolverScript.checksum(metabolites) * 0.001
		+ ExpressionSystemScript.checksum(expression_state) * 0.017
	)
	if genome != null:
		result += float(genome.checksum()) * 0.013
	return result
