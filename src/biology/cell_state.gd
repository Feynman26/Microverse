extends RefCounted
class_name CellState

var id: int
var parent_id: int
var generation: int
var birth_tick: int
var position: Vector2

var alive: bool = true
var death_reason: String = ""

var volume: float = 1.0
var internal_glucose: float = 0.0
var internal_oxygen: float = 0.0
var atp: float = 2.0
var precursor: float = 0.0
var ros: float = 0.0
var damage: float = 0.0
var energy_debt: float = 0.0

# M0-M2 ancestral traits. These become genome/proteome-derived values in M3+.
var glucose_transport_scale: float = 1.0
var oxygen_transport_scale: float = 1.0
var respiration_scale: float = 1.0
var growth_scale: float = 1.0
var repair_scale: float = 1.0

func _init(p_id: int = 0, p_parent_id: int = -1, p_generation: int = 0, p_birth_tick: int = 0, p_position: Vector2 = Vector2.ZERO, p_volume: float = 1.0) -> void:
	id = p_id
	parent_id = p_parent_id
	generation = p_generation
	birth_tick = p_birth_tick
	position = p_position
	volume = p_volume

# Transport is split into request/allocation/application phases by the engine.
# That avoids a hidden fitness advantage for whichever cell happens to be
# iterated first when several cells compete for the same finite grid resource.
func transport_requests(dt: float, world, config) -> Dictionary:
	if not alive:
		return {"glucose": 0.0, "oxygen": 0.0}

	var local_glucose: float = float(world.sample("glucose", position))
	var glucose_vmax: float = float(config.glucose_transport_vmax)
	var glucose_km: float = float(config.glucose_transport_km)
	var pool_capacity_per_volume: float = float(config.intracellular_pool_capacity_per_volume)
	var glucose_rate: float = glucose_vmax * glucose_transport_scale * local_glucose / (glucose_km + local_glucose)
	var glucose_capacity: float = maxf(0.0, pool_capacity_per_volume * volume - internal_glucose)
	var requested_glucose: float = minf(glucose_rate * dt, glucose_capacity)

	var local_oxygen: float = float(world.sample("oxygen", position))
	var oxygen_vmax: float = float(config.oxygen_transport_vmax)
	var oxygen_km: float = float(config.oxygen_transport_km)
	var oxygen_rate: float = oxygen_vmax * oxygen_transport_scale * local_oxygen / (oxygen_km + local_oxygen)
	var oxygen_capacity: float = maxf(0.0, pool_capacity_per_volume * volume - internal_oxygen)
	var requested_oxygen: float = minf(oxygen_rate * dt, oxygen_capacity)

	return {"glucose": requested_glucose, "oxygen": requested_oxygen}

func apply_uptake(glucose_amount: float, oxygen_amount: float) -> void:
	assert(glucose_amount >= 0.0 and oxygen_amount >= 0.0)
	internal_glucose += glucose_amount
	internal_oxygen += oxygen_amount

func step_intracellular(dt: float, config) -> void:
	if not alive:
		return
	_metabolize(dt, config)
	_pay_maintenance(dt, config)
	_update_damage_and_repair(dt, config)
	_grow(dt, config)
	_check_viability(config)
	_assert_state()

func _metabolize(dt: float, config) -> void:
	var catalytic_limit: float = float(config.respiration_vmax) * respiration_scale * volume * dt
	var oxygen_per_glucose: float = float(config.oxygen_per_glucose)
	var substrate_limit: float = minf(internal_glucose, internal_oxygen / oxygen_per_glucose)
	var flux: float = minf(catalytic_limit, substrate_limit)
	if flux <= 0.0:
		return
	internal_glucose -= flux
	internal_oxygen -= flux * oxygen_per_glucose
	atp += flux * float(config.atp_yield_per_glucose)
	precursor += flux * float(config.precursor_yield_per_glucose)
	ros += flux * float(config.ros_yield_per_glucose)

func _pay_maintenance(dt: float, config) -> void:
	var required: float = float(config.maintenance_atp_rate_per_volume) * volume * dt
	var paid: float = minf(atp, required)
	atp -= paid
	var unmet: float = required - paid
	if unmet > 0.0:
		energy_debt += unmet
	else:
		energy_debt = maxf(0.0, energy_debt - required * 0.25)

func _update_damage_and_repair(dt: float, config) -> void:
	var ros_decay_rate: float = float(config.basal_ros_decay_rate)
	var ros_damage_rate: float = float(config.ros_damage_rate)
	ros = maxf(0.0, ros - ros_decay_rate * ros * dt)
	damage += ros_damage_rate * ros * dt

	var possible_repair: float = float(config.basal_repair_rate) * repair_scale * dt
	possible_repair = minf(possible_repair, damage)
	var repair_cost: float = possible_repair * float(config.repair_atp_cost)
	if repair_cost > atp and repair_cost > 0.0:
		possible_repair *= atp / repair_cost
		repair_cost = atp
	atp -= repair_cost
	damage = maxf(0.0, damage - possible_repair)

func _grow(dt: float, config) -> void:
	var kinetic_limit: float = float(config.growth_vmax) * growth_scale * volume * dt
	var growth_atp_per_precursor: float = float(config.growth_atp_per_precursor)
	var energetic_limit: float = atp / growth_atp_per_precursor
	var flux: float = minf(precursor, minf(kinetic_limit, energetic_limit))
	if flux <= 0.0:
		return
	precursor -= flux
	atp -= flux * growth_atp_per_precursor
	volume += flux * float(config.volume_yield_per_precursor)

func _check_viability(config) -> void:
	if damage >= float(config.lethal_damage):
		alive = false
		death_reason = "damage"
	elif energy_debt >= float(config.lethal_energy_debt):
		alive = false
		death_reason = "energy_failure"

func ready_to_divide(config) -> bool:
	return alive and volume >= float(config.division_volume) and atp >= float(config.division_atp_cost)

func create_daughters(first_id: int, second_id: int, tick: int, rng, world, config) -> Array:
	assert(ready_to_divide(config))
	atp -= float(config.division_atp_cost)
	var partition_jitter: float = float(config.partition_jitter)
	var ratio: float = 0.5 + float(rng.randf_range(-partition_jitter, partition_jitter))
	var offset_scale: float = float(config.daughter_offset_grid)
	var first_offset: Vector2 = Vector2(float(rng.randf_range(-1.0, 1.0)), float(rng.randf_range(-1.0, 1.0))).normalized() * offset_scale
	if first_offset == Vector2.ZERO:
		first_offset = Vector2(offset_scale, 0.0)
	var second_offset: Vector2 = -first_offset

	var first = CellState.new(first_id, id, generation + 1, tick, world.clamp_position(position + first_offset), volume * ratio)
	var second = CellState.new(second_id, id, generation + 1, tick, world.clamp_position(position + second_offset), volume * (1.0 - ratio))
	_copy_partitioned_state(first, ratio)
	_copy_partitioned_state(second, 1.0 - ratio)
	_copy_traits(first)
	_copy_traits(second)
	alive = false
	death_reason = "division"
	return [first, second]

func _copy_partitioned_state(child, ratio: float) -> void:
	child.internal_glucose = internal_glucose * ratio
	child.internal_oxygen = internal_oxygen * ratio
	child.atp = atp * ratio
	child.precursor = precursor * ratio
	child.ros = ros * ratio
	child.damage = damage * ratio
	child.energy_debt = energy_debt * ratio

func _copy_traits(child) -> void:
	child.glucose_transport_scale = glucose_transport_scale
	child.oxygen_transport_scale = oxygen_transport_scale
	child.respiration_scale = respiration_scale
	child.growth_scale = growth_scale
	child.repair_scale = repair_scale

func releasable_pools() -> Dictionary:
	return {
		"glucose": maxf(0.0, internal_glucose),
		"oxygen": maxf(0.0, internal_oxygen)
	}

func _assert_state() -> void:
	assert(volume > 0.0)
	assert(internal_glucose >= -1e-10)
	assert(internal_oxygen >= -1e-10)
	assert(atp >= -1e-10)
	assert(precursor >= -1e-10)
	assert(ros >= -1e-10)
	assert(damage >= -1e-10)
	assert(energy_debt >= -1e-10)

func checksum() -> float:
	return (
		float(id) * 0.001
		+ position.x * 3.0
		+ position.y * 5.0
		+ volume * 7.0
		+ internal_glucose * 11.0
		+ internal_oxygen * 13.0
		+ atp * 17.0
		+ precursor * 19.0
		+ ros * 23.0
		+ damage * 29.0
		+ energy_debt * 31.0
	)
