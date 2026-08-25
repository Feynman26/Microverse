extends RefCounted
class_name MotorSystem

# Upper nibble marks a generic force-generating/membrane motor localization.
# No M10 ancestral protein has Fxxx, so motility is not gifted. The lower 12
# bits determine coupling to activated receptor/signalling signatures. A sequence
# bit supplies molecular coupling polarity, allowing the same generic mechanism
# to increase or decrease heading persistence without a chemotaxis API.
const MOTOR_LOCALIZATION_MASK: int = 0xF000
const MOTOR_LOCALIZATION_VALUE: int = 0xF000
const COUPLING_MASK: int = 0x0FFF
const POLARITY_BIT: int = 0x0800
const ACTIVE_MAX_DISTANCE: int = 4
const DISTANCE_DECAY: float = 0.70

static var _affinity_cache: Dictionary = {}

static func has_motor_localization(protein_signature: int) -> bool:
	return (protein_signature & MOTOR_LOCALIZATION_MASK) == MOTOR_LOCALIZATION_VALUE

static func hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & COUPLING_MASK
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

static func coupling_affinity(motor_signature: int, molecular_signature: int) -> float:
	if not has_motor_localization(motor_signature):
		return 0.0
	var motor_binding: int = motor_signature & COUPLING_MASK
	var molecular_binding: int = molecular_signature & COUPLING_MASK
	var key: int = (motor_binding << 12) | molecular_binding
	if _affinity_cache.has(key):
		return float(_affinity_cache[key])
	var distance: int = hamming_distance(motor_binding, molecular_binding)
	var result: float = 0.0 if distance > ACTIVE_MAX_DISTANCE else exp(-DISTANCE_DECAY * float(distance))
	_affinity_cache[key] = result
	return result

static func coupling_polarity(motor_signature: int) -> float:
	# Sequence-encoded sign; no semantic activating/inhibiting class exists.
	return 1.0 if (motor_signature & POLARITY_BIT) == 0 else -1.0

static func realized_motors(expression_state: Dictionary, reference_protein_count: float) -> Dictionary:
	assert(reference_protein_count > 0.0)
	var result: Dictionary = {}
	var loci: Array = expression_state.keys()
	loci.sort()
	for locus_variant in loci:
		var locus_id: int = int(locus_variant)
		var cohorts: Dictionary = expression_state[locus_id]["protein"]
		var signatures: Array = cohorts.keys()
		signatures.sort()
		for signature_variant in signatures:
			var signature: int = int(signature_variant)
			if not has_motor_localization(signature):
				continue
			var abundance: float = maxf(0.0, float(cohorts[signature_variant]))
			if abundance > 0.0:
				result[signature] = float(result.get(signature, 0.0)) + abundance / reference_protein_count
	return result

# Fast receptor occupancy and slower reversible molecular state are two physical
# pools. Their difference is a local temporal mismatch, not a gradient vector.
# Motor proteins couple to those molecular signatures by sequence affinity. A
# cell can therefore evolve different persistence responses without any goal or
# externally supplied direction-to-resource information.
static func control_drive(
	motors: Dictionary,
	occupied_by_signature: Dictionary,
	slow_state: Dictionary
) -> float:
	var result: float = 0.0
	var motor_signatures: Array = motors.keys()
	motor_signatures.sort()
	var molecular_signatures: Dictionary = {}
	for signature_variant in occupied_by_signature.keys():
		molecular_signatures[int(signature_variant)] = true
	for signature_variant in slow_state.keys():
		molecular_signatures[int(signature_variant)] = true
	var molecular_keys: Array = molecular_signatures.keys()
	molecular_keys.sort()
	for motor_variant in motor_signatures:
		var motor_signature: int = int(motor_variant)
		var motor_amount: float = maxf(0.0, float(motors[motor_variant]))
		var polarity: float = coupling_polarity(motor_signature)
		for molecular_variant in molecular_keys:
			var molecular_signature: int = int(molecular_variant)
			var affinity: float = coupling_affinity(motor_signature, molecular_signature)
			if affinity <= 0.0:
				continue
			var fast: float = maxf(0.0, float(occupied_by_signature.get(molecular_signature, 0.0)))
			var slow: float = maxf(0.0, float(slow_state.get(molecular_signature, 0.0)))
			result += motor_amount * polarity * affinity * (fast - slow)
	return result

static func normalize_heading(heading: Vector2) -> Vector2:
	if heading.length_squared() <= 1e-18:
		return Vector2.RIGHT
	return heading.normalized()

# Stochastic direction switching. Positive molecular drive reduces the tumble
# hazard for positive-polarity motors; negative drive increases it. Different
# sequence polarity reverses that mapping. RNG is the authoritative seeded stream.
static func update_heading(
	heading: Vector2,
	control: float,
	dt: float,
	rng,
	baseline_turn_rate_per_min: float,
	control_gain: float
) -> Dictionary:
	assert(dt >= 0.0 and baseline_turn_rate_per_min >= 0.0 and control_gain >= 0.0)
	var current: Vector2 = normalize_heading(heading)
	var hazard: float = baseline_turn_rate_per_min * exp(-control_gain * control)
	var turn_probability: float = 1.0 - exp(-hazard * dt)
	var turned: bool = false
	if turn_probability > 0.0 and float(rng.randf()) < turn_probability:
		var angle: float = float(rng.randf_range(-PI, PI))
		current = Vector2(cos(angle), sin(angle))
		turned = true
	return {
		"heading": current,
		"turned": turned,
		"turn_hazard_per_min": hazard,
		"turn_probability": turn_probability,
		"control_drive": control
	}

static func movement_request(
	motors: Dictionary,
	heading: Vector2,
	dt: float,
	speed_grid_per_min_per_activity: float
) -> Dictionary:
	assert(dt >= 0.0 and speed_grid_per_min_per_activity >= 0.0)
	var total_activity: float = 0.0
	for amount_variant in motors.values():
		total_activity += maxf(0.0, float(amount_variant))
	var distance: float = total_activity * speed_grid_per_min_per_activity * dt
	return {
		"activity": total_activity,
		"requested_distance": distance,
		"requested_displacement": normalize_heading(heading) * distance
	}

static func movement_cost(requested_distance: float, atp_cost_per_grid_distance: float) -> float:
	assert(requested_distance >= 0.0 and atp_cost_per_grid_distance >= 0.0)
	return requested_distance * atp_cost_per_grid_distance

# Energy allocation is deliberately separate from geometry. If ATP cannot fund
# the requested displacement, movement is reduced proportionally; zero ATP means
# exactly zero active displacement.
static func funded_displacement(
	requested_displacement: Vector2,
	available_atp: float,
	atp_cost_per_grid_distance: float
) -> Dictionary:
	assert(available_atp >= 0.0 and atp_cost_per_grid_distance >= 0.0)
	var requested_distance: float = requested_displacement.length()
	var requested_cost: float = movement_cost(requested_distance, atp_cost_per_grid_distance)
	var scale: float = 1.0
	if requested_cost > 0.0:
		scale = minf(1.0, available_atp / requested_cost)
	var displacement: Vector2 = requested_displacement * scale
	var spent: float = requested_cost * scale
	return {
		"displacement": displacement,
		"requested_distance": requested_distance,
		"actual_distance": displacement.length(),
		"atp_spent": spent,
		"energy_scale": scale
	}
