extends RefCounted
class_name SignallingSystem

# Authoritative reversible molecular activation. Keys are ordinary 16-bit
# molecular signatures inherited from the receptor cohort; there is no memory,
# signal, ligand-purpose or behavioral label. Occupied receptor molecules drive
# activation and the activated pool relaxes/decays when occupancy disappears.
static func step(
	state: Dictionary,
	occupied_by_signature: Dictionary,
	dt: float,
	activation_rate_per_min: float,
	decay_rate_per_min: float
) -> Dictionary:
	assert(dt >= 0.0)
	assert(activation_rate_per_min >= 0.0)
	assert(decay_rate_per_min >= 0.0)
	var prior: Dictionary = state.duplicate(true)
	var signatures: Dictionary = {}
	for signature_variant in prior.keys():
		signatures[int(signature_variant)] = true
	for signature_variant in occupied_by_signature.keys():
		signatures[int(signature_variant)] = true

	var next_state: Dictionary = {}
	var activated_total: float = 0.0
	var deactivated_total: float = 0.0
	var signature_keys: Array = signatures.keys()
	signature_keys.sort()
	for signature_variant in signature_keys:
		var signature: int = int(signature_variant)
		var old_active: float = maxf(0.0, float(prior.get(signature, 0.0)))
		var occupied: float = maxf(0.0, float(occupied_by_signature.get(signature, 0.0)))
		# Bounded first-order activation toward the currently occupied receptor
		# population, followed by independent first-order deactivation. Exact
		# exponential factors make the update stable for arbitrary dt.
		var activation_fraction: float = 1.0 - exp(-activation_rate_per_min * dt)
		var driven: float = old_active + maxf(0.0, occupied - old_active) * activation_fraction
		var decay_fraction: float = 1.0 - exp(-decay_rate_per_min * dt)
		var deactivated: float = driven * decay_fraction
		var next_active: float = maxf(0.0, driven - deactivated)
		activated_total += maxf(0.0, driven - old_active)
		deactivated_total += deactivated
		if next_active > 1e-12:
			next_state[signature] = next_active

	state.clear()
	for signature_variant in next_state.keys():
		state[signature_variant] = next_state[signature_variant]
	return {
		"activated": activated_total,
		"deactivated": deactivated_total,
		"active_total": total_active(state),
		"state": state.duplicate(true)
	}

static func total_active(state: Dictionary) -> float:
	var result: float = 0.0
	for amount_variant in state.values():
		result += maxf(0.0, float(amount_variant))
	return result

static func partition(state: Dictionary, first_fraction: float) -> Array:
	assert(first_fraction >= 0.0 and first_fraction <= 1.0)
	var first: Dictionary = {}
	var second: Dictionary = {}
	var signatures: Array = state.keys()
	signatures.sort()
	for signature_variant in signatures:
		var signature: int = int(signature_variant)
		var amount: float = maxf(0.0, float(state[signature_variant]))
		first[signature] = amount * first_fraction
		second[signature] = amount * (1.0 - first_fraction)
	return [first, second]

static func checksum(state: Dictionary) -> float:
	var result: float = 0.0
	var signatures: Array = state.keys()
	signatures.sort()
	for index in range(signatures.size()):
		var signature: int = int(signatures[index])
		result += float(signature) * 0.000001 * float(index + 1)
		result += float(state[signature]) * 0.019 * float(index + 1)
	return result
