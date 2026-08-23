extends RefCounted
class_name M8RegulatoryConfirmation

const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const M5CRegulatorySelectionScript = preload("res://src/experiments/m5c_regulatory_selection.gd")

# Frozen before confirmatory execution. These seeds are distinct from the prior
# M5-C panels and must not be changed in response to their observed outcomes.
const CONFIRMATORY_SEEDS: Array[int] = [730101, 730102, 730103, 730104, 730105, 730106, 730107, 730108]
const HORIZON_TICKS: int = 1200
const PRIMARY_ENDPOINT: String = "paired responsive-minus-constitutive cell-count difference"
const SUPPORT_MEDIAN_EFFECT: float = 1.0
const SUPPORT_MIN_POSITIVE_REPLICATES: int = 6
const AGAINST_MEDIAN_EFFECT: float = -1.0
const AGAINST_MIN_NEGATIVE_REPLICATES: int = 6

static func run_confirmatory_panel() -> Dictionary:
	var rows: Array = []
	var effects: Array = []
	var positive_count: int = 0
	var negative_count: int = 0
	var zero_count: int = 0
	var stable_extinctions: int = 0
	var fluctuating_extinctions: int = 0

	for seed in CONFIRMATORY_SEEDS:
		var stable: Dictionary = run_condition(seed, M5CRegulatorySelectionScript.CONDITION_STABLE)
		var fluctuating: Dictionary = run_condition(seed, M5CRegulatorySelectionScript.CONDITION_FLUCTUATING)
		var paired_effect: float = float(fluctuating["responsive_minus_constitutive"]) - float(stable["responsive_minus_constitutive"])
		effects.append(paired_effect)
		if paired_effect > 0.0:
			positive_count += 1
		elif paired_effect < 0.0:
			negative_count += 1
		else:
			zero_count += 1
		if bool(stable["extinct"]):
			stable_extinctions += 1
		if bool(fluctuating["extinct"]):
			fluctuating_extinctions += 1
		rows.append({
			"seed": seed,
			"stable": stable,
			"fluctuating": fluctuating,
			"paired_effect": paired_effect
		})

	var sorted_effects: Array = effects.duplicate()
	sorted_effects.sort()
	var median_effect: float = _quantile_sorted(sorted_effects, 0.50)
	var interpretation: String = "inconclusive"
	if median_effect >= SUPPORT_MEDIAN_EFFECT and positive_count >= SUPPORT_MIN_POSITIVE_REPLICATES:
		interpretation = "supports_environment_dependent_responsive_advantage"
	elif median_effect <= AGAINST_MEDIAN_EFFECT and negative_count >= AGAINST_MIN_NEGATIVE_REPLICATES:
		interpretation = "evidence_against_environment_dependent_responsive_advantage"

	return {
		"panel_status": "frozen_confirmatory",
		"seeds": CONFIRMATORY_SEEDS.duplicate(),
		"horizon_ticks": HORIZON_TICKS,
		"primary_endpoint": PRIMARY_ENDPOINT,
		"support_rule": {
			"median_effect_at_least": SUPPORT_MEDIAN_EFFECT,
			"positive_replicates_at_least": SUPPORT_MIN_POSITIVE_REPLICATES
		},
		"against_rule": {
			"median_effect_at_most": AGAINST_MEDIAN_EFFECT,
			"negative_replicates_at_least": AGAINST_MIN_NEGATIVE_REPLICATES
		},
		"rows": rows,
		"effects": effects,
		"median_effect": median_effect,
		"q25_effect": _quantile_sorted(sorted_effects, 0.25),
		"q75_effect": _quantile_sorted(sorted_effects, 0.75),
		"positive_replicates": positive_count,
		"negative_replicates": negative_count,
		"zero_replicates": zero_count,
		"stable_extinctions": stable_extinctions,
		"fluctuating_extinctions": fluctuating_extinctions,
		"interpretation": interpretation
	}

static func run_condition(seed: int, condition: String) -> Dictionary:
	assert(condition == M5CRegulatorySelectionScript.CONDITION_STABLE or condition == M5CRegulatorySelectionScript.CONDITION_FLUCTUATING)
	var config = M5CRegulatorySelectionScript.create_config(seed)
	var sim = SimulationEngineScript.new(config)
	var genomes: Dictionary = M5CRegulatorySelectionScript.create_competitor_genomes()
	var responsive = genomes["responsive"]
	var constitutive = genomes["constitutive"]
	var responsive_key: String = responsive.canonical_key()
	var constitutive_key: String = constitutive.canonical_key()
	M5CRegulatorySelectionScript._seed_founders(sim, responsive, constitutive, seed)

	var realized_ticks: int = 0
	for tick in range(HORIZON_TICKS):
		M5CRegulatorySelectionScript._maintain_environment(
			sim,
			M5CRegulatorySelectionScript.oxygen_for_tick(condition, tick)
		)
		sim.step(1)
		realized_ticks += 1
		if sim.cells.is_empty():
			# Extinction is a terminal fixed-horizon state. Stop computing but retain
			# the zero-cell endpoint rather than dropping the replicate.
			break

	var responsive_count: int = 0
	var constitutive_count: int = 0
	for cell in sim.cells:
		var key: String = cell.genome.canonical_key()
		if key == responsive_key:
			responsive_count += 1
		elif key == constitutive_key:
			constitutive_count += 1
		else:
			assert(false, "Mutation-free M8 confirmation produced an unexpected genotype")
	return {
		"seed": seed,
		"condition": condition,
		"planned_horizon_ticks": HORIZON_TICKS,
		"realized_ticks": realized_ticks,
		"population": sim.population_size(),
		"responsive": responsive_count,
		"constitutive": constitutive_count,
		"responsive_minus_constitutive": responsive_count - constitutive_count,
		"division_events": M5CRegulatorySelectionScript._division_event_count(sim),
		"max_generation": sim.maximum_generation(),
		"extinct": sim.cells.is_empty(),
		"final_checksum": sim.checksum()
	}

static func _quantile_sorted(sorted: Array, fraction: float) -> float:
	assert(not sorted.is_empty() and fraction >= 0.0 and fraction <= 1.0)
	if sorted.size() == 1:
		return float(sorted[0])
	var position: float = fraction * float(sorted.size() - 1)
	var low: int = int(floor(position))
	var high: int = int(ceil(position))
	if low == high:
		return float(sorted[low])
	var weight: float = position - float(low)
	return lerpf(float(sorted[low]), float(sorted[high]), weight)
