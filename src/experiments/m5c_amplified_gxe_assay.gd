extends RefCounted
class_name M5CAmplifiedGxeAssay

const BaseGxe = preload("res://src/experiments/m5c_gxe_assay.gd")
const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")
const AmplifiedCircuit = preload("res://src/experiments/m5c_amplified_circuit.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")

# Stronger validation construct for the same generic M5-B grammar. Both
# genotypes have identical high basal R03 and regulator promoters and therefore
# identical expression capacity/costs in the absence of regulation. They differ
# only in the inherited R03 promoter motif. This makes the regulated allocation
# a material fraction of the finite proteome without adding a fitness term.

static func run_lineage(
	seed: int,
	genotype_name: String,
	condition: String,
	phase_offset_ticks: int = 0,
	max_ticks: int = BaseGxe.MAX_TICKS
) -> Dictionary:
	assert(genotype_name == "constitutive" or genotype_name == "responsive")
	assert(condition == BaseGxe.CONDITION_HIGH_STABLE or condition == BaseGxe.CONDITION_HIGH_ANOXIC)
	var config = BaseExperiment.create_config(seed)
	config.max_cells = 1
	var sim = SimulationEngineScript.new(config)
	var genome = AmplifiedCircuit.create_genomes()[genotype_name]
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2(1.5, 1.5), config.ancestor_volume)
	cell.genome = genome.deep_copy()
	cell.initialize_molecular_state(config)
	sim.cells = [cell]
	sim.next_cell_id = 2

	var realized_ticks: int = 0
	var reached_division: bool = false
	for tick in range(max_ticks):
		_maintain_environment(sim, BaseGxe.oxygen_for_tick(condition, tick, phase_offset_ticks))
		sim.step(1)
		realized_ticks += 1
		if sim.cells.is_empty():
			break
		cell = sim.cells[0]
		if cell.ready_to_divide(config):
			reached_division = true
			break

	var elapsed_min: float = float(realized_ticks) * float(config.tick_dt_min)
	var growth_rate: float = 0.0
	if reached_division and elapsed_min > 0.0:
		growth_rate = log(2.0) / elapsed_min
	var has_cell: bool = not sim.cells.is_empty()
	if has_cell:
		cell = sim.cells[0]
	return {
		"seed": seed,
		"genotype": genotype_name,
		"condition": condition,
		"phase_offset_ticks": phase_offset_ticks,
		"ticks": realized_ticks,
		"time_min": elapsed_min,
		"reached_division": reached_division,
		"growth_rate": growth_rate,
		"alive": has_cell and bool(cell.alive),
		"death_reason": "" if not has_cell else String(cell.death_reason),
		"volume": 0.0 if not has_cell else float(cell.volume),
		"atp": 0.0 if not has_cell else float(cell.pool("ATP")),
		"damage": 0.0 if not has_cell else float(cell.damage),
		"energy_debt": 0.0 if not has_cell else float(cell.energy_debt)
	}

static func run_panel(seeds: Array, max_ticks: int = BaseGxe.MAX_TICKS) -> Dictionary:
	assert(not seeds.is_empty())
	var offsets: Array = [0, BaseGxe.PHASE_TICKS / 2, BaseGxe.PHASE_TICKS, (3 * BaseGxe.PHASE_TICKS) / 2]
	var values: Array = []
	var seed_results: Array = []
	var all_divided: bool = true
	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var stable_c: Dictionary = run_lineage(seed, "constitutive", BaseGxe.CONDITION_HIGH_STABLE, 0, max_ticks)
		var stable_r: Dictionary = run_lineage(seed, "responsive", BaseGxe.CONDITION_HIGH_STABLE, 0, max_ticks)
		all_divided = all_divided and bool(stable_c["reached_division"]) and bool(stable_r["reached_division"])
		var stable_adv: float = float(stable_r["growth_rate"]) - float(stable_c["growth_rate"])
		var reference_rate: float = maxf(1e-12, 0.5 * (float(stable_r["growth_rate"]) + float(stable_c["growth_rate"])))
		var fluct_advantages: Array = []
		var runs: Array = []
		for offset_variant in offsets:
			var offset: int = int(offset_variant)
			var fc: Dictionary = run_lineage(seed, "constitutive", BaseGxe.CONDITION_HIGH_ANOXIC, offset, max_ticks)
			var fr: Dictionary = run_lineage(seed, "responsive", BaseGxe.CONDITION_HIGH_ANOXIC, offset, max_ticks)
			all_divided = all_divided and bool(fc["reached_division"]) and bool(fr["reached_division"])
			fluct_advantages.append(float(fr["growth_rate"]) - float(fc["growth_rate"]))
			runs.append({"offset": offset, "constitutive": fc, "responsive": fr})
		var fluct_adv: float = _mean(fluct_advantages)
		var normalized: float = (fluct_adv - stable_adv) / reference_rate
		values.append(normalized)
		seed_results.append({
			"seed": seed,
			"stable_advantage": stable_adv,
			"fluctuating_advantage": fluct_adv,
			"normalized_differential": normalized,
			"stable_constitutive": stable_c,
			"stable_responsive": stable_r,
			"runs": runs
		})
	return {"all_divided": all_divided, "seed_results": seed_results, "normalized_differentials": values, "mean_normalized_differential": _mean(values), "sample_sd": _sample_sd(values)}

static func _maintain_environment(sim, oxygen_value: float) -> void:
	_fill_field(sim.world.get_field("glucose"), BaseExperiment.GLUCOSE_RESERVOIR)
	_fill_field(sim.world.get_field("oxygen"), oxygen_value)
	_fill_field(sim.world.get_field("nitrogen"), BaseExperiment.NITROGEN_RESERVOIR)
	_fill_field(sim.world.get_field("phosphorus"), BaseExperiment.PHOSPHORUS_RESERVOIR)

static func _fill_field(field, value: float) -> void:
	for y in range(field.height):
		for x in range(field.width):
			field.set_value(x, y, value)

static func _mean(values: Array) -> float:
	if values.is_empty(): return 0.0
	var total: float = 0.0
	for value_variant in values: total += float(value_variant)
	return total / float(values.size())

static func _sample_sd(values: Array) -> float:
	if values.size() < 2: return 0.0
	var mean_value: float = _mean(values)
	var total: float = 0.0
	for value_variant in values:
		var delta: float = float(value_variant) - mean_value
		total += delta * delta
	return sqrt(total / float(values.size() - 1))
