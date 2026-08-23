extends RefCounted
class_name M5CTimescaleSpike

const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")

const CONDITION_STABLE: String = "stable"
const CONDITION_FLUCTUATING: String = "fluctuating"

# Research-only schedule sweep. This file does not change production simulation
# semantics. It asks whether the null result at 20-minute O2 phases is specific
# to the relationship between environmental timescale and the existing 20-minute
# protein turnover time constant.

static func run_sweep(
	seeds: Array,
	phase_ticks_values: Array,
	max_ticks: int = BaseExperiment.DEFAULT_MAX_TICKS
) -> Dictionary:
	assert(not seeds.is_empty() and not phase_ticks_values.is_empty())
	var by_phase: Dictionary = {}
	for phase_variant in phase_ticks_values:
		var phase_ticks: int = int(phase_variant)
		assert(phase_ticks > 0 and phase_ticks % 2 == 0)
		by_phase[phase_ticks] = {
			"phase_ticks": phase_ticks,
			"phase_minutes": float(phase_ticks) * 0.10,
			"normalized_differentials": [],
			"seed_results": []
		}

	var all_divided: bool = true
	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var stable_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_STABLE, BaseExperiment.PHASE_TICKS, 0, max_ticks)
		var stable_r: Dictionary = run_lineage(seed, "responsive", CONDITION_STABLE, BaseExperiment.PHASE_TICKS, 0, max_ticks)
		all_divided = all_divided and bool(stable_c["reached_division"]) and bool(stable_r["reached_division"])
		var stable_advantage: float = float(stable_r["growth_rate"]) - float(stable_c["growth_rate"])
		var reference_rate: float = maxf(
			1e-12,
			0.5 * (float(stable_r["growth_rate"]) + float(stable_c["growth_rate"]))
		)

		for phase_variant in phase_ticks_values:
			var phase_ticks: int = int(phase_variant)
			var offsets: Array = [0, phase_ticks / 2, phase_ticks, (3 * phase_ticks) / 2]
			var fluctuating_advantages: Array = []
			var runs: Array = []
			for offset_variant in offsets:
				var offset: int = int(offset_variant)
				var fluct_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_FLUCTUATING, phase_ticks, offset, max_ticks)
				var fluct_r: Dictionary = run_lineage(seed, "responsive", CONDITION_FLUCTUATING, phase_ticks, offset, max_ticks)
				all_divided = all_divided and bool(fluct_c["reached_division"]) and bool(fluct_r["reached_division"])
				fluctuating_advantages.append(float(fluct_r["growth_rate"]) - float(fluct_c["growth_rate"]))
				runs.append({"offset": offset, "constitutive": fluct_c, "responsive": fluct_r})
			var fluctuating_advantage: float = _mean(fluctuating_advantages)
			var normalized: float = (fluctuating_advantage - stable_advantage) / reference_rate
			var record: Dictionary = by_phase[phase_ticks]
			record["normalized_differentials"].append(normalized)
			record["seed_results"].append({
				"seed": seed,
				"stable_advantage": stable_advantage,
				"fluctuating_advantage": fluctuating_advantage,
				"normalized_differential": normalized,
				"runs": runs
			})
			by_phase[phase_ticks] = record

	for phase_variant in phase_ticks_values:
		var phase_ticks: int = int(phase_variant)
		var record: Dictionary = by_phase[phase_ticks]
		record["mean_normalized_differential"] = _mean(record["normalized_differentials"])
		record["sample_sd"] = _sample_sd(record["normalized_differentials"])
		by_phase[phase_ticks] = record

	return {"all_divided": all_divided, "by_phase": by_phase}

static func run_lineage(
	seed: int,
	genotype_name: String,
	condition: String,
	phase_ticks: int,
	phase_offset_ticks: int = 0,
	max_ticks: int = BaseExperiment.DEFAULT_MAX_TICKS
) -> Dictionary:
	assert(genotype_name == "constitutive" or genotype_name == "responsive")
	assert(condition == CONDITION_STABLE or condition == CONDITION_FLUCTUATING)
	assert(phase_ticks > 0 and phase_offset_ticks >= 0 and max_ticks > 0)
	var config = BaseExperiment.create_config(seed)
	config.max_cells = 1
	var sim = SimulationEngineScript.new(config)
	var genome = BaseExperiment.create_competitor_genomes()[genotype_name]
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2(1.5, 1.5), config.ancestor_volume)
	cell.genome = genome.deep_copy()
	cell.initialize_molecular_state(config)
	sim.cells = [cell]
	sim.next_cell_id = 2

	var realized_ticks: int = 0
	var reached_division: bool = false
	for tick in range(max_ticks):
		var oxygen_value: float = BaseExperiment.STABLE_OXYGEN
		if condition == CONDITION_FLUCTUATING:
			oxygen_value = oxygen_for_tick(tick + phase_offset_ticks, phase_ticks)
		_maintain_environment(sim, oxygen_value)
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
	return {
		"seed": seed,
		"genotype": genotype_name,
		"condition": condition,
		"phase_ticks": phase_ticks,
		"phase_offset_ticks": phase_offset_ticks,
		"ticks": realized_ticks,
		"time_min": elapsed_min,
		"reached_division": reached_division,
		"growth_rate": growth_rate
	}

static func oxygen_for_tick(tick: int, phase_ticks: int) -> float:
	assert(tick >= 0 and phase_ticks > 0)
	var phase: int = (tick / phase_ticks) % 2
	return BaseExperiment.FLUCTUATING_OXYGEN_HIGH if phase == 0 else BaseExperiment.FLUCTUATING_OXYGEN_LOW

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
	if values.is_empty():
		return 0.0
	var total: float = 0.0
	for value_variant in values:
		total += float(value_variant)
	return total / float(values.size())

static func _sample_sd(values: Array) -> float:
	if values.size() < 2:
		return 0.0
	var mean_value: float = _mean(values)
	var total: float = 0.0
	for value_variant in values:
		var delta: float = float(value_variant) - mean_value
		total += delta * delta
	return sqrt(total / float(values.size() - 1))
