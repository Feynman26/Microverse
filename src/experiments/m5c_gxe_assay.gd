extends RefCounted
class_name M5CGxeAssay

const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")

# Final M5-C assay: the production biology is unchanged. The experiment asks a
# cleaner genotype-by-environment question than the earlier mean-transport-
# matched schedule. Stable and fluctuating conditions share the same oxygen-rich
# state (O2=6); only the fluctuating condition periodically removes O2. Thus the
# responsive circuit can be almost unrepressed when respiration is useful and
# can down-regulate the costly R03 proteome allocation during sustained anoxia.

const CONDITION_HIGH_STABLE: String = "high_stable"
const CONDITION_HIGH_ANOXIC: String = "high_anoxic"
const HIGH_OXYGEN: float = 6.0
const LOW_OXYGEN: float = 0.0
const PHASE_TICKS: int = 400 # 40 biological min = 2 x protein lifetime tau.
const MAX_TICKS: int = 4800

static func oxygen_for_tick(condition: String, tick: int, phase_offset_ticks: int = 0) -> float:
	assert(condition == CONDITION_HIGH_STABLE or condition == CONDITION_HIGH_ANOXIC)
	if condition == CONDITION_HIGH_STABLE:
		return HIGH_OXYGEN
	var phase: int = ((tick + phase_offset_ticks) / PHASE_TICKS) % 2
	return HIGH_OXYGEN if phase == 0 else LOW_OXYGEN

static func run_lineage(
	seed: int,
	genotype_name: String,
	condition: String,
	phase_offset_ticks: int = 0,
	max_ticks: int = MAX_TICKS
) -> Dictionary:
	assert(genotype_name == "constitutive" or genotype_name == "responsive")
	assert(condition == CONDITION_HIGH_STABLE or condition == CONDITION_HIGH_ANOXIC)
	assert(phase_offset_ticks >= 0 and max_ticks > 0)
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
		_maintain_environment(sim, oxygen_for_tick(condition, tick, phase_offset_ticks))
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
		"phase_offset_ticks": phase_offset_ticks,
		"ticks": realized_ticks,
		"time_min": elapsed_min,
		"reached_division": reached_division,
		"growth_rate": growth_rate,
		"alive": not sim.cells.is_empty() and sim.cells[0].alive,
		"r03": 0.0 if sim.cells.is_empty() else _current_r03(sim.cells[0])
	}

static func run_panel(seeds: Array, max_ticks: int = MAX_TICKS) -> Dictionary:
	assert(not seeds.is_empty())
	var phase_offsets: Array = [0, PHASE_TICKS / 2, PHASE_TICKS, (3 * PHASE_TICKS) / 2]
	var seed_results: Array = []
	var normalized_differentials: Array = []
	var all_divided: bool = true

	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var stable_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_HIGH_STABLE, 0, max_ticks)
		var stable_r: Dictionary = run_lineage(seed, "responsive", CONDITION_HIGH_STABLE, 0, max_ticks)
		all_divided = all_divided and bool(stable_c["reached_division"]) and bool(stable_r["reached_division"])
		var stable_advantage: float = float(stable_r["growth_rate"]) - float(stable_c["growth_rate"])
		var reference_rate: float = maxf(
			1e-12,
			0.5 * (float(stable_r["growth_rate"]) + float(stable_c["growth_rate"]))
		)

		var fluctuating_advantages: Array = []
		var runs: Array = []
		for offset_variant in phase_offsets:
			var offset: int = int(offset_variant)
			var fluct_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_HIGH_ANOXIC, offset, max_ticks)
			var fluct_r: Dictionary = run_lineage(seed, "responsive", CONDITION_HIGH_ANOXIC, offset, max_ticks)
			all_divided = all_divided and bool(fluct_c["reached_division"]) and bool(fluct_r["reached_division"])
			var advantage: float = float(fluct_r["growth_rate"]) - float(fluct_c["growth_rate"])
			fluctuating_advantages.append(advantage)
			runs.append({"offset": offset, "constitutive": fluct_c, "responsive": fluct_r})

		var fluctuating_advantage: float = _mean(fluctuating_advantages)
		var normalized: float = (fluctuating_advantage - stable_advantage) / reference_rate
		normalized_differentials.append(normalized)
		seed_results.append({
			"seed": seed,
			"stable_constitutive": stable_c,
			"stable_responsive": stable_r,
			"stable_advantage": stable_advantage,
			"fluctuating_advantage": fluctuating_advantage,
			"normalized_differential": normalized,
			"runs": runs
		})

	return {
		"all_divided": all_divided,
		"seed_results": seed_results,
		"normalized_differentials": normalized_differentials,
		"mean_normalized_differential": _mean(normalized_differentials),
		"sample_sd": _sample_sd(normalized_differentials)
	}

static func _current_r03(cell) -> float:
	var gene = cell.genome.get_gene_by_locus(3)
	return BaseExperiment.ExpressionSystemScript.current_gene_protein(cell.expression_state, gene)

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
