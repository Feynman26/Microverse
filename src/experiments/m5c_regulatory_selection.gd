extends RefCounted
class_name M5CRegulatorySelection

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

# M5-C is a deliberately narrow evolutionary-selection assay. It does not add
# an environment scheduler to production biology and does not assign fitness.
# Cells grow/divide/die through ordinary Microverse physiology; the harness only
# imposes controlled extracellular resource concentrations and measures outcomes.

const CONDITION_STABLE: String = "stable"
const CONDITION_FLUCTUATING: String = "fluctuating"

const DORMANT_MOTIF: int = 0xF139
const OXYGEN_LIGAND_SIGNATURE: int = 0xCCCC
const OXPHOS_PROTEIN_SIGNATURE: int = 0x369C

# With oxygen Km = 0.60, these values match the time-average empty-cell oxygen
# transport propensity exactly for a 50/50 square wave:
#   0.5/(0.6+0.5) = 0.5 * 6/(0.6+6) = 0.454545...
const STABLE_OXYGEN: float = 0.5
const FLUCTUATING_OXYGEN_HIGH: float = 6.0
const FLUCTUATING_OXYGEN_LOW: float = 0.0
const PHASE_TICKS: int = 200

const GLUCOSE_RESERVOIR: float = 4.0
const NITROGEN_RESERVOIR: float = 3.0
const PHOSPHORUS_RESERVOIR: float = 2.0
const FOUNDERS_PER_GENOTYPE: int = 2
const DEFAULT_TARGET_DIVISIONS: int = 20
const DEFAULT_MAX_TICKS: int = 3600

static func create_config(seed: int):
	var config = SimConfigScript.new()
	config.seed = seed
	config.world_width = 4
	config.world_height = 4
	config.max_cells = 128
	config.glucose_diffusion = 0.0
	config.oxygen_diffusion = 0.0
	config.nitrogen_diffusion = 0.0
	config.phosphorus_diffusion = 0.0
	config.mutation_enabled = false

	# Strong enough to create a measurable regulatory trade-off while remaining
	# inside the same bounded M5-B molecular grammar.
	config.regulatory_gain = 1.5
	config.regulatory_min_factor = 0.10
	config.regulatory_max_factor = 1.75
	config.allosteric_gain = 1.0
	config.allosteric_min_factor = 0.05
	config.allosteric_max_factor = 2.0
	config.validate()
	return config

static func create_competitor_genomes() -> Dictionary:
	# Minimal viable compressed network plus one catalytically inactive regulator.
	# The competitors are genetically identical except for locus 3's promoter
	# motif: responsive binds the O2-compatible regulator; constitutive does not.
	var common: Array = [
		GeneScript.new(1, 6200, 0x1357, 101, DORMANT_MOTIF),
		GeneScript.new(2, 5400, 0x2468, 102, DORMANT_MOTIF),
		GeneScript.new(3, 7100, OXPHOS_PROTEIN_SIGNATURE, 103, DORMANT_MOTIF),
		GeneScript.new(4, 4800, 0x48AD, 104, DORMANT_MOTIF),
		GeneScript.new(5, 3900, 0x5ACE, 105, DORMANT_MOTIF),
		GeneScript.new(6, 6600, 0x6BDF, 106, DORMANT_MOTIF),
		GeneScript.new(7, 5700, 0x7CE1, 107, DORMANT_MOTIF),
		GeneScript.new(8, 4500, 0x8DF2, 108, DORMANT_MOTIF),
		GeneScript.new(9, 6000, 0xAF14, 109, DORMANT_MOTIF),
		# 0xCCCC exactly matches O2 ligand identity. Its sequence bits make it a
		# ligand-inhibited repressor in the generic M5-B grammar and it is outside
		# every M4 catalytic radius.
		GeneScript.new(10, 5000, OXYGEN_LIGAND_SIGNATURE, 110, DORMANT_MOTIF)
	]
	var constitutive = GenomeScript.new(common)
	var responsive = constitutive.deep_copy()
	responsive.get_gene_by_locus(3).regulatory_signature = OXYGEN_LIGAND_SIGNATURE
	responsive.validate()
	return {"constitutive": constitutive, "responsive": responsive}

static func empty_cell_oxygen_transport_factor(external_oxygen: float, oxygen_km: float) -> float:
	assert(external_oxygen >= 0.0 and oxygen_km > 0.0)
	return external_oxygen / (oxygen_km + external_oxygen)

static func schedule_transport_means(config) -> Dictionary:
	var km: float = float(config.oxygen_transport_km)
	var stable: float = empty_cell_oxygen_transport_factor(STABLE_OXYGEN, km)
	var fluctuating: float = 0.5 * (
		empty_cell_oxygen_transport_factor(FLUCTUATING_OXYGEN_HIGH, km)
		+ empty_cell_oxygen_transport_factor(FLUCTUATING_OXYGEN_LOW, km)
	)
	return {"stable": stable, "fluctuating": fluctuating}

static func oxygen_for_tick(condition: String, tick: int) -> float:
	assert(condition == CONDITION_STABLE or condition == CONDITION_FLUCTUATING)
	if condition == CONDITION_STABLE:
		return STABLE_OXYGEN
	var phase: int = (tick / PHASE_TICKS) % 2
	return FLUCTUATING_OXYGEN_HIGH if phase == 0 else FLUCTUATING_OXYGEN_LOW

# Population assay retained as characterization of drift and frequency dynamics.
# It is intentionally not the final M5-C causal gate: at very small N the
# descendant-count signal can be dominated by stochastic lineage branching.
static func run_replicate(
	seed: int,
	condition: String,
	max_ticks: int = DEFAULT_MAX_TICKS,
	target_divisions: int = DEFAULT_TARGET_DIVISIONS
) -> Dictionary:
	assert(max_ticks > 0 and target_divisions > 0)
	var config = create_config(seed)
	var sim = SimulationEngineScript.new(config)
	var genomes: Dictionary = create_competitor_genomes()
	var constitutive = genomes["constitutive"]
	var responsive = genomes["responsive"]
	var constitutive_key: String = constitutive.canonical_key()
	var responsive_key: String = responsive.canonical_key()

	_seed_founders(sim, responsive, constitutive, seed)
	var realized_ticks: int = 0
	var division_events: int = 0
	for tick in range(max_ticks):
		_maintain_environment(sim, oxygen_for_tick(condition, tick))
		sim.step(1)
		realized_ticks += 1
		division_events = _division_event_count(sim)
		if division_events >= target_divisions or sim.cells.is_empty():
			break

	var responsive_count: int = 0
	var constitutive_count: int = 0
	var responsive_r03: float = 0.0
	var constitutive_r03: float = 0.0
	for cell in sim.cells:
		var key: String = cell.genome.canonical_key()
		var r03_gene = cell.genome.get_gene_by_locus(3)
		var r03_amount: float = ExpressionSystemScript.current_gene_protein(cell.expression_state, r03_gene)
		if key == responsive_key:
			responsive_count += 1
			responsive_r03 += r03_amount
		elif key == constitutive_key:
			constitutive_count += 1
			constitutive_r03 += r03_amount

	var classified: int = responsive_count + constitutive_count
	assert(classified == sim.population_size(), "M5-C mutation-free competition produced an unexpected genotype")
	assert(classified < int(config.max_cells), "M5-C endpoint touched population cap and is not a valid selection assay")
	var responsive_fraction: float = 0.0 if classified == 0 else float(responsive_count) / float(classified)
	var log_ratio: float = log((float(responsive_count) + 0.5) / (float(constitutive_count) + 0.5))
	return {
		"seed": seed,
		"condition": condition,
		"ticks": realized_ticks,
		"division_events": division_events,
		"target_divisions": target_divisions,
		"reached_target": division_events >= target_divisions,
		"population": classified,
		"responsive": responsive_count,
		"constitutive": constitutive_count,
		"responsive_fraction": responsive_fraction,
		"log_ratio": log_ratio,
		"max_generation": sim.maximum_generation(),
		"mean_r03_responsive": 0.0 if responsive_count == 0 else responsive_r03 / float(responsive_count),
		"mean_r03_constitutive": 0.0 if constitutive_count == 0 else constitutive_r03 / float(constitutive_count)
	}

static func run_paired_replicates(
	seeds: Array,
	max_ticks: int = DEFAULT_MAX_TICKS,
	target_divisions: int = DEFAULT_TARGET_DIVISIONS
) -> Dictionary:
	var stable_results: Array = []
	var fluctuating_results: Array = []
	var paired_differences: Array = []
	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var stable: Dictionary = run_replicate(seed, CONDITION_STABLE, max_ticks, target_divisions)
		var fluctuating: Dictionary = run_replicate(seed, CONDITION_FLUCTUATING, max_ticks, target_divisions)
		stable_results.append(stable)
		fluctuating_results.append(fluctuating)
		paired_differences.append(float(fluctuating["log_ratio"]) - float(stable["log_ratio"]))
	return {
		"stable": stable_results,
		"fluctuating": fluctuating_results,
		"paired_differences": paired_differences,
		"mean_stable_log_ratio": _mean_metric(stable_results, "log_ratio"),
		"mean_fluctuating_log_ratio": _mean_metric(fluctuating_results, "log_ratio"),
		"mean_paired_difference": _mean_values(paired_differences)
	}

# Causal M5-C gate: measure reproductive rate of one genotype at a time, using
# the production SimulationEngine/CellState/ExpressionSystem/MetabolicSolver.
# max_cells=1 prevents SimulationEngine from replacing the mature founder with
# daughters; the harness can therefore observe the first instant at which normal
# division criteria become satisfied without changing intracellular physics.
static func run_lineage(
	seed: int,
	genotype_name: String,
	condition: String,
	phase_offset_ticks: int = 0,
	max_ticks: int = DEFAULT_MAX_TICKS
) -> Dictionary:
	assert(genotype_name == "constitutive" or genotype_name == "responsive")
	assert(condition == CONDITION_STABLE or condition == CONDITION_FLUCTUATING)
	assert(phase_offset_ticks >= 0 and max_ticks > 0)
	var config = create_config(seed)
	config.max_cells = 1
	var sim = SimulationEngineScript.new(config)
	var genome = create_competitor_genomes()[genotype_name]
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2(1.5, 1.5), config.ancestor_volume)
	cell.genome = genome.deep_copy()
	cell.initialize_molecular_state(config)
	sim.cells = [cell]
	sim.next_cell_id = 2

	var realized_ticks: int = 0
	var reached_division: bool = false
	for tick in range(max_ticks):
		_maintain_environment(sim, oxygen_for_tick(condition, tick + phase_offset_ticks))
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
		"alive": not sim.cells.is_empty() and sim.cells[0].alive,
		"growth_rate": growth_rate,
		"volume": 0.0 if sim.cells.is_empty() else float(sim.cells[0].volume),
		"atp": 0.0 if sim.cells.is_empty() else float(sim.cells[0].pool("ATP"))
	}

# For each seed, stable growth advantage is measured once. Fluctuating growth
# advantage is averaged across four equally spaced phase offsets so the result
# cannot depend on always starting the founder at the same square-wave phase.
# The output is a dimensionless environment-selection differential normalized by
# the mean stable division rate for that seed.
static func run_lineage_selection_panel(
	seeds: Array,
	phase_offsets: Array = [0, 100, 200, 300],
	max_ticks: int = DEFAULT_MAX_TICKS
) -> Dictionary:
	assert(not seeds.is_empty() and not phase_offsets.is_empty())
	var seed_results: Array = []
	var normalized_differentials: Array = []
	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var stable_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_STABLE, 0, max_ticks)
		var stable_r: Dictionary = run_lineage(seed, "responsive", CONDITION_STABLE, 0, max_ticks)
		var stable_advantage: float = float(stable_r["growth_rate"]) - float(stable_c["growth_rate"])
		var fluctuating_advantages: Array = []
		var fluctuating_runs: Array = []
		for offset_variant in phase_offsets:
			var offset: int = int(offset_variant)
			var fluct_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_FLUCTUATING, offset, max_ticks)
			var fluct_r: Dictionary = run_lineage(seed, "responsive", CONDITION_FLUCTUATING, offset, max_ticks)
			fluctuating_runs.append({"offset": offset, "constitutive": fluct_c, "responsive": fluct_r})
			fluctuating_advantages.append(float(fluct_r["growth_rate"]) - float(fluct_c["growth_rate"]))
		var fluctuating_advantage: float = _mean_values(fluctuating_advantages)
		var differential: float = fluctuating_advantage - stable_advantage
		var reference_rate: float = maxf(
			1e-12,
			0.5 * (float(stable_r["growth_rate"]) + float(stable_c["growth_rate"]))
		)
		var normalized: float = differential / reference_rate
		normalized_differentials.append(normalized)
		seed_results.append({
			"seed": seed,
			"stable_constitutive": stable_c,
			"stable_responsive": stable_r,
			"stable_advantage": stable_advantage,
			"fluctuating_runs": fluctuating_runs,
			"fluctuating_advantage": fluctuating_advantage,
			"differential": differential,
			"normalized_differential": normalized
		})
	return {
		"seeds": seed_results,
		"normalized_differentials": normalized_differentials,
		"mean_normalized_differential": _mean_values(normalized_differentials)
	}

static func _seed_founders(sim, responsive, constitutive, seed: int) -> void:
	var order: Array = []
	for _i in range(FOUNDERS_PER_GENOTYPE):
		order.append(responsive)
		order.append(constitutive)
	if seed % 2 != 0:
		order.reverse()

	var positions: Array[Vector2] = [
		Vector2(1.25, 1.25), Vector2(2.75, 1.25),
		Vector2(1.25, 2.75), Vector2(2.75, 2.75)
	]
	var next_id: int = 1
	for i in range(order.size()):
		var cell = CellStateScript.new(next_id, -1, 0, 0, positions[i % positions.size()], sim.config.ancestor_volume)
		cell.genome = order[i].deep_copy()
		cell.initialize_molecular_state(sim.config)
		sim.cells.append(cell)
		next_id += 1
	sim.next_cell_id = next_id

static func _division_event_count(sim) -> int:
	var result: int = 0
	for event in sim.event_log:
		if String(event.get("kind", "")) == "division":
			result += 1
	return result

static func _maintain_environment(sim, oxygen_value: float) -> void:
	_fill_field(sim.world.get_field("glucose"), GLUCOSE_RESERVOIR)
	_fill_field(sim.world.get_field("oxygen"), oxygen_value)
	_fill_field(sim.world.get_field("nitrogen"), NITROGEN_RESERVOIR)
	_fill_field(sim.world.get_field("phosphorus"), PHOSPHORUS_RESERVOIR)

static func _fill_field(field, value: float) -> void:
	assert(value >= 0.0)
	for y in range(field.height):
		for x in range(field.width):
			field.set_value(x, y, value)

static func _mean_metric(results: Array, key: String) -> float:
	if results.is_empty():
		return 0.0
	var total: float = 0.0
	for result in results:
		total += float(result[key])
	return total / float(results.size())

static func _mean_values(values: Array) -> float:
	if values.is_empty():
		return 0.0
	var total: float = 0.0
	for value in values:
		total += float(value)
	return total / float(values.size())
