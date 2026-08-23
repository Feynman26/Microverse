extends RefCounted
class_name M5CAmplifiedCircuit

const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")

const CONDITION_STABLE: String = "stable"
const CONDITION_FLUCTUATING: String = "fluctuating"
const AMPLIFIED_PROMOTER: int = 10000

# This validation construct changes only the shared basal abundance of R03 and
# its O2-compatible regulator. Constitutive and responsive genomes still differ
# in exactly one inherited field: the R03 regulatory motif. No fitness term,
# genotype-specific cost, or named response code is introduced.
static func create_genomes() -> Dictionary:
	var common: Array = [
		GeneScript.new(1, 6200, 0x1357, 101, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(2, 5400, 0x2468, 102, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(3, AMPLIFIED_PROMOTER, BaseExperiment.OXPHOS_PROTEIN_SIGNATURE, 103, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(4, 4800, 0x48AD, 104, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(5, 3900, 0x5ACE, 105, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(6, 6600, 0x6BDF, 106, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(7, 5700, 0x7CE1, 107, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(8, 4500, 0x8DF2, 108, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(9, 6000, 0xAF14, 109, BaseExperiment.DORMANT_MOTIF),
		GeneScript.new(10, AMPLIFIED_PROMOTER, BaseExperiment.OXYGEN_LIGAND_SIGNATURE, 110, BaseExperiment.DORMANT_MOTIF)
	]
	var constitutive = GenomeScript.new(common)
	var responsive = constitutive.deep_copy()
	responsive.get_gene_by_locus(3).regulatory_signature = BaseExperiment.OXYGEN_LIGAND_SIGNATURE
	responsive.validate()
	return {"constitutive": constitutive, "responsive": responsive}

static func run_lineage(
	seed: int,
	genotype_name: String,
	condition: String,
	phase_offset_ticks: int = 0,
	max_ticks: int = BaseExperiment.DEFAULT_MAX_TICKS
) -> Dictionary:
	assert(genotype_name == "constitutive" or genotype_name == "responsive")
	assert(condition == CONDITION_STABLE or condition == CONDITION_FLUCTUATING)
	var config = BaseExperiment.create_config(seed)
	config.max_cells = 1
	var sim = SimulationEngineScript.new(config)
	var genome = create_genomes()[genotype_name]
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
			var phase: int = ((tick + phase_offset_ticks) / BaseExperiment.PHASE_TICKS) % 2
			oxygen_value = BaseExperiment.FLUCTUATING_OXYGEN_HIGH if phase == 0 else BaseExperiment.FLUCTUATING_OXYGEN_LOW
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
		"ticks": realized_ticks,
		"time_min": elapsed_min,
		"reached_division": reached_division,
		"growth_rate": growth_rate
	}

static func paired_differential(seed: int) -> Dictionary:
	var stable_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_STABLE)
	var stable_r: Dictionary = run_lineage(seed, "responsive", CONDITION_STABLE)
	var stable_adv: float = float(stable_r["growth_rate"]) - float(stable_c["growth_rate"])
	var reference_rate: float = maxf(1e-12, 0.5 * (float(stable_r["growth_rate"]) + float(stable_c["growth_rate"])))
	var offsets: Array = [0, BaseExperiment.PHASE_TICKS / 2, BaseExperiment.PHASE_TICKS, (3 * BaseExperiment.PHASE_TICKS) / 2]
	var advantages: Array = []
	var all_divided: bool = bool(stable_c["reached_division"]) and bool(stable_r["reached_division"])
	for offset_variant in offsets:
		var offset: int = int(offset_variant)
		var fluct_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_FLUCTUATING, offset)
		var fluct_r: Dictionary = run_lineage(seed, "responsive", CONDITION_FLUCTUATING, offset)
		all_divided = all_divided and bool(fluct_c["reached_division"]) and bool(fluct_r["reached_division"])
		advantages.append(float(fluct_r["growth_rate"]) - float(fluct_c["growth_rate"]))
	var fluct_adv: float = _mean(advantages)
	return {
		"seed": seed,
		"all_divided": all_divided,
		"stable_advantage": stable_adv,
		"fluctuating_advantage": fluct_adv,
		"normalized_differential": (fluct_adv - stable_adv) / reference_rate
	}

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
