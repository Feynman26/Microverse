extends RefCounted
class_name M5CBurdenCircuit

const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")
const BaseGxe = preload("res://src/experiments/m5c_gxe_assay.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")

const CONDITION_HIGH_STABLE: String = BaseGxe.CONDITION_HIGH_STABLE
const CONDITION_HIGH_ANOXIC: String = BaseGxe.CONDITION_HIGH_ANOXIC
const OXYGEN_REPRESSOR_SIGNATURE: int = 0xCCCC
const NEUTRAL_BURDEN_SIGNATURE: int = 0xC63F
const HIGH_PROMOTER: int = 10000
const DORMANT_MOTIF: int = BaseExperiment.DORMANT_MOTIF
const BURDEN_LOCI: Array = [11, 12, 13]
const MAX_TICKS: int = 7200

# Capability assay for the M5 exit gate. Three extra loci encode the same
# deliberately non-catalytic protein. Both architectures pay identical basal
# expression costs and have identical coding sequences. Only the responsive
# architecture connects the burden promoters to the generic O2-compatible
# repressor. In high O2 allosteric inhibition releases repression; during
# anoxia repression lowers unnecessary expression and frees finite proteome,
# ribosome throughput and ATP/material. No fitness or genotype-specific cost is
# introduced.
static func create_genomes() -> Dictionary:
	var common: Array = [
		GeneScript.new(1, 6200, 0x1357, 101, DORMANT_MOTIF),
		GeneScript.new(2, 5400, 0x2468, 102, DORMANT_MOTIF),
		GeneScript.new(3, 7100, 0x369C, 103, DORMANT_MOTIF),
		GeneScript.new(4, 4800, 0x48AD, 104, DORMANT_MOTIF),
		GeneScript.new(5, 3900, 0x5ACE, 105, DORMANT_MOTIF),
		GeneScript.new(6, 6600, 0x6BDF, 106, DORMANT_MOTIF),
		GeneScript.new(7, 5700, 0x7CE1, 107, DORMANT_MOTIF),
		GeneScript.new(8, 4500, 0x8DF2, 108, DORMANT_MOTIF),
		GeneScript.new(9, 6000, 0xAF14, 109, DORMANT_MOTIF),
		GeneScript.new(10, HIGH_PROMOTER, OXYGEN_REPRESSOR_SIGNATURE, 110, DORMANT_MOTIF),
		GeneScript.new(11, HIGH_PROMOTER, NEUTRAL_BURDEN_SIGNATURE, 111, DORMANT_MOTIF),
		GeneScript.new(12, HIGH_PROMOTER, NEUTRAL_BURDEN_SIGNATURE, 112, DORMANT_MOTIF),
		GeneScript.new(13, HIGH_PROMOTER, NEUTRAL_BURDEN_SIGNATURE, 113, DORMANT_MOTIF)
	]
	var constitutive = GenomeScript.new(common)
	var responsive = constitutive.deep_copy()
	for locus_variant in BURDEN_LOCI:
		responsive.get_gene_by_locus(int(locus_variant)).regulatory_signature = OXYGEN_REPRESSOR_SIGNATURE
	responsive.validate()
	return {"constitutive": constitutive, "responsive": responsive}

static func run_lineage(
	seed: int,
	genotype_name: String,
	condition: String,
	phase_offset_ticks: int = 0,
	max_ticks: int = MAX_TICKS
) -> Dictionary:
	assert(genotype_name == "constitutive" or genotype_name == "responsive")
	var config = BaseExperiment.create_config(seed)
	config.max_cells = 1
	var sim = SimulationEngineScript.new(config)
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2(1.5, 1.5), config.ancestor_volume)
	cell.genome = create_genomes()[genotype_name].deep_copy()
	cell.initialize_molecular_state(config)
	sim.cells = [cell]
	sim.next_cell_id = 2
	var ticks: int = 0
	var reached_division: bool = false
	for tick in range(max_ticks):
		var oxygen_value: float = BaseGxe.HIGH_OXYGEN
		if condition == CONDITION_HIGH_ANOXIC:
			var phase: int = ((tick + phase_offset_ticks) / BaseGxe.PHASE_TICKS) % 2
			oxygen_value = BaseGxe.HIGH_OXYGEN if phase == 0 else BaseGxe.LOW_OXYGEN
		_maintain_environment(sim, oxygen_value)
		sim.step(1)
		ticks += 1
		if sim.cells.is_empty():
			break
		cell = sim.cells[0]
		if cell.ready_to_divide(config):
			reached_division = true
			break
	var elapsed_min: float = float(ticks) * float(config.tick_dt_min)
	var growth_rate: float = 0.0
	if reached_division and elapsed_min > 0.0:
		growth_rate = log(2.0) / elapsed_min
	return {
		"seed": seed,
		"genotype": genotype_name,
		"condition": condition,
		"ticks": ticks,
		"time_min": elapsed_min,
		"reached_division": reached_division,
		"growth_rate": growth_rate
	}

static func paired_differential(seed: int) -> Dictionary:
	var stable_c: Dictionary = run_lineage(seed, "constitutive", CONDITION_HIGH_STABLE)
	var stable_r: Dictionary = run_lineage(seed, "responsive", CONDITION_HIGH_STABLE)
	var stable_adv: float = float(stable_r["growth_rate"]) - float(stable_c["growth_rate"])
	var reference_rate: float = maxf(1e-12, 0.5 * (float(stable_r["growth_rate"]) + float(stable_c["growth_rate"])))
	var offsets: Array = [0, BaseGxe.PHASE_TICKS / 2, BaseGxe.PHASE_TICKS, (3 * BaseGxe.PHASE_TICKS) / 2]
	var advantages: Array = []
	var all_divided: bool = bool(stable_c["reached_division"]) and bool(stable_r["reached_division"])
	for offset_variant in offsets:
		var offset: int = int(offset_variant)
		var fc: Dictionary = run_lineage(seed, "constitutive", CONDITION_HIGH_ANOXIC, offset)
		var fr: Dictionary = run_lineage(seed, "responsive", CONDITION_HIGH_ANOXIC, offset)
		all_divided = all_divided and bool(fc["reached_division"]) and bool(fr["reached_division"])
		advantages.append(float(fr["growth_rate"]) - float(fc["growth_rate"]))
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
	var total: float = 0.0
	for value_variant in values:
		total += float(value_variant)
	return 0.0 if values.is_empty() else total / float(values.size())
