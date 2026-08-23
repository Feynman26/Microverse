extends RefCounted
class_name M5CRosCircuit

const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")
const BaseGxe = preload("res://src/experiments/m5c_gxe_assay.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")

const CONDITION_HIGH_STABLE: String = BaseGxe.CONDITION_HIGH_STABLE
const CONDITION_HIGH_ANOXIC: String = BaseGxe.CONDITION_HIGH_ANOXIC
const OXPHOS_LOCUS: int = 3
const ROS_CONTROL_LOCUS: int = 8
const REGULATOR_LOCUS: int = 10
const OXPHOS_PROTEIN_SIGNATURE: int = 0x369C
const ROS_CONTROL_PROTEIN_SIGNATURE: int = 0x8DF2
const OXYGEN_REPRESSOR_SIGNATURE: int = 0xCCCC
const DORMANT_MOTIF: int = BaseExperiment.DORMANT_MOTIF
const HIGH_PROMOTER: int = 10000

# Final M5-C circuit redesign after the R03-targeted hypotheses were falsified.
# Both genomes pay for the same proteins and differ in one inherited field only:
# the regulatory motif of the ROS-control enzyme locus. O2 already causes R03
# to generate ROS, and R09 already consumes ROS + ATP. No new fitness term or
# assay-specific physiology is introduced here.
static func create_genomes() -> Dictionary:
	var common: Array = [
		GeneScript.new(1, 6200, 0x1357, 101, DORMANT_MOTIF),
		GeneScript.new(2, 5400, 0x2468, 102, DORMANT_MOTIF),
		GeneScript.new(OXPHOS_LOCUS, HIGH_PROMOTER, OXPHOS_PROTEIN_SIGNATURE, 103, DORMANT_MOTIF),
		GeneScript.new(4, 4800, 0x48AD, 104, DORMANT_MOTIF),
		GeneScript.new(5, 3900, 0x5ACE, 105, DORMANT_MOTIF),
		GeneScript.new(6, 6600, 0x6BDF, 106, DORMANT_MOTIF),
		GeneScript.new(7, 5700, 0x7CE1, 107, DORMANT_MOTIF),
		GeneScript.new(ROS_CONTROL_LOCUS, HIGH_PROMOTER, ROS_CONTROL_PROTEIN_SIGNATURE, 108, DORMANT_MOTIF),
		GeneScript.new(9, 6000, 0xAF14, 109, DORMANT_MOTIF),
		GeneScript.new(REGULATOR_LOCUS, HIGH_PROMOTER, OXYGEN_REPRESSOR_SIGNATURE, 110, DORMANT_MOTIF)
	]
	var constitutive = GenomeScript.new(common)
	var responsive = constitutive.deep_copy()
	responsive.get_gene_by_locus(ROS_CONTROL_LOCUS).regulatory_signature = OXYGEN_REPRESSOR_SIGNATURE
	responsive.validate()
	return {"constitutive": constitutive, "responsive": responsive}

static func run_lineage(
	seed: int,
	genotype_name: String,
	condition: String,
	phase_offset_ticks: int = 0,
	max_ticks: int = 7200
) -> Dictionary:
	assert(genotype_name == "constitutive" or genotype_name == "responsive")
	assert(condition == CONDITION_HIGH_STABLE or condition == CONDITION_HIGH_ANOXIC)
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
	var peak_damage: float = 0.0
	for tick in range(max_ticks):
		var oxygen_value: float = BaseGxe.HIGH_OXYGEN
		if condition == CONDITION_HIGH_ANOXIC:
			var phase: int = ((tick + phase_offset_ticks) / BaseGxe.PHASE_TICKS) % 2
			oxygen_value = BaseGxe.HIGH_OXYGEN if phase == 0 else BaseGxe.ANOXIC_OXYGEN
		_maintain_environment(sim, oxygen_value)
		sim.step(1)
		realized_ticks += 1
		if sim.cells.is_empty():
			break
		cell = sim.cells[0]
		peak_damage = maxf(peak_damage, float(cell.damage))
		if cell.ready_to_divide(config):
			reached_division = true
			break

	var elapsed_min: float = float(realized_ticks) * float(config.tick_dt_min)
	var growth_rate: float = 0.0
	var final_damage: float = peak_damage
	var death_reason: String = ""
	if not sim.cells.is_empty():
		cell = sim.cells[0]
		final_damage = float(cell.damage)
		death_reason = String(cell.death_reason)
	if reached_division and elapsed_min > 0.0:
		growth_rate = log(2.0) / elapsed_min
	return {
		"seed": seed,
		"genotype": genotype_name,
		"condition": condition,
		"ticks": realized_ticks,
		"time_min": elapsed_min,
		"reached_division": reached_division,
		"growth_rate": growth_rate,
		"peak_damage": peak_damage,
		"final_damage": final_damage,
		"death_reason": death_reason
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
