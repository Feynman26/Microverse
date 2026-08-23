extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const ExtracellularCatalysisScript = preload("res://src/chemistry/extracellular_catalysis.gd")
const ExtracellularReactionCatalogScript = preload("res://src/chemistry/extracellular_reaction_catalog.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

const DAMAGE_PROTEIN_SIGNATURE: int = 0xDE03
const INACTIVE_SIGNATURE: int = 0x0000

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_damaging_chemistry_is_sequence_derived_and_ancestrally_dormant()
	_test_secreted_damaging_chemistry_harms_permeable_neighbor()
	_test_ordinary_intracellular_detox_provides_self_resistance_at_atp_cost()
	_test_biosynthetic_loss_creates_external_metabolite_dependence()
	_test_secondary_metabolite_uptake_can_restore_growth_and_division_readiness()
	_test_closure_scenarios_replay_exactly()

	if failures == 0:
		print("PASS: %d compact M7 closure tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d compact M7 closure tests failed" % [failures, tests_run])
		quit(1)

func _config():
	var config = SimConfigScript.new()
	config.world_width = 16
	config.world_height = 12
	config.max_cells = 8
	config.mutation_enabled = false
	config.initial_glucose = 0.0
	config.initial_oxygen = 0.0
	config.initial_nitrogen = 0.0
	config.initial_phosphorus = 0.0
	config.glucose_diffusion = 0.0
	config.oxygen_diffusion = 0.0
	config.nitrogen_diffusion = 0.0
	config.phosphorus_diffusion = 0.0
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		config.secondary_extracellular_initial[metabolite_id] = 0.0
		config.secondary_extracellular_diffusion[metabolite_id] = 0.0
	config.extracellular_protein_diffusion = 0.0
	config.extracellular_protein_secretion_fraction_per_min = 1.0
	config.extracellular_catalysis_rate_scale = 4.0
	config.spontaneous_ros_decay_rate = 0.0
	config.basal_repair_rate = 0.0
	config.ros_damage_rate = 1.0
	return config

func _test_damaging_chemistry_is_sequence_derived_and_ancestrally_dormant() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var locus9 = ancestor.get_gene_by_locus(9)
	_assert_true(int(locus9.protein_signature) == 0x9E03, "antagonism precursor starts from canonical ancestral locus-9 sequence")
	_assert_true(CatalyticLandscapeScript.hamming_distance(int(locus9.protein_signature), DAMAGE_PROTEIN_SIGNATURE) == 1, "one coding bit creates secretion-compatible DE03")
	_assert_true(ExtracellularCatalysisScript.has_secretion_signal(DAMAGE_PROTEIN_SIGNATURE), "DE03 uses the existing generic Dxxx secretion motif")

	var reactions: Array = ExtracellularReactionCatalogScript.create_m7_candidate()
	var e03 = ExtracellularReactionCatalogScript.by_id(reactions, "E03")
	_assert_true(e03 != null, "M7 catalog contains E03 extracellular oxidant generation")
	_assert_true(CatalyticLandscapeScript.hamming_distance(DAMAGE_PROTEIN_SIGNATURE, int(e03.signature)) == 4, "DE03 lies inside the ordinary E03 catalytic radius")
	var ancestral_dormant: bool = true
	for gene in ancestor.genes:
		if CatalyticLandscapeScript.hamming_distance(int(gene.protein_signature), int(e03.signature)) <= CatalyticLandscapeScript.ACTIVE_MAX_DISTANCE:
			ancestral_dormant = false
			break
	_assert_true(ancestral_dormant, "unmodified ancestor has no E03 extracellular damaging activity")
	var e01 = ExtracellularReactionCatalogScript.by_id(reactions, "E01")
	var e02 = ExtracellularReactionCatalogScript.by_id(reactions, "E02")
	_assert_true(CatalyticLandscapeScript.hamming_distance(DAMAGE_PROTEIN_SIGNATURE, int(e01.signature)) > CatalyticLandscapeScript.ACTIVE_MAX_DISTANCE, "DE03 does not receive public-resource E01 activity for free")
	_assert_true(CatalyticLandscapeScript.hamming_distance(DAMAGE_PROTEIN_SIGNATURE, int(e02.signature)) > CatalyticLandscapeScript.ACTIVE_MAX_DISTANCE, "DE03 does not receive detox E02 activity for free")

func _test_secreted_damaging_chemistry_harms_permeable_neighbor() -> void:
	var exposure: Dictionary = _run_antagonism(true, false)
	var control: Dictionary = _run_antagonism(false, false)
	_assert_true(float(exposure["secreted"]) > 0.0, "producer secretes a real DE03 protein cohort")
	_assert_true(float(exposure["secretion_atp"]) > 0.0, "producer pays ATP for damaging-protein secretion")
	_assert_true(float(exposure["e03_flux"]) > 0.0, "secreted DE03 converts neutral X into extracellular oxidant")
	_assert_close(float(exposure["x_consumed"]), float(exposure["oxidant_created"]), 1e-10, "E03 preserves the tracked one-to-one X/oxidant ledger")
	_assert_true(float(exposure["victim_ros_import"]) > 0.0, "permeable neighbor imports producer-generated extracellular oxidant")
	_assert_true(float(exposure["victim_damage"]) > 0.0, "imported oxidant produces ordinary intracellular ROS damage")
	_assert_close(float(control["e03_flux"]), 0.0, 1e-12, "without damaging catalyst the same neutral-X environment produces no E03 flux")
	_assert_close(float(control["victim_damage"]), 0.0, 1e-12, "without producer chemistry the neighbor receives no damage")

func _test_ordinary_intracellular_detox_provides_self_resistance_at_atp_cost() -> void:
	var resistant: Dictionary = _run_antagonism(true, true)
	_assert_true(float(resistant["producer_ros_import"]) > 0.0, "producer can be exposed to its own oxidant through ordinary ROS permeability")
	_assert_true(float(resistant["producer_r09_flux"]) > 0.0, "ordinary R09 chemistry detoxifies imported ROS inside the producer")
	_assert_true(float(resistant["producer_r09_atp_cost"]) > 0.0, "self-resistance spends ATP through R09 stoichiometry")
	_assert_true(float(resistant["producer_damage"]) < float(resistant["victim_damage"]), "R09-positive producer accumulates less damage than equally permeable R09-disabled neighbor")

func _run_antagonism(with_damage_catalyst: bool, producer_permeable: bool) -> Dictionary:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var producer = sim.seed_ancestor(Vector2(6.0, 5.51))
	var victim = sim.seed_ancestor(Vector2(6.0, 6.49))
	producer.set_pool("ATP", 12.0)
	producer.set_pool("ADP", 4.0)
	victim.set_pool("ATP", 12.0)
	victim.set_pool("ADP", 4.0)
	producer.set_pool("ROS", 0.0)
	victim.set_pool("ROS", 0.0)
	producer.damage = 0.0
	victim.damage = 0.0

	if with_damage_catalyst:
		_install_exact_protein(producer, 9, DAMAGE_PROTEIN_SIGNATURE, config)
	_install_exact_protein(victim, 11, MembraneTransportScript.target_signature("ROS"), config)
	_disable_reaction_protein(victim, 8, config)
	if producer_permeable:
		_install_exact_protein(producer, 11, MembraneTransportScript.target_signature("ROS"), config)

	sim.world.release("neutral_x", Vector2(6.0, 6.0), 2.0)
	var x_before: float = sim.world.get_field("neutral_x").total_amount()
	var oxidant_before: float = sim.world.get_field("oxidant").total_amount()
	var secretion: Dictionary = sim._secrete_extracellular_proteins(1.0)
	var e03 = ExtracellularReactionCatalogScript.by_id(sim.extracellular_reactions, "E03")
	var catalysis: Dictionary = ExtracellularCatalysisScript.step(sim.world, [e03], 1.0, config)
	var x_after: float = sim.world.get_field("neutral_x").total_amount()
	var oxidant_after: float = sim.world.get_field("oxidant").total_amount()
	var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var victim_import: float = maxf(0.0, float(transport["by_cell"][victim.id]["exchange"]["ROS"]))
	var producer_import: float = maxf(0.0, float(transport["by_cell"][producer.id]["exchange"]["ROS"]))

	var r09 = ReactionCatalogScript.by_id(sim.reactions, "R09")
	var producer_atp_before_r09: float = producer.pool("ATP")
	var producer_flux: Dictionary = MetabolicSolverScript.step(producer.metabolites, producer.genome, producer.expression_state, [r09], 1.0, producer.volume, config)
	var victim_flux: Dictionary = MetabolicSolverScript.step(victim.metabolites, victim.genome, victim.expression_state, [r09], 1.0, victim.volume, config)
	producer._update_damage_and_repair(1.0, config)
	victim._update_damage_and_repair(1.0, config)
	var producer_secretion: Dictionary = secretion["by_cell"].get(producer.id, {"total_secreted": 0.0, "atp_spent": 0.0})
	return {
		"secreted": float(producer_secretion.get("total_secreted", 0.0)),
		"secretion_atp": float(producer_secretion.get("atp_spent", 0.0)),
		"e03_flux": float(catalysis["fluxes"]["E03"]),
		"x_consumed": x_before - x_after,
		"oxidant_created": oxidant_after - oxidant_before,
		"victim_ros_import": victim_import,
		"producer_ros_import": producer_import,
		"producer_r09_flux": float(producer_flux["R09"]),
		"victim_r09_flux": float(victim_flux["R09"]),
		"producer_r09_atp_cost": producer_atp_before_r09 - producer.pool("ATP"),
		"producer_damage": producer.damage,
		"victim_damage": victim.damage,
		"checksum": sim.checksum()
	}

func _test_biosynthetic_loss_creates_external_metabolite_dependence() -> void:
	var intact: Dictionary = _run_aa_growth(true, false, false)
	var lost: Dictionary = _run_aa_growth(false, false, false)
	var external_without_transport: Dictionary = _run_aa_growth(false, true, false)
	_assert_true(float(intact["bio_gain"]) > 0.10, "R06-positive cell can synthesize AA and grow without extracellular AA")
	_assert_close(float(lost["bio_gain"]), 0.0, 1e-10, "loss of R06 blocks biomass growth when AA is absent")
	_assert_close(float(external_without_transport["bio_gain"]), 0.0, 1e-10, "extracellular AA alone cannot rescue biosynthetic loss without matching membrane machinery")

func _test_secondary_metabolite_uptake_can_restore_growth_and_division_readiness() -> void:
	var rescued: Dictionary = _run_aa_growth(false, true, true)
	_assert_true(float(rescued["aa_imported"]) > 0.0, "matching generic AA transporter imports extracellular secondary metabolite")
	_assert_true(float(rescued["bio_gain"]) > 0.50, "imported AA is consumed by ordinary R12 and restores structural biomass growth")
	_assert_true(bool(rescued["ready_to_divide"]), "secondary-resource rescue can restore ordinary division readiness")

func _run_aa_growth(biosynthesis_intact: bool, external_aa: bool, aa_transport: bool) -> Dictionary:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	var r06 = ReactionCatalogScript.by_id(sim.reactions, "R06")
	var r12 = ReactionCatalogScript.by_id(sim.reactions, "R12")
	if biosynthesis_intact:
		_install_exact_protein(cell, 5, int(r06.signature), config)
	else:
		_disable_reaction_protein(cell, 5, config)
	_install_exact_protein(cell, 10, int(r12.signature), config)
	if aa_transport:
		_install_exact_protein(cell, 11, MembraneTransportScript.target_signature("AA"), config)

	cell.set_pool("BIO", 1.0)
	cell._sync_volume_from_biomass(config)
	cell.set_pool("ATP", 40.0)
	cell.set_pool("ADP", 10.0)
	cell.set_pool("C2", 20.0)
	cell.set_pool("NH4", 20.0)
	cell.set_pool("AA", 0.0)
	cell.set_pool("LIP", 10.0)
	cell.set_pool("NUC", 20.0)
	if external_aa:
		sim.world.release("amino_acids", cell.position, 10.0)

	var bio_before: float = cell.pool("BIO")
	var aa_imported: float = 0.0
	for _i in range(80):
		var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
		aa_imported += maxf(0.0, float(transport["by_cell"][cell.id]["exchange"]["AA"]))
		MetabolicSolverScript.step(cell.metabolites, cell.genome, cell.expression_state, [r06, r12], config.tick_dt_min, cell.volume, config)
		cell._sync_volume_from_biomass(config)
	return {
		"bio_gain": cell.pool("BIO") - bio_before,
		"aa_imported": aa_imported,
		"ready_to_divide": cell.ready_to_divide(config),
		"checksum": sim.checksum()
	}

func _test_closure_scenarios_replay_exactly() -> void:
	var first_antagonism: Dictionary = _run_antagonism(true, true)
	var second_antagonism: Dictionary = _run_antagonism(true, true)
	_assert_true(first_antagonism == second_antagonism, "same state reproduces exact antagonism/self-resistance outcome")
	var first_rescue: Dictionary = _run_aa_growth(false, true, true)
	var second_rescue: Dictionary = _run_aa_growth(false, true, true)
	_assert_true(first_rescue == second_rescue, "same state reproduces exact auxotrophic rescue/growth outcome")

func _install_exact_protein(cell, locus_id: int, signature: int, config) -> void:
	var gene = cell.genome.get_gene_by_locus(locus_id)
	assert(gene != null)
	gene.protein_signature = signature & 0xFFFF
	cell.expression_state[locus_id]["protein"] = {signature & 0xFFFF: float(config.expression_reference_protein_count)}

func _disable_reaction_protein(cell, locus_id: int, config) -> void:
	_install_exact_protein(cell, locus_id, INACTIVE_SIGNATURE, config)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
