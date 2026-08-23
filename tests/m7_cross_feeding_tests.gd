extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_w1_one_way_cross_feeding_chain()
	_test_consumer_without_donor_has_no_w1_recovery()
	_test_donor_transport_is_required_for_transfer()
	_test_diffusion_distance_changes_transfer_magnitude()
	_test_same_state_replays_cross_feeding_exactly()

	if failures == 0:
		print("PASS: %d M7 one-way cross-feeding tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M7 one-way cross-feeding tests failed" % [failures, tests_run])
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
	return config

func _reactions() -> Array:
	return ReactionCatalogScript.create_m4_candidate()

func _test_w1_one_way_cross_feeding_chain() -> void:
	var config = _config()
	config.secondary_extracellular_diffusion["W1"] = 0.0
	var sim = SimulationEngineScript.new(config)
	# These positions occupy distinct physical disks but round to the same lattice
	# site, isolating biochemical transfer from diffusion in this first assay.
	var donor = sim.seed_ancestor(Vector2(6.0, 5.51))
	var recipient = sim.seed_ancestor(Vector2(6.0, 6.49))
	_prepare_donor(donor, config)
	_prepare_recipient(recipient, config)

	var reactions: Array = _reactions()
	var r04 = ReactionCatalogScript.by_id(reactions, "R04")
	var r05 = ReactionCatalogScript.by_id(reactions, "R05")
	var donor_flux: Dictionary = MetabolicSolverScript.step(donor.metabolites, donor.genome, donor.expression_state, [r04], 1.0, donor.volume, config)
	_assert_true(float(donor_flux["R04"]) > 0.0, "donor converts intracellular C3 into W1 through ordinary R04 chemistry")
	_assert_true(donor.pool("W1") > 0.0, "R04 creates a physically explicit intracellular waste pool")

	var w1_field: String = MetaboliteCatalogScript.extracellular_field("W1")
	var conserved_before_transfer: float = _system_molecule_total(sim, "W1")
	var first_transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var donor_export: float = float(first_transport["by_cell"][donor.id]["exchange"]["W1"])
	var first_recipient_exchange: float = float(first_transport["by_cell"][recipient.id]["exchange"]["W1"])
	_assert_true(donor_export < 0.0, "W1-producing donor exports W1 down its concentration gradient")
	_assert_close(first_recipient_exchange, 0.0, 1e-12, "recipient cannot consume donor export in zero time from the same snapshot")
	_assert_true(sim.world.get_field(w1_field).total_amount() > 0.0, "donor export creates extracellular W1 without a secretion behavior API")

	var second_transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var recipient_import: float = float(second_transport["by_cell"][recipient.id]["exchange"]["W1"])
	_assert_true(recipient_import > 0.0 and recipient.pool("W1") > 0.0, "recipient imports donor-derived extracellular W1 through compatible protein")
	_assert_close(_system_molecule_total(sim, "W1"), conserved_before_transfer, 1e-9, "two membrane phases conserve W1 across donor, recipient, and world")

	var c2_before: float = recipient.pool("C2")
	var atp_before: float = recipient.pool("ATP")
	var w1_before_recovery: float = recipient.pool("W1")
	var recipient_flux: Dictionary = MetabolicSolverScript.step(recipient.metabolites, recipient.genome, recipient.expression_state, [r05], 1.0, recipient.volume, config)
	_assert_true(float(recipient_flux["R05"]) > 0.0, "recipient recovers imported donor W1 through ordinary R05 catalysis")
	_assert_true(recipient.pool("W1") < w1_before_recovery, "recipient R05 consumes the transferred molecule")
	_assert_true(recipient.pool("C2") > c2_before, "cross-fed W1 becomes downstream intracellular carbon rather than an ecological score")
	_assert_true(recipient.pool("ATP") > atp_before, "cross-fed W1 yields an explicit energetic benefit through reaction stoichiometry")

func _test_consumer_without_donor_has_no_w1_recovery() -> void:
	var config = _config()
	config.secondary_extracellular_diffusion["W1"] = 0.0
	var sim = SimulationEngineScript.new(config)
	var recipient = sim.seed_ancestor(Vector2(6.0, 6.0))
	_prepare_recipient(recipient, config)
	var reactions: Array = _reactions()
	var r05 = ReactionCatalogScript.by_id(reactions, "R05")
	var c2_before: float = recipient.pool("C2")
	var atp_before: float = recipient.pool("ATP")

	sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var flux: Dictionary = MetabolicSolverScript.step(recipient.metabolites, recipient.genome, recipient.expression_state, [r05], 1.0, recipient.volume, config)
	_assert_close(float(flux["R05"]), 0.0, 1e-12, "same recipient without donor-derived W1 has zero recovery flux")
	_assert_close(recipient.pool("C2"), c2_before, 1e-12, "recipient-alone control gains no downstream carbon from absent W1")
	_assert_close(recipient.pool("ATP"), atp_before, 1e-12, "recipient-alone control gains no W1-derived ATP")

func _test_donor_transport_is_required_for_transfer() -> void:
	var config = _config()
	config.secondary_extracellular_diffusion["W1"] = 0.0
	var sim = SimulationEngineScript.new(config)
	var donor = sim.seed_ancestor(Vector2(6.0, 5.51))
	var recipient = sim.seed_ancestor(Vector2(6.0, 6.49))
	_prepare_donor(donor, config, false)
	_prepare_recipient(recipient, config)
	var reactions: Array = _reactions()
	var r04 = ReactionCatalogScript.by_id(reactions, "R04")
	var r05 = ReactionCatalogScript.by_id(reactions, "R05")
	var donor_flux: Dictionary = MetabolicSolverScript.step(donor.metabolites, donor.genome, donor.expression_state, [r04], 1.0, donor.volume, config)
	_assert_true(float(donor_flux["R04"]) > 0.0 and donor.pool("W1") > 0.0, "transport-disabled donor still produces intracellular W1")

	var field_name: String = MetaboliteCatalogScript.extracellular_field("W1")
	sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	_assert_close(sim.world.get_field(field_name).total_amount(), 0.0, 1e-12, "without compatible donor membrane protein, produced W1 remains intracellular")
	_assert_close(recipient.pool("W1"), 0.0, 1e-12, "recipient cannot receive W1 that donor cannot export")
	var recipient_flux: Dictionary = MetabolicSolverScript.step(recipient.metabolites, recipient.genome, recipient.expression_state, [r05], 1.0, recipient.volume, config)
	_assert_close(float(recipient_flux["R05"]), 0.0, 1e-12, "disabling donor transport causally abolishes recipient W1 recovery")

func _test_diffusion_distance_changes_transfer_magnitude() -> void:
	var near_result: Dictionary = _run_spatial_transfer(Vector2(6.0, 6.0))
	var far_result: Dictionary = _run_spatial_transfer(Vector2(12.0, 6.0))
	_assert_true(float(near_result["imported"]) > 0.0, "near recipient receives donor-derived W1 after ordinary field diffusion")
	_assert_true(float(near_result["imported"]) > float(far_result["imported"]) + 1e-8, "greater physical separation weakens transfer through diffusion alone")
	_assert_true(float(near_result["local_w1"]) > float(far_result["local_w1"]), "distance effect is explained by measured local W1 concentration rather than a distance bonus rule")

func _run_spatial_transfer(recipient_position: Vector2) -> Dictionary:
	var config = _config()
	config.secondary_extracellular_diffusion["W1"] = 0.60
	var sim = SimulationEngineScript.new(config)
	var donor = sim.seed_ancestor(Vector2(4.0, 6.0))
	var recipient = sim.seed_ancestor(recipient_position)
	_prepare_donor(donor, config)
	_prepare_recipient(recipient, config)
	var r04 = ReactionCatalogScript.by_id(_reactions(), "R04")
	MetabolicSolverScript.step(donor.metabolites, donor.genome, donor.expression_state, [r04], 1.0, donor.volume, config)
	var export_summary: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	assert(float(export_summary["by_cell"][donor.id]["exchange"]["W1"]) < 0.0)
	for _i in range(16):
		sim.world.diffuse(config.tick_dt_min)
	var field_name: String = MetaboliteCatalogScript.extracellular_field("W1")
	var local_before: float = sim.world.sample(field_name, recipient.position)
	var import_summary: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	return {
		"imported": maxf(0.0, float(import_summary["by_cell"][recipient.id]["exchange"]["W1"])),
		"local_w1": local_before,
		"checksum": sim.checksum()
	}

func _test_same_state_replays_cross_feeding_exactly() -> void:
	var first: Dictionary = _run_replay_chain(771331)
	var second: Dictionary = _run_replay_chain(771331)
	_assert_true(first["transport_1"] == second["transport_1"], "same state reproduces exact donor export ledger")
	_assert_true(first["transport_2"] == second["transport_2"], "same state reproduces exact recipient import ledger")
	_assert_true(first["donor_flux"] == second["donor_flux"] and first["recipient_flux"] == second["recipient_flux"], "same state reproduces exact cross-feeding metabolic fluxes")
	_assert_close(float(first["checksum"]), float(second["checksum"]), 1e-12, "same state and seed reproduce exact post-cross-feeding checksum")

func _run_replay_chain(seed: int) -> Dictionary:
	var config = _config()
	config.seed = seed
	config.secondary_extracellular_diffusion["W1"] = 0.0
	var sim = SimulationEngineScript.new(config)
	var donor = sim.seed_ancestor(Vector2(6.0, 5.51))
	var recipient = sim.seed_ancestor(Vector2(6.0, 6.49))
	_prepare_donor(donor, config)
	_prepare_recipient(recipient, config)
	var reactions: Array = _reactions()
	var r04 = ReactionCatalogScript.by_id(reactions, "R04")
	var r05 = ReactionCatalogScript.by_id(reactions, "R05")
	var donor_flux: Dictionary = MetabolicSolverScript.step(donor.metabolites, donor.genome, donor.expression_state, [r04], 1.0, donor.volume, config)
	var transport_1: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var transport_2: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var recipient_flux: Dictionary = MetabolicSolverScript.step(recipient.metabolites, recipient.genome, recipient.expression_state, [r05], 1.0, recipient.volume, config)
	return {
		"donor_flux": donor_flux,
		"transport_1": transport_1,
		"transport_2": transport_2,
		"recipient_flux": recipient_flux,
		"checksum": sim.checksum()
	}

func _prepare_donor(cell, config, install_transport: bool = true) -> void:
	var r04 = ReactionCatalogScript.by_id(_reactions(), "R04")
	_install_exact_protein(cell, 4, int(r04.signature), config)
	if install_transport:
		_install_exact_protein(cell, 12, MembraneTransportScript.target_signature("W1"), config)
	else:
		_clear_secondary_transport_activity(cell, config)
	cell.set_pool("C3", 4.0)
	cell.set_pool("W1", 0.0)
	cell.set_pool("NADH", 4.0)
	cell.set_pool("NAD", 0.0)
	cell.set_pool("ADP", 10.0)
	cell.set_pool("ATP", 4.0)

func _prepare_recipient(cell, config) -> void:
	var r05 = ReactionCatalogScript.by_id(_reactions(), "R05")
	_install_exact_protein(cell, 4, int(r05.signature), config)
	_install_exact_protein(cell, 12, MembraneTransportScript.target_signature("W1"), config)
	cell.set_pool("W1", 0.0)
	cell.set_pool("C2", 0.0)
	cell.set_pool("NAD", 4.0)
	cell.set_pool("NADH", 0.0)
	cell.set_pool("ADP", 10.0)
	cell.set_pool("ATP", 4.0)

func _install_exact_protein(cell, locus_id: int, signature: int, config) -> void:
	var gene = cell.genome.get_gene_by_locus(locus_id)
	assert(gene != null)
	gene.protein_signature = signature & 0xFFFF
	var amount: float = maxf(1.0, float(config.expression_reference_protein_count) * float(gene.promoter_strength()))
	cell.expression_state[locus_id]["protein"] = {signature & 0xFFFF: amount}

func _clear_secondary_transport_activity(cell, config) -> void:
	# Restore every realized protein cohort to the canonical ancestral coding
	# sequence for that locus; M7-B proves this proteome has zero secondary
	# transport activity across all ten molecules.
	var ancestor = preload("res://src/genetics/genome.gd").create_ancestor()
	for gene in cell.genome.genes:
		var ancestral_gene = ancestor.get_gene_by_locus(int(gene.locus_id))
		gene.protein_signature = int(ancestral_gene.protein_signature)
		var amount: float = maxf(1.0, float(config.expression_reference_protein_count) * float(gene.promoter_strength()))
		cell.expression_state[int(gene.locus_id)]["protein"] = {int(ancestral_gene.protein_signature): amount}
	# Reinstall R04 after restoring the donor proteome, while leaving all other
	# loci canonical. R04 catalytic activity itself must not imply W1 transport.
	var r04 = ReactionCatalogScript.by_id(_reactions(), "R04")
	_install_exact_protein(cell, 4, int(r04.signature), config)

func _system_molecule_total(sim, metabolite_id: String) -> float:
	var total: float = sim.world.get_field(MetaboliteCatalogScript.extracellular_field(metabolite_id)).total_amount()
	for cell in sim.cells:
		total += cell.pool(metabolite_id)
	return total

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
