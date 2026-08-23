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
	_test_reciprocal_w1_w2_exchange_benefits_both_lineages()
	_test_each_partner_product_is_causally_required_for_other_recovery()
	_test_reciprocal_transport_conserves_both_exchanged_molecules()
	_test_same_state_replays_reciprocal_exchange_exactly()

	if failures == 0:
		print("PASS: %d M7 reciprocal cross-feeding tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M7 reciprocal cross-feeding tests failed" % [failures, tests_run])
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
	config.secondary_extracellular_diffusion["W1"] = 0.0
	config.secondary_extracellular_diffusion["W2"] = 0.0
	return config

func _test_reciprocal_w1_w2_exchange_benefits_both_lineages() -> void:
	var result: Dictionary = _run_chain(551109, true, true)
	_assert_true(float(result["a_production_flux"]) > 0.0, "lineage A produces W1 through R04")
	_assert_true(float(result["b_production_flux"]) > 0.0, "lineage B produces W2 through R10")
	_assert_true(float(result["a_export_w1"]) < 0.0, "lineage A exports its W1 product through generic reversible transport")
	_assert_true(float(result["b_export_w2"]) < 0.0, "lineage B exports its W2 product through generic reversible transport")
	_assert_close(float(result["a_first_w2_exchange"]), 0.0, 1e-12, "A cannot consume B's zero-time W2 export from the same transport snapshot")
	_assert_close(float(result["b_first_w1_exchange"]), 0.0, 1e-12, "B cannot consume A's zero-time W1 export from the same transport snapshot")
	_assert_true(float(result["a_import_w2"]) > 0.0, "lineage A imports B-derived W2 on the subsequent membrane phase")
	_assert_true(float(result["b_import_w1"]) > 0.0, "lineage B imports A-derived W1 on the subsequent membrane phase")
	_assert_true(float(result["a_recovery_flux"]) > 0.0, "lineage A converts partner-derived W2 through R11")
	_assert_true(float(result["b_recovery_flux"]) > 0.0, "lineage B converts partner-derived W1 through R05")
	_assert_true(float(result["a_c2_gain"]) > 0.0 and float(result["b_c2_gain"]) > 0.0, "both lineages gain explicit downstream carbon from the other's metabolic product")
	_assert_true(float(result["a_atp_gain"]) > 0.0 and float(result["b_atp_gain"]) > 0.0, "both lineages gain explicit ATP through ordinary recovery stoichiometry")
	_assert_true(float(result["a_ros_gain"]) > 0.0, "W2 oxidative recovery retains its explicit ROS side effect instead of granting a free reciprocal benefit")

func _test_each_partner_product_is_causally_required_for_other_recovery() -> void:
	var without_w2: Dictionary = _run_chain(551109, true, false)
	_assert_true(float(without_w2["b_recovery_flux"]) > 0.0, "B can still recover A-derived W1 when B does not produce W2")
	_assert_close(float(without_w2["a_import_w2"]), 0.0, 1e-12, "A receives no W2 when B's W2-producing reaction is omitted")
	_assert_close(float(without_w2["a_recovery_flux"]), 0.0, 1e-12, "removing B's W2 production specifically abolishes A's reciprocal recovery")
	_assert_close(float(without_w2["a_c2_gain"]), 0.0, 1e-12, "A gains no partner-derived C2 when the corresponding partner product is absent")

	var without_w1: Dictionary = _run_chain(551109, false, true)
	_assert_true(float(without_w1["a_recovery_flux"]) > 0.0, "A can still recover B-derived W2 when A does not produce W1")
	_assert_close(float(without_w1["b_import_w1"]), 0.0, 1e-12, "B receives no W1 when A's W1-producing reaction is omitted")
	_assert_close(float(without_w1["b_recovery_flux"]), 0.0, 1e-12, "removing A's W1 production specifically abolishes B's reciprocal recovery")
	_assert_close(float(without_w1["b_c2_gain"]), 0.0, 1e-12, "B gains no partner-derived C2 when the corresponding partner product is absent")

func _test_reciprocal_transport_conserves_both_exchanged_molecules() -> void:
	var result: Dictionary = _run_chain(771005, true, true)
	_assert_close(float(result["w1_after_transport"]), float(result["w1_before_transport"]), 1e-9, "reciprocal membrane exchange conserves total W1 before recovery")
	_assert_close(float(result["w2_after_transport"]), float(result["w2_before_transport"]), 1e-9, "reciprocal membrane exchange conserves total W2 before recovery")

func _test_same_state_replays_reciprocal_exchange_exactly() -> void:
	var first: Dictionary = _run_chain(902177, true, true)
	var second: Dictionary = _run_chain(902177, true, true)
	_assert_true(first["transport_1"] == second["transport_1"], "same state reproduces exact reciprocal export ledger")
	_assert_true(first["transport_2"] == second["transport_2"], "same state reproduces exact reciprocal import ledger")
	_assert_close(float(first["checksum"]), float(second["checksum"]), 1e-12, "same state and seed reproduce exact reciprocal cross-feeding checksum")

func _run_chain(seed: int, produce_w1: bool, produce_w2: bool) -> Dictionary:
	var config = _config()
	config.seed = seed
	var sim = SimulationEngineScript.new(config)
	# Separate finite disks sampling one common lattice site: this isolates the
	# reciprocal molecular handoff. Spatial attenuation was already gated in M7-C.
	var lineage_a = sim.seed_ancestor(Vector2(6.0, 5.51))
	var lineage_b = sim.seed_ancestor(Vector2(6.0, 6.49))
	_prepare_lineage_a(lineage_a, config)
	_prepare_lineage_b(lineage_b, config)
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var r04 = ReactionCatalogScript.by_id(reactions, "R04")
	var r05 = ReactionCatalogScript.by_id(reactions, "R05")
	var r10 = ReactionCatalogScript.by_id(reactions, "R10")
	var r11 = ReactionCatalogScript.by_id(reactions, "R11")

	var a_production_flux: float = 0.0
	var b_production_flux: float = 0.0
	if produce_w1:
		a_production_flux = float(MetabolicSolverScript.step(lineage_a.metabolites, lineage_a.genome, lineage_a.expression_state, [r04], 1.0, lineage_a.volume, config)["R04"])
	if produce_w2:
		b_production_flux = float(MetabolicSolverScript.step(lineage_b.metabolites, lineage_b.genome, lineage_b.expression_state, [r10], 1.0, lineage_b.volume, config)["R10"])

	var w1_before_transport: float = _system_molecule_total(sim, "W1")
	var w2_before_transport: float = _system_molecule_total(sim, "W2")
	var transport_1: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var a_export_w1: float = float(transport_1["by_cell"][lineage_a.id]["exchange"]["W1"])
	var b_export_w2: float = float(transport_1["by_cell"][lineage_b.id]["exchange"]["W2"])
	var a_first_w2_exchange: float = float(transport_1["by_cell"][lineage_a.id]["exchange"]["W2"])
	var b_first_w1_exchange: float = float(transport_1["by_cell"][lineage_b.id]["exchange"]["W1"])

	var transport_2: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var a_import_w2: float = maxf(0.0, float(transport_2["by_cell"][lineage_a.id]["exchange"]["W2"]))
	var b_import_w1: float = maxf(0.0, float(transport_2["by_cell"][lineage_b.id]["exchange"]["W1"]))
	var w1_after_transport: float = _system_molecule_total(sim, "W1")
	var w2_after_transport: float = _system_molecule_total(sim, "W2")

	var a_c2_before: float = lineage_a.pool("C2")
	var b_c2_before: float = lineage_b.pool("C2")
	var a_atp_before: float = lineage_a.pool("ATP")
	var b_atp_before: float = lineage_b.pool("ATP")
	var a_ros_before: float = lineage_a.pool("ROS")
	var a_recovery_flux: float = float(MetabolicSolverScript.step(lineage_a.metabolites, lineage_a.genome, lineage_a.expression_state, [r11], 1.0, lineage_a.volume, config)["R11"])
	var b_recovery_flux: float = float(MetabolicSolverScript.step(lineage_b.metabolites, lineage_b.genome, lineage_b.expression_state, [r05], 1.0, lineage_b.volume, config)["R05"])

	return {
		"a_production_flux": a_production_flux,
		"b_production_flux": b_production_flux,
		"a_export_w1": a_export_w1,
		"b_export_w2": b_export_w2,
		"a_first_w2_exchange": a_first_w2_exchange,
		"b_first_w1_exchange": b_first_w1_exchange,
		"a_import_w2": a_import_w2,
		"b_import_w1": b_import_w1,
		"a_recovery_flux": a_recovery_flux,
		"b_recovery_flux": b_recovery_flux,
		"a_c2_gain": lineage_a.pool("C2") - a_c2_before,
		"b_c2_gain": lineage_b.pool("C2") - b_c2_before,
		"a_atp_gain": lineage_a.pool("ATP") - a_atp_before,
		"b_atp_gain": lineage_b.pool("ATP") - b_atp_before,
		"a_ros_gain": lineage_a.pool("ROS") - a_ros_before,
		"w1_before_transport": w1_before_transport,
		"w1_after_transport": w1_after_transport,
		"w2_before_transport": w2_before_transport,
		"w2_after_transport": w2_after_transport,
		"transport_1": transport_1,
		"transport_2": transport_2,
		"checksum": sim.checksum()
	}

func _prepare_lineage_a(cell, config) -> void:
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	_install_exact_protein(cell, 4, int(ReactionCatalogScript.by_id(reactions, "R04").signature), config)
	_install_exact_protein(cell, 9, int(ReactionCatalogScript.by_id(reactions, "R11").signature), config)
	_install_exact_protein(cell, 11, MembraneTransportScript.target_signature("W1"), config)
	_install_exact_protein(cell, 12, MembraneTransportScript.target_signature("W2"), config)
	cell.set_pool("C3", 4.0)
	cell.set_pool("W1", 0.0)
	cell.set_pool("W2", 0.0)
	cell.set_pool("C2", 0.0)
	cell.set_pool("NADH", 4.0)
	cell.set_pool("NAD", 2.0)
	cell.set_pool("ADP", 10.0)
	cell.set_pool("ATP", 4.0)
	cell.set_pool("O2", 4.0)
	cell.set_pool("ROS", 0.0)

func _prepare_lineage_b(cell, config) -> void:
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	_install_exact_protein(cell, 4, int(ReactionCatalogScript.by_id(reactions, "R05").signature), config)
	_install_exact_protein(cell, 9, int(ReactionCatalogScript.by_id(reactions, "R10").signature), config)
	_install_exact_protein(cell, 11, MembraneTransportScript.target_signature("W1"), config)
	_install_exact_protein(cell, 12, MembraneTransportScript.target_signature("W2"), config)
	cell.set_pool("C3", 4.0)
	cell.set_pool("W1", 0.0)
	cell.set_pool("W2", 0.0)
	cell.set_pool("C2", 0.0)
	cell.set_pool("NADH", 4.0)
	cell.set_pool("NAD", 2.0)
	cell.set_pool("ADP", 10.0)
	cell.set_pool("ATP", 4.0)
	cell.set_pool("O2", 4.0)
	cell.set_pool("ROS", 0.0)

func _install_exact_protein(cell, locus_id: int, signature: int, config) -> void:
	var gene = cell.genome.get_gene_by_locus(locus_id)
	assert(gene != null)
	gene.protein_signature = signature & 0xFFFF
	var amount: float = maxf(1.0, float(config.expression_reference_protein_count) * float(gene.promoter_strength()))
	cell.expression_state[locus_id]["protein"] = {signature & 0xFFFF: amount}

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
