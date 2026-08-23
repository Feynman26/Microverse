extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const ExtracellularCatalysisScript = preload("res://src/chemistry/extracellular_catalysis.gd")
const ExtracellularReactionCatalogScript = preload("res://src/chemistry/extracellular_reaction_catalog.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

const DETOX_PROTEIN_SIGNATURE: int = 0xDACE

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_detox_secretion_and_catalysis_are_one_bit_accessible()
	_test_e02_is_bounded_nonnegative_and_catalytic()
	_test_detox_producer_reduces_neighbor_ros_import_at_explicit_cost()
	_test_reduced_ros_import_reduces_neighbor_damage()
	_test_spatial_separation_weakens_detoxification_benefit_through_field_gradient()
	_test_same_state_replays_detoxification_exactly()

	if failures == 0:
		print("PASS: %d M7 neighborhood-detoxification tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M7 neighborhood-detoxification tests failed" % [failures, tests_run])
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
	config.extracellular_protein_diffusion = 0.0
	config.extracellular_protein_secretion_fraction_per_min = 1.0
	config.extracellular_catalysis_rate_scale = 4.0
	config.spontaneous_ros_decay_rate = 0.0
	config.basal_repair_rate = 0.0
	config.ros_damage_rate = 1.0
	return config

func _test_detox_secretion_and_catalysis_are_one_bit_accessible() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var locus5 = ancestor.get_gene_by_locus(5)
	_assert_true(int(locus5.protein_signature) == 0x5ACE, "detox precursor starts from canonical ancestral locus-5 sequence")
	_assert_true(CatalyticLandscapeScript.hamming_distance(int(locus5.protein_signature), DETOX_PROTEIN_SIGNATURE) == 1, "one ordinary coding bit flip creates DACE from ancestral 5ACE")
	_assert_true(ExtracellularCatalysisScript.has_secretion_signal(DETOX_PROTEIN_SIGNATURE), "DACE enters the same generic Dxxx secretion class used by M7-E")

	var e02 = ExtracellularReactionCatalogScript.by_id(ExtracellularReactionCatalogScript.create_m7_candidate(), "E02")
	_assert_true(e02 != null, "M7 extracellular catalog contains E02 oxidant neutralization")
	_assert_true(CatalyticLandscapeScript.hamming_distance(DETOX_PROTEIN_SIGNATURE, int(e02.signature)) == 4, "secreted DACE lies inside ordinary catalytic radius of E02")
	var ancestral_dormant: bool = true
	for gene in ancestor.genes:
		if CatalyticLandscapeScript.hamming_distance(int(gene.protein_signature), int(e02.signature)) <= CatalyticLandscapeScript.ACTIVE_MAX_DISTANCE:
			ancestral_dormant = false
			break
	_assert_true(ancestral_dormant, "E02 is outside catalytic radius of every ancestral protein")

	var e01 = ExtracellularReactionCatalogScript.by_id(ExtracellularReactionCatalogScript.create_m7_candidate(), "E01")
	_assert_true(CatalyticLandscapeScript.hamming_distance(DETOX_PROTEIN_SIGNATURE, int(e01.signature)) > CatalyticLandscapeScript.ACTIVE_MAX_DISTANCE, "DACE does not receive M7-E lipid-hydrolase activity for free")

func _test_e02_is_bounded_nonnegative_and_catalytic() -> void:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var site := Vector2(6.0, 6.0)
	sim.world.release_protein(DETOX_PROTEIN_SIGNATURE, site, 160.0, config.extracellular_protein_diffusion)
	sim.world.release("oxidant", site, 1.0)
	var e02 = ExtracellularReactionCatalogScript.by_id(sim.extracellular_reactions, "E02")
	var enzyme_before: float = sim.world.total_extracellular_protein()
	var oxidant_before: float = sim.world.get_field("oxidant").total_amount()
	var neutral_before: float = sim.world.get_field("neutral_x").total_amount()
	var summary: Dictionary = ExtracellularCatalysisScript.step(sim.world, [e02], 1.0, config)
	var flux: float = float(summary["fluxes"]["E02"])
	var oxidant_consumed: float = oxidant_before - sim.world.get_field("oxidant").total_amount()
	var neutral_created: float = sim.world.get_field("neutral_x").total_amount() - neutral_before
	_assert_true(flux > 0.0 and flux <= oxidant_before, "E02 flux is positive but bounded by extracellular oxidant inventory")
	_assert_close(oxidant_consumed, flux, 1e-10, "E02 consumes exactly one tracked oxidant unit per unit flux")
	_assert_close(neutral_created, flux, 1e-10, "E02 creates exactly one neutral X unit per consumed oxidant unit")
	_assert_close(oxidant_consumed, neutral_created, 1e-10, "E02 conserves the tracked one-to-one detox species ledger")
	_assert_close(sim.world.total_extracellular_protein(), enzyme_before, 1e-12, "secreted DACE acts catalytically and is not consumed by E02")
	_assert_true(sim.world.get_field("oxidant").minimum_value() >= -1e-12 and sim.world.get_field("neutral_x").minimum_value() >= -1e-12, "E02 leaves extracellular fields nonnegative")

func _test_detox_producer_reduces_neighbor_ros_import_at_explicit_cost() -> void:
	var detox: Dictionary = _shared_site_exposure(true)
	var control: Dictionary = _shared_site_exposure(false)
	_assert_true(float(detox["secreted"]) > 0.0, "detox producer secretes a physically explicit DACE cohort")
	_assert_true(float(detox["producer_atp_spent"]) > 0.0, "detox producer pays explicit ATP for protein secretion")
	_assert_true(float(detox["producer_protein_lost"]) > 0.0, "detox producer loses the secreted molecules from its intracellular proteome")
	_assert_true(float(detox["e02_flux"]) > 0.0, "secreted DACE removes extracellular oxidant through E02")
	_assert_true(float(detox["local_oxidant_before_import"]) < float(control["local_oxidant_before_import"]), "detoxification lowers the chemical concentration encountered by the neighboring cell")
	_assert_true(float(detox["neighbor_ros_import"]) < float(control["neighbor_ros_import"]), "neighbor imports less ROS only because local extracellular oxidant was chemically reduced")

func _shared_site_exposure(with_detox: bool) -> Dictionary:
	var config = _config()
	config.secondary_extracellular_diffusion["ROS"] = 0.0
	var sim = SimulationEngineScript.new(config)
	var producer = sim.seed_ancestor(Vector2(6.0, 5.51))
	var neighbor = sim.seed_ancestor(Vector2(6.0, 6.49))
	_prepare_vulnerable_neighbor(neighbor)
	producer.set_pool("ATP", 8.0)
	producer.set_pool("ADP", 2.0)
	if with_detox:
		_install_exact_protein(producer, 5, DETOX_PROTEIN_SIGNATURE, 160.0)

	var protein_before: float = producer.total_protein()
	var atp_before: float = producer.pool("ATP")
	sim.world.release("oxidant", Vector2(6.0, 6.0), 2.0)
	var secretion: Dictionary = sim._secrete_extracellular_proteins(1.0)
	var e02 = ExtracellularReactionCatalogScript.by_id(sim.extracellular_reactions, "E02")
	var catalysis: Dictionary = ExtracellularCatalysisScript.step(sim.world, [e02], 1.0, config)
	var local_oxidant: float = sim.world.sample("oxidant", neighbor.position)
	var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var imported: float = maxf(0.0, float(transport["by_cell"][neighbor.id]["exchange"]["ROS"]))
	var producer_summary: Dictionary = secretion["by_cell"].get(producer.id, {"total_secreted": 0.0, "atp_spent": 0.0})
	return {
		"secreted": float(producer_summary.get("total_secreted", 0.0)),
		"producer_atp_spent": atp_before - producer.pool("ATP"),
		"producer_protein_lost": protein_before - producer.total_protein(),
		"e02_flux": float(catalysis["fluxes"]["E02"]),
		"local_oxidant_before_import": local_oxidant,
		"neighbor_ros_import": imported,
		"neighbor": neighbor,
		"sim": sim
	}

func _test_reduced_ros_import_reduces_neighbor_damage() -> void:
	var detox: Dictionary = _shared_site_exposure(true)
	var control: Dictionary = _shared_site_exposure(false)
	var detox_neighbor = detox["neighbor"]
	var control_neighbor = control["neighbor"]
	var config = _config()
	detox_neighbor.damage = 0.0
	control_neighbor.damage = 0.0
	detox_neighbor._update_damage_and_repair(1.0, config)
	control_neighbor._update_damage_and_repair(1.0, config)
	_assert_true(detox_neighbor.damage < control_neighbor.damage, "lower neighbor ROS burden produces less ordinary ROS-driven cellular damage")
	_assert_close(detox_neighbor.damage, detox_neighbor.pool("ROS"), 1e-10, "controlled no-decay/no-repair assay maps imported ROS into damage without a neighbor-benefit rule")
	_assert_close(control_neighbor.damage, control_neighbor.pool("ROS"), 1e-10, "control damage follows the same ordinary intracellular damage law")

func _test_spatial_separation_weakens_detoxification_benefit_through_field_gradient() -> void:
	var result: Dictionary = _spatial_detox_exposure()
	_assert_true(float(result["near_local_oxidant"]) < float(result["far_local_oxidant"]), "local detoxification creates a spatial oxidant gradient through ordinary diffusion")
	_assert_true(float(result["near_import"]) < float(result["far_import"]), "near neighbor receives greater detoxification benefit than distant neighbor")
	_assert_true(float(result["far_local_oxidant"]) <= 1.0 + 1e-12, "distance effect does not create oxidant or use a distance penalty")

func _spatial_detox_exposure() -> Dictionary:
	var config = _config()
	config.secondary_extracellular_diffusion["ROS"] = 1.0
	var sim = SimulationEngineScript.new(config)
	var producer = sim.seed_ancestor(Vector2(4.0, 6.0))
	var near_neighbor = sim.seed_ancestor(Vector2(5.0, 6.0))
	var far_neighbor = sim.seed_ancestor(Vector2(12.0, 6.0))
	_install_exact_protein(producer, 5, DETOX_PROTEIN_SIGNATURE, 160.0)
	_prepare_vulnerable_neighbor(near_neighbor)
	_prepare_vulnerable_neighbor(far_neighbor)
	producer.set_pool("ATP", 8.0)

	var oxidant_field = sim.world.get_field("oxidant")
	for y in range(oxidant_field.height):
		for x in range(oxidant_field.width):
			oxidant_field.set_value(x, y, 1.0)

	sim._secrete_extracellular_proteins(1.0)
	var e02 = ExtracellularReactionCatalogScript.by_id(sim.extracellular_reactions, "E02")
	ExtracellularCatalysisScript.step(sim.world, [e02], 1.0, config)
	for _i in range(10):
		sim.world.diffuse(config.tick_dt_min)
	var near_local: float = sim.world.sample("oxidant", near_neighbor.position)
	var far_local: float = sim.world.sample("oxidant", far_neighbor.position)
	var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	return {
		"near_local_oxidant": near_local,
		"far_local_oxidant": far_local,
		"near_import": maxf(0.0, float(transport["by_cell"][near_neighbor.id]["exchange"]["ROS"])),
		"far_import": maxf(0.0, float(transport["by_cell"][far_neighbor.id]["exchange"]["ROS"]))
	}

func _test_same_state_replays_detoxification_exactly() -> void:
	var first: Dictionary = _replay_chain(992731)
	var second: Dictionary = _replay_chain(992731)
	_assert_true(first["secretion"] == second["secretion"], "same state reproduces exact detox-protein secretion ledger")
	_assert_true(first["catalysis"] == second["catalysis"], "same state reproduces exact extracellular detox catalytic ledger")
	_assert_true(first["transport"] == second["transport"], "same state reproduces exact neighbor ROS transport ledger")
	_assert_close(float(first["checksum"]), float(second["checksum"]), 1e-12, "same state and seed reproduce exact post-detoxification checksum")

func _replay_chain(seed: int) -> Dictionary:
	var config = _config()
	config.seed = seed
	config.secondary_extracellular_diffusion["ROS"] = 0.0
	var sim = SimulationEngineScript.new(config)
	var producer = sim.seed_ancestor(Vector2(6.0, 5.51))
	var neighbor = sim.seed_ancestor(Vector2(6.0, 6.49))
	_install_exact_protein(producer, 5, DETOX_PROTEIN_SIGNATURE, 160.0)
	_prepare_vulnerable_neighbor(neighbor)
	producer.set_pool("ATP", 8.0)
	sim.world.release("oxidant", Vector2(6.0, 6.0), 2.0)
	var secretion: Dictionary = sim._secrete_extracellular_proteins(1.0)
	var e02 = ExtracellularReactionCatalogScript.by_id(sim.extracellular_reactions, "E02")
	var catalysis: Dictionary = ExtracellularCatalysisScript.step(sim.world, [e02], 1.0, config)
	var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	return {
		"secretion": secretion,
		"catalysis": catalysis,
		"transport": transport,
		"checksum": sim.checksum()
	}

func _prepare_vulnerable_neighbor(cell) -> void:
	_install_exact_protein(cell, 11, MembraneTransportScript.target_signature("ROS"), 160.0)
	cell.set_pool("ROS", 0.0)
	cell.damage = 0.0

func _install_exact_protein(cell, locus_id: int, signature: int, amount: float) -> void:
	var gene = cell.genome.get_gene_by_locus(locus_id)
	assert(gene != null)
	gene.protein_signature = signature & 0xFFFF
	cell.expression_state[locus_id]["protein"] = {signature & 0xFFFF: maxf(0.0, amount)}

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
