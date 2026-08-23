extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const ExtracellularCatalysisScript = preload("res://src/chemistry/extracellular_catalysis.gd")
const ExtracellularReactionCatalogScript = preload("res://src/chemistry/extracellular_reaction_catalog.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

const SECRETED_TEST_SIGNATURE: int = 0xD136

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_ancestor_has_no_protein_secretion()
	_test_secretion_and_extracellular_catalysis_are_mutationally_accessible()
	_test_secretion_moves_real_protein_and_spends_atp()
	_test_extracellular_reaction_is_bounded_and_materially_balanced()
	_test_nonproducer_exploits_public_product_in_shared_site()
	_test_spatial_local_capture_prevents_universal_nonproducer_advantage()
	_test_same_state_replays_public_good_history_exactly()

	if failures == 0:
		print("PASS: %d M7 extracellular-public-good tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M7 extracellular-public-good tests failed" % [failures, tests_run])
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

func _test_ancestor_has_no_protein_secretion() -> void:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	var intracellular_before: float = cell.total_protein()
	var summary: Dictionary = sim._secrete_extracellular_proteins(1.0)
	_assert_close(float(summary["total_secreted"]), 0.0, 1e-12, "ancestral realized proteome has no secretion-signal cohort")
	_assert_close(sim.world.total_extracellular_protein(), 0.0, 1e-12, "ancestral cell creates no extracellular protein field")
	_assert_close(cell.total_protein(), intracellular_before, 1e-12, "absence of secretion signal leaves ancestral proteome untouched")

func _test_secretion_and_extracellular_catalysis_are_mutationally_accessible() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var locus12 = ancestor.get_gene_by_locus(12)
	_assert_true(int(locus12.protein_signature) == 0xC136, "controlled precursor starts from canonical ancestral locus-12 sequence")
	_assert_true(CatalyticLandscapeScript.hamming_distance(int(locus12.protein_signature), SECRETED_TEST_SIGNATURE) == 1, "one ordinary coding bit flip creates the Dxxx secretion motif")
	_assert_true(ExtracellularCatalysisScript.has_secretion_signal(SECRETED_TEST_SIGNATURE), "one-bit mutant protein is recognized by generic secretion motif")

	var e01 = ExtracellularReactionCatalogScript.by_id(ExtracellularReactionCatalogScript.create_m7_candidate(), "E01")
	_assert_true(CatalyticLandscapeScript.hamming_distance(SECRETED_TEST_SIGNATURE, int(e01.signature)) == 4, "same secreted one-bit mutant lies inside ordinary extracellular catalytic radius")
	var ancestral_dormant: bool = true
	for gene in ancestor.genes:
		if CatalyticLandscapeScript.hamming_distance(int(gene.protein_signature), int(e01.signature)) <= CatalyticLandscapeScript.ACTIVE_MAX_DISTANCE:
			ancestral_dormant = false
			break
	_assert_true(ancestral_dormant, "E01 is outside catalytic radius of every ancestral protein")

func _test_secretion_moves_real_protein_and_spends_atp() -> void:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	_install_exact_protein(cell, 12, SECRETED_TEST_SIGNATURE, 160.0)
	cell.set_pool("ATP", 8.0)
	cell.set_pool("ADP", 2.0)
	var protein_before: float = cell.total_protein()
	var adenylate_before: float = cell.total_adenylate()
	var atp_before: float = cell.pool("ATP")
	var summary: Dictionary = sim._secrete_extracellular_proteins(1.0)
	var secreted: float = float(summary["by_cell"][cell.id]["total_secreted"])
	var spent: float = float(summary["by_cell"][cell.id]["atp_spent"])
	_assert_true(secreted > 0.0, "secretion-compatible realized protein leaves the cell")
	_assert_close(protein_before - cell.total_protein(), secreted, 1e-10, "secreted amount is removed from authoritative intracellular protein cohorts")
	_assert_close(sim.world.total_extracellular_protein(), secreted, 1e-10, "same physical protein amount appears in extracellular sequence field")
	_assert_close(spent, secreted * config.extracellular_protein_secretion_atp_cost_per_unit, 1e-12, "protein secretion ATP debit equals actual secreted molecules times unit cost")
	_assert_close(atp_before - cell.pool("ATP"), spent, 1e-12, "secretion directly reduces producer ATP")
	_assert_close(cell.total_adenylate(), adenylate_before, 1e-12, "secretion converts ATP to ADP without destroying adenylate")

func _test_extracellular_reaction_is_bounded_and_materially_balanced() -> void:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var site := Vector2(6.0, 6.0)
	sim.world.release_protein(SECRETED_TEST_SIGNATURE, site, 160.0, config.extracellular_protein_diffusion)
	sim.world.release("lipids", site, 1.0)
	var enzyme_before: float = sim.world.total_extracellular_protein()
	var lip_before: float = sim.world.get_field("lipids").total_amount()
	var c2_before: float = sim.world.get_field("carbon_c2").total_amount()
	var summary: Dictionary = ExtracellularCatalysisScript.step(sim.world, sim.extracellular_reactions, 1.0, config)
	var flux: float = float(summary["fluxes"]["E01"])
	var lip_consumed: float = lip_before - sim.world.get_field("lipids").total_amount()
	var c2_created: float = sim.world.get_field("carbon_c2").total_amount() - c2_before
	_assert_true(flux > 0.0 and flux <= lip_before, "extracellular E01 flux is positive but bounded by available substrate")
	_assert_close(lip_consumed, flux, 1e-10, "E01 consumes one extracellular LIP per unit flux")
	_assert_close(c2_created, 2.0 * flux, 1e-10, "E01 creates exactly two C2 per consumed C4 lipid precursor")
	_assert_close(4.0 * lip_consumed, 2.0 * c2_created, 1e-10, "extracellular lipid hydrolysis conserves modeled structural carbon")
	_assert_close(sim.world.total_extracellular_protein(), enzyme_before, 1e-12, "extracellular protein acts catalytically rather than being consumed as substrate")
	_assert_true(sim.world.get_field("lipids").minimum_value() >= -1e-12 and sim.world.get_field("carbon_c2").minimum_value() >= -1e-12, "extracellular reaction leaves every chemical field nonnegative")

func _test_nonproducer_exploits_public_product_in_shared_site() -> void:
	var config = _config()
	config.extracellular_protein_diffusion = 0.0
	config.secondary_extracellular_diffusion["C2"] = 0.0
	var sim = SimulationEngineScript.new(config)
	# Distinct non-overlapping disks round to one lattice site. Both consumers
	# therefore see the same public-product concentration and are allocated fairly.
	var producer = sim.seed_ancestor(Vector2(6.0, 5.51))
	var nonproducer = sim.seed_ancestor(Vector2(6.0, 6.49))
	_prepare_public_consumer(producer, config)
	_prepare_public_consumer(nonproducer, config)
	_install_exact_protein(producer, 12, SECRETED_TEST_SIGNATURE, 160.0)
	producer.set_pool("ATP", 8.0)
	producer.set_pool("ADP", 2.0)
	nonproducer.set_pool("ATP", 8.0)
	nonproducer.set_pool("ADP", 2.0)
	sim.world.release("lipids", Vector2(6.0, 6.0), 2.0)

	var producer_atp_before: float = producer.pool("ATP")
	var nonproducer_atp_before: float = nonproducer.pool("ATP")
	var secretion: Dictionary = sim._secrete_extracellular_proteins(1.0)
	var producer_secreted: float = float(secretion["by_cell"][producer.id]["total_secreted"])
	var nonproducer_secreted: float = float(secretion["by_cell"][nonproducer.id]["total_secreted"])
	_assert_true(producer_secreted > 0.0, "producer pays to place catalytic protein into shared extracellular space")
	_assert_close(nonproducer_secreted, 0.0, 1e-12, "nonproducer contributes no catalytic protein")

	var catalysis: Dictionary = ExtracellularCatalysisScript.step(sim.world, sim.extracellular_reactions, 1.0, config)
	_assert_true(float(catalysis["fluxes"]["E01"]) > 0.0 and sim.world.get_field("carbon_c2").total_amount() > 0.0, "producer enzyme converts shared LIP into public C2")

	var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var producer_import: float = float(transport["by_cell"][producer.id]["exchange"]["C2"])
	var nonproducer_import: float = float(transport["by_cell"][nonproducer.id]["exchange"]["C2"])
	_assert_true(producer_import > 0.0 and nonproducer_import > 0.0, "both cells can capture product using only matching generic C2 transport machinery")
	_assert_close(producer_import, nonproducer_import, 1e-12, "equal colocated transporters share public product without producer ownership bonus")
	_assert_true((producer_atp_before - producer.pool("ATP")) > (nonproducer_atp_before - nonproducer.pool("ATP")), "well-mixed nonproducer retains an energetic advantage by avoiding secretion cost while capturing equal product")

func _test_spatial_local_capture_prevents_universal_nonproducer_advantage() -> void:
	var result: Dictionary = _run_spatial_public_good(Vector2(12.0, 6.0))
	_assert_true(float(result["producer_import"]) > 0.0, "producer can recapture locally generated C2 after ordinary diffusion")
	_assert_true(float(result["producer_import"]) > float(result["nonproducer_import"]) + 1e-8, "distant nonproducer captures less public product than local producer")
	_assert_true(float(result["producer_local_c2"]) > float(result["nonproducer_local_c2"]), "spatial capture difference is explained by local C2 concentration rather than an ownership or distance rule")

func _run_spatial_public_good(nonproducer_position: Vector2) -> Dictionary:
	var config = _config()
	config.extracellular_protein_diffusion = 0.0
	config.secondary_extracellular_diffusion["C2"] = 0.65
	var sim = SimulationEngineScript.new(config)
	var producer = sim.seed_ancestor(Vector2(4.0, 6.0))
	var nonproducer = sim.seed_ancestor(nonproducer_position)
	_prepare_public_consumer(producer, config)
	_prepare_public_consumer(nonproducer, config)
	_install_exact_protein(producer, 12, SECRETED_TEST_SIGNATURE, 160.0)
	producer.set_pool("ATP", 8.0)
	nonproducer.set_pool("ATP", 8.0)
	sim.world.release("lipids", producer.position, 2.0)
	sim._secrete_extracellular_proteins(1.0)
	ExtracellularCatalysisScript.step(sim.world, sim.extracellular_reactions, 1.0, config)
	for _i in range(12):
		sim.world.diffuse(config.tick_dt_min)
	var producer_local: float = sim.world.sample("carbon_c2", producer.position)
	var nonproducer_local: float = sim.world.sample("carbon_c2", nonproducer.position)
	var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	return {
		"producer_import": maxf(0.0, float(transport["by_cell"][producer.id]["exchange"]["C2"])),
		"nonproducer_import": maxf(0.0, float(transport["by_cell"][nonproducer.id]["exchange"]["C2"])),
		"producer_local_c2": producer_local,
		"nonproducer_local_c2": nonproducer_local
	}

func _test_same_state_replays_public_good_history_exactly() -> void:
	var first: Dictionary = _run_replay_chain(884217)
	var second: Dictionary = _run_replay_chain(884217)
	_assert_true(first["secretion"] == second["secretion"], "same state reproduces exact protein-secretion ledger")
	_assert_true(first["catalysis"] == second["catalysis"], "same state reproduces exact extracellular catalytic flux ledger")
	_assert_true(first["transport"] == second["transport"], "same state reproduces exact public-product uptake ledger")
	_assert_close(float(first["checksum"]), float(second["checksum"]), 1e-12, "same state and seed reproduce exact extracellular-protein/public-good checksum")

func _run_replay_chain(seed: int) -> Dictionary:
	var config = _config()
	config.seed = seed
	config.extracellular_protein_diffusion = 0.0
	config.secondary_extracellular_diffusion["C2"] = 0.0
	var sim = SimulationEngineScript.new(config)
	var producer = sim.seed_ancestor(Vector2(6.0, 5.51))
	var nonproducer = sim.seed_ancestor(Vector2(6.0, 6.49))
	_prepare_public_consumer(producer, config)
	_prepare_public_consumer(nonproducer, config)
	_install_exact_protein(producer, 12, SECRETED_TEST_SIGNATURE, 160.0)
	producer.set_pool("ATP", 8.0)
	nonproducer.set_pool("ATP", 8.0)
	sim.world.release("lipids", Vector2(6.0, 6.0), 2.0)
	var secretion: Dictionary = sim._secrete_extracellular_proteins(1.0)
	var catalysis: Dictionary = ExtracellularCatalysisScript.step(sim.world, sim.extracellular_reactions, 1.0, config)
	var transport: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	return {
		"secretion": secretion,
		"catalysis": catalysis,
		"transport": transport,
		"checksum": sim.checksum()
	}

func _prepare_public_consumer(cell, config) -> void:
	_install_exact_protein(cell, 11, MembraneTransportScript.target_signature("C2"), 160.0)
	cell.set_pool("C2", 0.0)

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
