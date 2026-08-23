extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const MembraneTransportScript = preload("res://src/transport/membrane_transport.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_ancestor_has_no_secondary_transport()
	_test_every_secondary_transport_is_one_coding_step_reachable()
	_test_realized_protein_not_dna_alone_controls_transport()
	_test_same_generic_mechanism_imports_and_exports_by_gradient()
	_test_transport_conserves_molecule_and_adenylate()
	_test_scarce_external_substrate_is_shared_proportionally()
	_test_atp_scarcity_scales_all_secondary_exchange()
	_test_same_state_replays_secondary_exchange_exactly()

	if failures == 0:
		print("PASS: %d M7 membrane-transport tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M7 membrane-transport tests failed" % [failures, tests_run])
		quit(1)

func _config():
	var config = SimConfigScript.new()
	config.world_width = 12
	config.world_height = 12
	config.initial_glucose = 0.0
	config.initial_oxygen = 0.0
	config.initial_nitrogen = 0.0
	config.initial_phosphorus = 0.0
	config.mutation_enabled = false
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		config.secondary_extracellular_initial[metabolite_id] = 0.0
	return config

func _test_ancestor_has_no_secondary_transport() -> void:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	var all_zero: bool = true
	for metabolite_id in config.SECONDARY_EXTRACELLULAR_IDS:
		var activity: float = MembraneTransportScript.proteome_activity(cell.expression_state, metabolite_id, config)
		if activity > 1e-15:
			all_zero = false
			break
	_assert_true(all_zero, "M6 ancestor proteome has zero secondary membrane-transport activity")

	var w1_field: String = MetaboliteCatalogScript.extracellular_field("W1")
	sim.world.release(w1_field, cell.position, 10.0)
	var before_external: float = sim.world.get_field(w1_field).total_amount()
	var before_internal: float = cell.pool("W1")
	var summary: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	_assert_close(cell.pool("W1"), before_internal, 1e-12, "secondary field presence does not give ancestor W1 uptake")
	_assert_close(sim.world.get_field(w1_field).total_amount(), before_external, 1e-12, "ancestor leaves inaccessible W1 field untouched")
	_assert_close(float(summary["total_moved"]), 0.0, 1e-12, "ancestor performs no secondary exchange")

func _test_every_secondary_transport_is_one_coding_step_reachable() -> void:
	var genome = GenomeScript.create_ancestor()
	for metabolite_id in SimConfigScript.SECONDARY_EXTRACELLULAR_IDS:
		var reach: Dictionary = _nearest_ancestor_signature(genome, metabolite_id)
		_assert_true(int(reach["distance"]) == 5, "%s transport target is dormant exactly one coding bit beyond active radius" % metabolite_id)
		var one_step: int = _one_bit_toward(int(reach["signature"]), MembraneTransportScript.target_signature(metabolite_id))
		_assert_true(MembraneTransportScript.hamming_distance(one_step, int(reach["signature"])) == 1, "%s transport capability is reached by one ordinary coding bit flip" % metabolite_id)
		_assert_true(MembraneTransportScript.affinity(one_step, metabolite_id) > 0.0, "%s one-step mutant protein enters active transport radius" % metabolite_id)

func _test_realized_protein_not_dna_alone_controls_transport() -> void:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	var reach: Dictionary = _nearest_ancestor_signature(cell.genome, "W1")
	var gene = cell.genome.get_gene_by_locus(int(reach["locus_id"]))
	var mutant_signature: int = _one_bit_toward(int(gene.protein_signature), MembraneTransportScript.target_signature("W1"))
	gene.protein_signature = mutant_signature
	_assert_close(MembraneTransportScript.proteome_activity(cell.expression_state, "W1", config), 0.0, 1e-12, "coding mutation alone does not rewrite inherited proteome into a transporter")

	var locus_state: Dictionary = cell.expression_state[int(gene.locus_id)]
	var inherited_amount: float = 0.0
	for amount in locus_state["protein"].values():
		inherited_amount += float(amount)
	locus_state["protein"] = {mutant_signature: inherited_amount}
	_assert_true(MembraneTransportScript.proteome_activity(cell.expression_state, "W1", config) > 0.0, "transport appears only when compatible protein physically exists")

func _test_same_generic_mechanism_imports_and_exports_by_gradient() -> void:
	var import_config = _config()
	var import_sim = SimulationEngineScript.new(import_config)
	var importer = import_sim.seed_ancestor(Vector2(6.0, 6.0))
	_install_exact_transporter(importer, "W1")
	importer.set_pool("W1", 0.0)
	var field_name: String = MetaboliteCatalogScript.extracellular_field("W1")
	import_sim.world.release(field_name, importer.position, 5.0)
	var import_summary: Dictionary = import_sim._allocate_secondary_membrane_transport(import_config.tick_dt_min)
	var imported: float = float(import_summary["by_cell"][importer.id]["exchange"]["W1"])
	_assert_true(imported > 0.0 and importer.pool("W1") > 0.0, "compatible protein imports W1 when extracellular concentration is higher")

	var export_config = _config()
	var export_sim = SimulationEngineScript.new(export_config)
	var exporter = export_sim.seed_ancestor(Vector2(6.0, 6.0))
	_install_exact_transporter(exporter, "W1")
	exporter.set_pool("W1", 5.0)
	var export_summary: Dictionary = export_sim._allocate_secondary_membrane_transport(export_config.tick_dt_min)
	var exported_signed: float = float(export_summary["by_cell"][exporter.id]["exchange"]["W1"])
	_assert_true(exported_signed < 0.0 and export_sim.world.get_field(field_name).total_amount() > 0.0, "same compatible protein exports W1 when intracellular concentration is higher")
	_assert_true(not import_summary["by_cell"][importer.id].has("mode") and not export_summary["by_cell"][exporter.id].has("mode"), "transport direction requires no import/export behavioral state")

func _test_transport_conserves_molecule_and_adenylate() -> void:
	var config = _config()
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	_install_exact_transporter(cell, "W1")
	cell.set_pool("W1", 0.75)
	var field_name: String = MetaboliteCatalogScript.extracellular_field("W1")
	sim.world.release(field_name, cell.position, 3.0)
	var molecule_before: float = cell.pool("W1") + sim.world.get_field(field_name).total_amount()
	var adenylate_before: float = cell.total_adenylate()
	var atp_before: float = cell.pool("ATP")
	var summary: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var moved: float = float(summary["by_cell"][cell.id]["moved"])
	var spent: float = float(summary["by_cell"][cell.id]["atp_spent"])
	_assert_true(moved > 0.0 and spent > 0.0 and cell.pool("ATP") < atp_before, "secondary exchange moves material and spends explicit ATP")
	_assert_close(cell.pool("W1") + sim.world.get_field(field_name).total_amount(), molecule_before, 1e-10, "membrane exchange conserves W1 across cell plus world")
	_assert_close(cell.total_adenylate(), adenylate_before, 1e-10, "transport converts ATP to ADP without creating/destroying adenylate")
	_assert_close(spent, moved * config.secondary_transport_atp_cost_per_unit, 1e-12, "transport ATP debit equals actual moved units times configured unit cost")

func _test_scarce_external_substrate_is_shared_proportionally() -> void:
	var config = _config()
	config.secondary_transport_vmax_per_reference_protein = 100.0
	var sim = SimulationEngineScript.new(config)
	var first = sim.seed_ancestor(Vector2(6.05, 6.05))
	var second = sim.seed_ancestor(Vector2(6.15, 6.15))
	_install_exact_transporter(first, "W1")
	_install_exact_transporter(second, "W1")
	first.set_pool("W1", 0.0)
	second.set_pool("W1", 0.0)
	var field_name: String = MetaboliteCatalogScript.extracellular_field("W1")
	sim.world.release(field_name, Vector2(6.0, 6.0), 0.01)
	var summary: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var first_import: float = float(summary["by_cell"][first.id]["exchange"]["W1"])
	var second_import: float = float(summary["by_cell"][second.id]["exchange"]["W1"])
	_assert_true(first_import > 0.0 and second_import > 0.0, "both compatible cells receive scarce shared substrate")
	_assert_close(first_import, second_import, 1e-12, "equal simultaneous importers receive equal proportional allocation independent of iteration priority")
	_assert_close(first_import + second_import, 0.01, 1e-10, "scarce site inventory is exhausted but never overdrawn")
	_assert_close(sim.world.get_field(field_name).total_amount(), 0.0, 1e-10, "proportional competition cannot make extracellular field negative")

func _test_atp_scarcity_scales_all_secondary_exchange() -> void:
	var config = _config()
	config.secondary_transport_vmax_per_reference_protein = 100.0
	config.secondary_transport_atp_cost_per_unit = 1.0
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(6.0, 6.0))
	_install_exact_transporter(cell, "W1")
	_install_exact_transporter(cell, "W2")
	cell.set_pool("W1", 0.0)
	cell.set_pool("W2", 0.0)
	cell.set_pool("ATP", 0.001)
	var adenylate_before: float = cell.total_adenylate()
	sim.world.release(MetaboliteCatalogScript.extracellular_field("W1"), cell.position, 10.0)
	sim.world.release(MetaboliteCatalogScript.extracellular_field("W2"), cell.position, 10.0)
	var summary: Dictionary = sim._allocate_secondary_membrane_transport(config.tick_dt_min)
	var cell_summary: Dictionary = summary["by_cell"][cell.id]
	_assert_true(float(cell_summary["energy_scale"]) < 1.0, "ATP scarcity proportionally limits simultaneous secondary transport proposals")
	_assert_close(float(cell_summary["moved"]), 0.001, 1e-10, "ATP-limited cell cannot move more secondary material than it can energetically fund")
	_assert_close(cell.pool("ATP"), 0.0, 1e-10, "ATP-limited transport spends available ATP but never goes negative")
	_assert_close(cell.total_adenylate(), adenylate_before, 1e-10, "ATP-limited transport still conserves ATP plus ADP")
	var w1: float = float(cell_summary["exchange"]["W1"])
	var w2: float = float(cell_summary["exchange"]["W2"])
	_assert_true(w1 > 0.0 and w2 > 0.0, "ATP allocation does not privilege one simultaneously transportable molecule")

func _test_same_state_replays_secondary_exchange_exactly() -> void:
	var first_config = _config()
	var second_config = _config()
	first_config.seed = 994411
	second_config.seed = 994411
	var first = SimulationEngineScript.new(first_config)
	var second = SimulationEngineScript.new(second_config)
	var first_cell = first.seed_ancestor(Vector2(5.0, 5.0))
	var second_cell = second.seed_ancestor(Vector2(5.0, 5.0))
	for metabolite_id in ["W1", "W2", "C2"]:
		_install_exact_transporter(first_cell, metabolite_id)
		_install_exact_transporter(second_cell, metabolite_id)
		first.world.release(MetaboliteCatalogScript.extracellular_field(metabolite_id), first_cell.position, 2.0)
		second.world.release(MetaboliteCatalogScript.extracellular_field(metabolite_id), second_cell.position, 2.0)
	var first_summary: Dictionary = first._allocate_secondary_membrane_transport(first_config.tick_dt_min)
	var second_summary: Dictionary = second._allocate_secondary_membrane_transport(second_config.tick_dt_min)
	_assert_true(first_summary == second_summary, "same state reproduces exact secondary transport ledger")
	_assert_close(first.checksum(), second.checksum(), 1e-12, "same state and seed reproduce exact post-transport world/cell checksum")

func _nearest_ancestor_signature(genome, metabolite_id: String) -> Dictionary:
	var target: int = MembraneTransportScript.target_signature(metabolite_id)
	var best_distance: int = 99
	var best_locus: int = -1
	var best_signature: int = -1
	for gene in genome.genes:
		var distance: int = MembraneTransportScript.hamming_distance(int(gene.protein_signature), target)
		if distance < best_distance:
			best_distance = distance
			best_locus = int(gene.locus_id)
			best_signature = int(gene.protein_signature)
	return {"distance": best_distance, "locus_id": best_locus, "signature": best_signature}

func _one_bit_toward(source: int, target: int) -> int:
	var difference: int = (source ^ target) & 0xFFFF
	assert(difference != 0)
	for bit_index in range(16):
		var bit: int = 1 << bit_index
		if (difference & bit) != 0:
			return source ^ bit
	return source

func _install_exact_transporter(cell, metabolite_id: String) -> void:
	var reach: Dictionary = _nearest_ancestor_signature(cell.genome, metabolite_id)
	var locus_id: int = int(reach["locus_id"])
	var target: int = MembraneTransportScript.target_signature(metabolite_id)
	var gene = cell.genome.get_gene_by_locus(locus_id)
	gene.protein_signature = target
	var locus_state: Dictionary = cell.expression_state[locus_id]
	var amount: float = 0.0
	for cohort_amount in locus_state["protein"].values():
		amount += float(cohort_amount)
	locus_state["protein"] = {target: amount}

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
