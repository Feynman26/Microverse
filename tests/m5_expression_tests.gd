extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ExpressionSolverScript = preload("res://src/expression/expression_solver.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_initial_proteome_matches_m4_scale()
	_test_zero_proteome_means_zero_flux()
	_test_expression_has_energy_and_material_cost()
	_test_expression_turnover_conserves_structural_material()
	_test_generic_activation_and_repression()
	_test_generic_metabolite_sensing_changes_regulation()
	_test_expression_noise_is_seed_reproducible()
	_test_expression_partition_conserves_state()
	_test_regulatory_motif_mutation_changes_network_genotype()
	_test_coding_mutation_preserves_old_protein_identity_until_turnover()
	_test_full_simulation_expression_is_reproducible()
	if failures == 0:
		print("PASS: %d M5 expression/regulation tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5 tests failed" % [failures, tests_run])
		quit(1)

func _test_initial_proteome_matches_m4_scale() -> void:
	var config = SimConfigScript.new(); config.regulation_enabled = false
	var genome = GenomeScript.create_ancestor()
	var state = ExpressionSolverScript.initialize(genome, config)
	for gene in genome.genes:
		_assert_close(state.protein_for(int(gene.locus_id)), float(gene.promoter_strength()), 1e-10, "initial protein abundance preserves M4 promoter scale")
		_assert_close(state.protein_for_signature(int(gene.locus_id), int(gene.protein_signature)), float(gene.promoter_strength()), 1e-10, "initial cohort carries encoded ancestral signature")

func _test_zero_proteome_means_zero_flux() -> void:
	var config = SimConfigScript.new()
	var genome = GenomeScript.create_ancestor()
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["G"] = 5.0; pools["O2"] = 5.0; pools["NH4"] = 5.0; pools["P"] = 5.0
	var proteins: Dictionary = {}
	for gene in genome.genes: proteins[int(gene.locus_id)] = 0.0
	var fluxes: Dictionary = MetabolicSolverScript.step(pools, genome, reactions, 1.0, 1.0, config, proteins)
	var total: float = 0.0
	for value in fluxes.values(): total += float(value)
	_assert_close(total, 0.0, 1e-12, "explicit zero proteome produces zero catalytic flux")

func _test_expression_has_energy_and_material_cost() -> void:
	var config = SimConfigScript.new(); config.regulation_enabled = false; config.expression_noise_fraction = 0.0
	var genome = GenomeScript.create_ancestor()
	var state = ExpressionSolverScript.initialize(genome, config)
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 20.0; pools["ADP"] = 0.0; pools["AA"] = 10.0; pools["NUC"] = 10.0
	var before_atp: float = float(pools["ATP"])
	var stats: Dictionary = ExpressionSolverScript.step(state, genome, pools, 0.1, DeterministicRngScript.new(11), config)
	_assert_true(float(stats["atp_cost"]) > 0.0, "gene expression consumes ATP")
	_assert_true(float(pools["ATP"]) < before_atp and float(pools["ADP"]) > 0.0, "expression converts ATP to ADP")
	_assert_true(float(stats["aa_committed"]) > 0.0, "translation commits amino precursor material to proteins")
	_assert_true(float(stats["nuc_committed"]) > 0.0, "transcription commits nucleotide precursor material to mRNA")

func _test_expression_turnover_conserves_structural_material() -> void:
	var config = SimConfigScript.new(); config.regulation_enabled = false; config.expression_noise_fraction = 0.0
	var genome = GenomeScript.create_ancestor()
	var state = ExpressionSolverScript.initialize(genome, config)
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 20.0; pools["AA"] = 10.0; pools["NUC"] = 10.0
	var before_met: Dictionary = MetabolicSolverScript.structural_totals(pools)
	var before_expr: Dictionary = ExpressionSolverScript.structural_totals(state, config)
	ExpressionSolverScript.step(state, genome, pools, 0.3, DeterministicRngScript.new(19), config)
	var after_met: Dictionary = MetabolicSolverScript.structural_totals(pools)
	var after_expr: Dictionary = ExpressionSolverScript.structural_totals(state, config)
	for element in ["C", "N", "P"]:
		_assert_close(float(before_met[element]) + float(before_expr[element]), float(after_met[element]) + float(after_expr[element]), 1e-9, "expression turnover conserves represented %s" % element)

func _test_generic_activation_and_repression() -> void:
	var config = SimConfigScript.new(); config.expression_noise_fraction = 0.0; config.allostery_enabled = false; config.regulatory_gain = 0.5
	var activator = GeneScript.new(1, 5000, 0x1234, 1, 0x0000)
	var repressor = GeneScript.new(2, 5000, 0x9234, 2, 0x0000)
	var target_a = GeneScript.new(3, 5000, 0x3333, 3, 0x1234)
	var target_r = GeneScript.new(4, 5000, 0x4444, 4, 0x9234)
	var genome = GenomeScript.new([activator, repressor, target_a, target_r])
	var state = ExpressionSolverScript.initialize(genome, config)
	state.set_single_protein_cohort(1, 0x1234, 2.0)
	state.set_single_protein_cohort(2, 0x9234, 2.0)
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 20.0; pools["AA"] = 20.0; pools["NUC"] = 20.0
	var factors: Dictionary = ExpressionSolverScript.step(state, genome, pools, 0.1, DeterministicRngScript.new(22), config)["regulation"]
	_assert_true(float(factors[3]) > 1.0, "matching low-high-bit protein generically activates target promoter")
	_assert_true(float(factors[4]) < 1.0, "matching high-bit protein generically represses target promoter")

func _test_generic_metabolite_sensing_changes_regulation() -> void:
	var config = SimConfigScript.new(); config.expression_noise_fraction = 0.0; config.regulatory_gain = 0.8; config.allosteric_gain = 1.0
	var regulator = GeneScript.new(1, 6000, 0x0000, 1, 0xBEEF)
	var target = GeneScript.new(2, 5000, 0x2222, 2, 0x0000)
	var genome = GenomeScript.new([regulator, target])
	var low_state = ExpressionSolverScript.initialize(genome, config)
	var high_state = low_state.deep_copy()
	low_state.set_single_protein_cohort(1, 0x0000, 1.0)
	high_state.set_single_protein_cohort(1, 0x0000, 1.0)
	var low_pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	var high_pools: Dictionary = low_pools.duplicate(true)
	for pools in [low_pools, high_pools]: pools["ATP"] = 20.0; pools["AA"] = 20.0; pools["NUC"] = 20.0
	low_pools["G"] = 0.0; high_pools["G"] = 5.0
	var low_stats: Dictionary = ExpressionSolverScript.step(low_state, genome, low_pools, 0.1, DeterministicRngScript.new(55), config)
	var high_stats: Dictionary = ExpressionSolverScript.step(high_state, genome, high_pools, 0.1, DeterministicRngScript.new(55), config)
	_assert_true(float(high_stats["regulation"][2]) > float(low_stats["regulation"][2]), "ordinary ligand binding makes the same regulatory circuit environment-conditioned")
	_assert_true(high_state.mrna_for(2) > low_state.mrna_for(2), "chemical sensing changes downstream molecular phenotype without a named response API")

func _test_expression_noise_is_seed_reproducible() -> void:
	var config = SimConfigScript.new(); config.regulation_enabled = false; config.expression_noise_fraction = 0.25
	var genome = GenomeScript.create_ancestor()
	var first = ExpressionSolverScript.initialize(genome, config); var second = ExpressionSolverScript.initialize(genome, config)
	var pools_a: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config); var pools_b: Dictionary = pools_a.duplicate(true)
	pools_a["ATP"] = 50.0; pools_b["ATP"] = 50.0; pools_a["AA"] = 50.0; pools_b["AA"] = 50.0; pools_a["NUC"] = 50.0; pools_b["NUC"] = 50.0
	var rng_a = DeterministicRngScript.new(773); var rng_b = DeterministicRngScript.new(773)
	for _i in range(20):
		ExpressionSolverScript.step(first, genome, pools_a, 0.1, rng_a, config)
		ExpressionSolverScript.step(second, genome, pools_b, 0.1, rng_b, config)
	_assert_close(first.checksum(), second.checksum(), 1e-12, "same seed reproduces expression-noise trajectory")

func _test_expression_partition_conserves_state() -> void:
	var config = SimConfigScript.new()
	var genome = GenomeScript.create_ancestor(); var state = ExpressionSolverScript.initialize(genome, config)
	var before_mrna: float = state.total_mrna(); var before_protein: float = state.total_protein()
	var parts: Array = ExpressionSolverScript.partition(state, 0.5, DeterministicRngScript.new(331), config)
	_assert_close(parts[0].total_mrna() + parts[1].total_mrna(), before_mrna, 1e-10, "division conserves total mRNA state")
	_assert_close(parts[0].total_protein() + parts[1].total_protein(), before_protein, 1e-10, "division conserves total protein state")
	_assert_true(absf(parts[0].total_protein() - parts[1].total_protein()) > 1e-9, "seeded partition noise can separate clonal daughter phenotypes")

func _test_regulatory_motif_mutation_changes_network_genotype() -> void:
	var config = SimConfigScript.new()
	config.promoter_mutation_rate_per_gene = 0.0; config.signature_mutation_rate_per_gene = 0.0; config.regulatory_signature_mutation_rate_per_gene = 1.0; config.neutral_marker_mutation_rate_per_gene = 0.0
	var sim = SimulationEngineScript.new(config)
	var parent = GenomeScript.create_ancestor(); var old_key: String = parent.canonical_key()
	var result: Dictionary = sim.mutation_engine.mutate_copy(parent, DeterministicRngScript.new(44), config)
	_assert_true(String(result["genome"].canonical_key()) != old_key, "regulatory motif mutation changes inherited network genotype")
	_assert_true(result["events"].size() == parent.gene_count(), "forced regulatory mutation emits one event per locus")

func _test_coding_mutation_preserves_old_protein_identity_until_turnover() -> void:
	var config = SimConfigScript.new(); config.regulation_enabled = false; config.expression_noise_fraction = 0.0
	var genome = GenomeScript.create_ancestor()
	var state = ExpressionSolverScript.initialize(genome, config)
	var gene = genome.get_gene_by_locus(1)
	var old_signature: int = int(gene.protein_signature)
	var new_signature: int = old_signature ^ 1
	var inherited_old: float = state.protein_for_signature(1, old_signature)
	gene.protein_signature = new_signature
	_assert_close(state.protein_for_signature(1, old_signature), inherited_old, 1e-12, "coding mutation does not rewrite inherited protein cohort")
	_assert_close(state.protein_for_signature(1, new_signature), 0.0, 1e-12, "new coding signature has no protein before translation")
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 20.0; pools["AA"] = 20.0; pools["NUC"] = 20.0
	ExpressionSolverScript.step(state, genome, pools, 0.5, DeterministicRngScript.new(99), config)
	_assert_true(state.protein_for_signature(1, old_signature) < inherited_old, "old protein cohort decays after mutation")
	_assert_true(state.protein_for_signature(1, new_signature) > 0.0, "translation gradually creates protein with mutated signature")

func _test_full_simulation_expression_is_reproducible() -> void:
	var a = SimConfigScript.new(); a.seed = 9182; a.max_cells = 12
	var b = SimConfigScript.new(); b.seed = 9182; b.max_cells = 12
	var first = SimulationEngineScript.new(a); var second = SimulationEngineScript.new(b)
	first.seed_ancestor(); second.seed_ancestor(); first.step(350); second.step(350)
	_assert_close(first.checksum(), second.checksum(), 1e-9, "same seed reproduces full M5 molecular state")
	_assert_true(first.event_log == second.event_log, "same seed reproduces M5 event history")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition: print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
