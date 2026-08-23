extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

var failures: int = 0
var tests_run: int = 0
var config
var reactions: Array

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	config = SimConfigScript.new()
	reactions = ReactionCatalogScript.create_m4_candidate()
	_test_equilibrium_initialization_matches_m4_abundance_scale()
	_test_same_seed_reproduces_expression_noise()
	_test_different_seeds_can_separate_clone_expression()
	_test_expression_conserves_structural_material_with_recycling()
	_test_expression_converts_atp_to_adp()
	_test_gene_array_order_does_not_change_expression_history()
	_test_protein_abundance_controls_catalytic_flux()
	_test_signature_mutation_does_not_rewrite_inherited_proteins()
	_test_new_allele_builds_new_protein_cohort_over_time()
	_test_stochastic_partition_conserves_molecules_and_separates_sisters()
	_test_poisson_sampler_matches_expected_mean()
	_test_resource_scarcity_scales_synthesis_without_negative_pools()

	if failures == 0:
		print("PASS: %d M5 expression-core tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5 expression-core tests failed" % [failures, tests_run])
		quit(1)

func _rich_expression_pools() -> Dictionary:
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 20.0
	pools["ADP"] = 10.0
	pools["NUC"] = 10.0
	pools["AA"] = 10.0
	return pools

func _combined_structural(pools: Dictionary, expression: Dictionary) -> Dictionary:
	var chemical: Dictionary = MetabolicSolverScript.structural_totals(pools)
	var stored: Dictionary = ExpressionSystemScript.structural_storage_totals(expression, config)
	return {
		"C": float(chemical["C"]) + float(stored["C"]),
		"N": float(chemical["N"]) + float(stored["N"]),
		"P": float(chemical["P"]) + float(stored["P"])
	}

func _test_equilibrium_initialization_matches_m4_abundance_scale() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var all_match: bool = true
	for gene in genome.genes:
		var protein: float = ExpressionSystemScript.current_gene_protein(expression, gene)
		var normalized: float = protein / float(config.expression_reference_protein_count)
		if absf(normalized - float(gene.promoter_strength())) > 1e-12:
			all_match = false
			break
	_assert_true(all_match, "M5 equilibrium initialization preserves the former M4 promoter-abundance scale")

func _test_same_seed_reproduces_expression_noise() -> void:
	var genome = GenomeScript.create_ancestor()
	var first: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var second: Dictionary = first.duplicate(true)
	var pools_a: Dictionary = _rich_expression_pools()
	var pools_b: Dictionary = pools_a.duplicate(true)
	var rng_a = DeterministicRngScript.new(192837)
	var rng_b = DeterministicRngScript.new(192837)
	for _i in range(150):
		ExpressionSystemScript.step(first, genome, pools_a, config.tick_dt_min, rng_a, config)
		ExpressionSystemScript.step(second, genome, pools_b, config.tick_dt_min, rng_b, config)
	_assert_close(ExpressionSystemScript.checksum(first), ExpressionSystemScript.checksum(second), 1e-12, "same seed reproduces exact expression-noise trajectory")
	_assert_true(pools_a == pools_b, "same expression seed reproduces exact synthesis resource history")

func _test_different_seeds_can_separate_clone_expression() -> void:
	var genome = GenomeScript.create_ancestor()
	var first: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var second: Dictionary = first.duplicate(true)
	var pools_a: Dictionary = _rich_expression_pools()
	var pools_b: Dictionary = pools_a.duplicate(true)
	var rng_a = DeterministicRngScript.new(101)
	var rng_b = DeterministicRngScript.new(202)
	for _i in range(80):
		ExpressionSystemScript.step(first, genome, pools_a, config.tick_dt_min, rng_a, config)
		ExpressionSystemScript.step(second, genome, pools_b, config.tick_dt_min, rng_b, config)
	_assert_true(absf(ExpressionSystemScript.checksum(first) - ExpressionSystemScript.checksum(second)) > 1e-8, "different seeds can create phenotypic expression divergence between identical genotypes")

func _test_expression_conserves_structural_material_with_recycling() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var pools: Dictionary = _rich_expression_pools()
	var before: Dictionary = _combined_structural(pools, expression)
	var rng = DeterministicRngScript.new(778899)
	for _i in range(100):
		ExpressionSystemScript.step(expression, genome, pools, config.tick_dt_min, rng, config)
	var after: Dictionary = _combined_structural(pools, expression)
	_assert_close(float(after["C"]), float(before["C"]), 1e-9, "expression transfers but does not create/destroy modeled structural carbon")
	_assert_close(float(after["N"]), float(before["N"]), 1e-9, "expression transfers but does not create/destroy modeled structural nitrogen")
	_assert_close(float(after["P"]), float(before["P"]), 1e-9, "expression transfers but does not create/destroy modeled structural phosphorus")

func _test_expression_converts_atp_to_adp() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var pools: Dictionary = _rich_expression_pools()
	var adenylate_before: float = float(pools["ATP"]) + float(pools["ADP"])
	var atp_before: float = float(pools["ATP"])
	var rng = DeterministicRngScript.new(31337)
	var spent: float = 0.0
	for _i in range(50):
		var summary: Dictionary = ExpressionSystemScript.step(expression, genome, pools, config.tick_dt_min, rng, config)
		spent += float(summary["atp_spent"])
	_assert_true(spent > 0.0 and float(pools["ATP"]) < atp_before, "transcription/translation spend explicit ATP")
	_assert_close(float(pools["ATP"]) + float(pools["ADP"]), adenylate_before, 1e-9, "expression spending conserves ATP+ADP currency")

func _test_gene_array_order_does_not_change_expression_history() -> void:
	var genome_a = GenomeScript.create_ancestor()
	var genome_b = genome_a.deep_copy()
	genome_b.genes.reverse()
	var state_a: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome_a, config)
	var state_b: Dictionary = state_a.duplicate(true)
	var pools_a: Dictionary = _rich_expression_pools()
	var pools_b: Dictionary = pools_a.duplicate(true)
	var rng_a = DeterministicRngScript.new(50505)
	var rng_b = DeterministicRngScript.new(50505)
	for _i in range(100):
		ExpressionSystemScript.step(state_a, genome_a, pools_a, config.tick_dt_min, rng_a, config)
		ExpressionSystemScript.step(state_b, genome_b, pools_b, config.tick_dt_min, rng_b, config)
	_assert_close(ExpressionSystemScript.checksum(state_a), ExpressionSystemScript.checksum(state_b), 1e-12, "reversing genome array order does not alter stochastic expression trajectory")
	_assert_true(pools_a == pools_b, "gene array order does not alter shared synthesis-resource allocation")

func _test_protein_abundance_controls_catalytic_flux() -> void:
	var genome = GenomeScript.create_ancestor()
	var low_expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var high_expression: Dictionary = low_expression.duplicate(true)
	for locus_id in high_expression.keys():
		for signature in high_expression[locus_id]["protein"].keys():
			high_expression[locus_id]["protein"][signature] = float(high_expression[locus_id]["protein"][signature]) * 2.0
	var low_pools: Dictionary = _rich_expression_pools()
	var high_pools: Dictionary = low_pools.duplicate(true)
	low_pools["G"] = 10.0
	high_pools["G"] = 10.0
	low_pools["NAD"] = 10.0
	high_pools["NAD"] = 10.0
	var low_flux: Dictionary = MetabolicSolverScript.step(low_pools, genome, low_expression, reactions, 0.25, 1.0, config)
	var high_flux: Dictionary = MetabolicSolverScript.step(high_pools, genome, high_expression, reactions, 0.25, 1.0, config)
	_assert_true(float(high_flux["R01"]) > float(low_flux["R01"]), "more matching protein produces more catalytic flux before substrate limitation dominates")

func _test_signature_mutation_does_not_rewrite_inherited_proteins() -> void:
	var ancestral_genome = GenomeScript.create_ancestor()
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(ancestral_genome, config)
	var mutant = ancestral_genome.deep_copy()
	var gene = mutant.get_gene_by_locus(1)
	var old_signature: int = int(gene.protein_signature)
	gene.protein_signature = old_signature ^ 1
	var cohorts: Dictionary = ExpressionSystemScript.protein_cohorts_for_locus(expression, 1)
	_assert_true(cohorts.has(old_signature), "inherited proteins retain ancestral coding signature after daughter DNA mutation")
	_assert_true(not cohorts.has(int(gene.protein_signature)), "DNA mutation does not instantaneously manufacture mutant protein")

func _test_new_allele_builds_new_protein_cohort_over_time() -> void:
	var ancestral_genome = GenomeScript.create_ancestor()
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(ancestral_genome, config)
	var mutant = ancestral_genome.deep_copy()
	var gene = mutant.get_gene_by_locus(1)
	var old_signature: int = int(gene.protein_signature)
	gene.protein_signature = old_signature ^ 1
	var new_signature: int = int(gene.protein_signature)
	var pools: Dictionary = _rich_expression_pools()
	var rng = DeterministicRngScript.new(424242)
	for _i in range(300):
		ExpressionSystemScript.step(expression, mutant, pools, config.tick_dt_min, rng, config)
	var mrna_cohorts: Dictionary = expression[1]["mrna"]
	var protein_cohorts: Dictionary = expression[1]["protein"]
	_assert_true(mrna_cohorts.has(new_signature) or protein_cohorts.has(new_signature), "mutant DNA creates a new-sequence molecular cohort through ordinary expression")
	_assert_true(not protein_cohorts.has(old_signature) or float(protein_cohorts[old_signature]) < float(config.expression_reference_protein_count), "ancestral protein cohort can decay after coding mutation instead of being rewritten")

func _test_stochastic_partition_conserves_molecules_and_separates_sisters() -> void:
	var genome = GenomeScript.create_ancestor()
	var state: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var mrna_before: float = ExpressionSystemScript.total_mrna(state)
	var protein_before: float = ExpressionSystemScript.total_protein(state)
	var parts: Array = ExpressionSystemScript.partition(state, 0.5, DeterministicRngScript.new(9090), config)
	_assert_close(ExpressionSystemScript.total_mrna(parts[0]) + ExpressionSystemScript.total_mrna(parts[1]), mrna_before, 1e-9, "division partition conserves mRNA")
	_assert_close(ExpressionSystemScript.total_protein(parts[0]) + ExpressionSystemScript.total_protein(parts[1]), protein_before, 1e-9, "division partition conserves protein")
	_assert_true(absf(ExpressionSystemScript.checksum(parts[0]) - ExpressionSystemScript.checksum(parts[1])) > 1e-8, "stochastic molecular partition can separate otherwise identical sisters")

func _test_poisson_sampler_matches_expected_mean() -> void:
	var rng = DeterministicRngScript.new(67890)
	var lambda_value: float = 2.5
	var trials: int = 20000
	var total: int = 0
	for _i in range(trials):
		total += rng.poisson(lambda_value)
	var observed_mean: float = float(total) / float(trials)
	var mean_sigma: float = sqrt(lambda_value / float(trials))
	_assert_true(absf(observed_mean - lambda_value) <= 5.0 * mean_sigma, "deterministic Poisson sampler has expected ensemble mean")

func _test_resource_scarcity_scales_synthesis_without_negative_pools() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var pools: Dictionary = _rich_expression_pools()
	pools["ATP"] = 0.0001
	pools["NUC"] = 0.0001
	pools["AA"] = 0.0001
	var rng = DeterministicRngScript.new(121212)
	var saw_scale: bool = false
	for _i in range(50):
		var summary: Dictionary = ExpressionSystemScript.step(expression, genome, pools, config.tick_dt_min, rng, config)
		if float(summary["tx_scale"]) < 1.0 or float(summary["translation_scale"]) < 1.0:
			saw_scale = true
		MetabolicSolverScript.assert_nonnegative(pools)
	_assert_true(saw_scale, "scarce ATP/material invokes proportional synthesis scaling")
	_assert_true(float(pools["ATP"]) >= -1e-12 and float(pools["NUC"]) >= -1e-12 and float(pools["AA"]) >= -1e-12, "expression scarcity never creates negative resource pools")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
