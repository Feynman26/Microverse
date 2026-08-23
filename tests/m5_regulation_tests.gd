extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

# Expected-value harness used only for causal M5-B controls. Production biology
# remains stochastic Poisson expression; M5-A already validates that sampler.
class MeanRng:
	func poisson(lambda_value: float) -> float:
		return maxf(0.0, lambda_value)

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_ancestor_has_no_active_regulatory_edges()
	_test_ancestor_has_no_accidental_allosteric_matches()
	_test_one_motif_bit_can_create_first_regulatory_edge()
	_test_generic_protein_signatures_can_activate_or_repress()
	_test_regulatory_mutation_operator_is_one_bit_and_genotype_only()
	_test_neutral_x_modulates_only_compatible_regulator()
	_test_x_information_propagates_to_protein_and_biomass_flux()
	_test_historical_protein_cohort_preserves_regulatory_effect_after_dna_change()
	_test_active_regulatory_network_is_order_independent_and_replayable()

	if failures == 0:
		print("PASS: %d M5-B regulation/sensing tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5-B tests failed" % [failures, tests_run])
		quit(1)

func _base_config():
	var config = SimConfigScript.new()
	config.regulatory_gain = 1.0
	config.allosteric_gain = 1.0
	return config

func _rich_pools(config) -> Dictionary:
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 100.0
	pools["ADP"] = 20.0
	pools["NUC"] = 50.0
	pools["AA"] = 50.0
	return pools

func _hamming(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & 0xFFFF
	var result: int = 0
	while value != 0:
		result += value & 1
		value >>= 1
	return result

func _regulation_summary(genome, expression: Dictionary, pools: Dictionary, config) -> Dictionary:
	var local_state: Dictionary = expression.duplicate(true)
	var local_pools: Dictionary = pools.duplicate(true)
	var summary: Dictionary = ExpressionSystemScript.step(local_state, genome, local_pools, 0.01, MeanRng.new(), config)
	return summary["regulation"]

func _test_ancestor_has_no_active_regulatory_edges() -> void:
	var config = _base_config()
	var genome = GenomeScript.create_ancestor()
	var proteins: Array = []
	for gene in genome.genes:
		proteins.append(int(gene.protein_signature))

	var minimum_distance: int = 99
	var all_dormant: bool = true
	for target_gene in genome.genes:
		for signature in proteins:
			var distance: int = _hamming(int(signature), int(target_gene.regulatory_signature))
			minimum_distance = mini(minimum_distance, distance)
			if distance <= int(config.regulatory_max_distance):
				all_dormant = false
	_assert_true(all_dormant, "ancestral promoter motifs have no active protein-binding edges")
	_assert_true(minimum_distance == config.regulatory_max_distance + 1, "nearest ancestral regulatory opportunity is exactly one motif mutation outside active radius")

	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var factors: Dictionary = _regulation_summary(genome, expression, _rich_pools(config), config)
	var neutral: bool = true
	for factor in factors.values():
		if absf(float(factor) - 1.0) > 1e-12:
			neutral = false
	_assert_true(neutral, "enabling M5-B regulation does not alter ancestral basal transcription")

func _test_ancestor_has_no_accidental_allosteric_matches() -> void:
	var config = _base_config()
	var genome = GenomeScript.create_ancestor()
	var minimum_distance: int = 99
	for gene in genome.genes:
		for metabolite_id in MetaboliteCatalogScript.ids():
			minimum_distance = mini(
				minimum_distance,
				_hamming(int(gene.protein_signature), MetaboliteCatalogScript.ligand_signature(metabolite_id))
			)
	_assert_true(minimum_distance > int(config.allosteric_max_distance), "no ancestral protein is accidentally inside any ligand allosteric radius")

func _test_one_motif_bit_can_create_first_regulatory_edge() -> void:
	var config = _base_config()
	config.allostery_enabled = false
	var genome = GenomeScript.create_ancestor()
	var target = genome.get_gene_by_locus(1)
	var regulator_signature: int = int(genome.get_gene_by_locus(1).protein_signature)
	var before_distance: int = _hamming(regulator_signature, int(target.regulatory_signature))
	var before_expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var before_factor: float = float(_regulation_summary(genome, before_expression, _rich_pools(config), config)[1])

	# Locus 1 ancestral motif differs from its protein at bits 0..3. Flipping one
	# of those bits moves distance 4 -> 3, crossing the generic binding radius.
	target.regulatory_signature = int(target.regulatory_signature) ^ 0x0001
	var after_distance: int = _hamming(regulator_signature, int(target.regulatory_signature))
	var after_expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var after_factor: float = float(_regulation_summary(genome, after_expression, _rich_pools(config), config)[1])

	_assert_true(before_distance == 4 and after_distance == 3, "one ordinary motif bit flip can cross the promoter-binding threshold")
	_assert_close(before_factor, 1.0, 1e-12, "pre-mutation dormant motif has no regulatory effect")
	_assert_true(after_factor > 1.0, "newly accessible low-high-bit protein creates activation without a named regulator class")

func _test_generic_protein_signatures_can_activate_or_repress() -> void:
	var config = _base_config()
	config.allostery_enabled = false
	var genome = GenomeScript.create_ancestor()
	var target = genome.get_gene_by_locus(1)
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var pools: Dictionary = _rich_pools(config)

	target.regulatory_signature = int(genome.get_gene_by_locus(1).protein_signature)
	var activation: float = float(_regulation_summary(genome, expression, pools, config)[1])
	target.regulatory_signature = int(genome.get_gene_by_locus(8).protein_signature)
	var repression: float = float(_regulation_summary(genome, expression, pools, config)[1])

	_assert_true(activation > 1.0, "exact compatible protein with high bit clear activates a promoter")
	_assert_true(repression < 1.0, "exact compatible protein with high bit set represses a promoter")

func _test_regulatory_mutation_operator_is_one_bit_and_genotype_only() -> void:
	var config = _base_config()
	config.promoter_mutation_rate_per_gene = 0.0
	config.signature_mutation_rate_per_gene = 0.0
	config.regulatory_signature_mutation_rate_per_gene = 1.0
	config.neutral_marker_mutation_rate_per_gene = 0.0
	var parent = GenomeScript.create_ancestor()
	var parent_key: String = parent.canonical_key()
	var result: Dictionary = MutationEngineScript.new().mutate_copy(parent, DeterministicRngScript.new(51051), config)
	var child = result["genome"]
	var events: Array = result["events"]
	var valid: bool = events.size() == parent.gene_count()
	for event in events:
		valid = valid and String(event["mutation_type"]) == "regulatory_signature_bit_flip"
		valid = valid and _hamming(int(event["old_value"]), int(event["new_value"])) == 1
	for parent_gene in parent.genes:
		var child_gene = child.get_gene_by_locus(int(parent_gene.locus_id))
		valid = valid and int(child_gene.promoter_code) == int(parent_gene.promoter_code)
		valid = valid and int(child_gene.protein_signature) == int(parent_gene.protein_signature)
		valid = valid and int(child_gene.neutral_marker) == int(parent_gene.neutral_marker)
	_assert_true(parent.canonical_key() == parent_key, "regulatory mutation never rewrites the parent genotype")
	_assert_true(child.canonical_key() != parent_key, "regulatory motif mutation changes inherited network genotype")
	_assert_true(valid, "forced M5-B mutation changes each promoter motif by exactly one bit and no other gene field")

func _create_x_regulated_biomass_genome():
	# X signature is 0x9696. Its high bit is set (repressor), while bit 14 is
	# clear, so exact X binding potentiates the repressor under the digital rule.
	var regulator = GeneScript.new(1, 6200, 0x9696, 1, 0x0000)
	var biomass_catalyst = GeneScript.new(2, 6000, 0xAF14, 2, 0x9696)
	return GenomeScript.new([regulator, biomass_catalyst])

func _test_neutral_x_modulates_only_compatible_regulator() -> void:
	var config = _base_config()
	var genome = _create_x_regulated_biomass_genome()
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var low_pools: Dictionary = _rich_pools(config)
	var high_pools: Dictionary = low_pools.duplicate(true)
	low_pools["X"] = 0.0
	high_pools["X"] = 10.0
	var low_factor: float = float(_regulation_summary(genome, expression, low_pools, config)[2])
	var high_factor: float = float(_regulation_summary(genome, expression, high_pools, config)[2])

	var no_allostery = _base_config()
	no_allostery.allostery_enabled = false
	var disabled_factor: float = float(_regulation_summary(genome, expression, high_pools, no_allostery)[2])
	_assert_true(high_factor < low_factor, "high neutral-X concentration potentiates its compatible repressor and lowers target transcription factor")
	_assert_close(disabled_factor, low_factor, 1e-12, "X has no special semantic effect when generic allostery is disabled")

func _test_x_information_propagates_to_protein_and_biomass_flux() -> void:
	var config = _base_config()
	var genome = _create_x_regulated_biomass_genome()
	var low_expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var high_expression: Dictionary = low_expression.duplicate(true)
	var low_pools: Dictionary = _rich_pools(config)
	var high_pools: Dictionary = low_pools.duplicate(true)
	low_pools["X"] = 0.0
	high_pools["X"] = 10.0
	var low_rng = MeanRng.new()
	var high_rng = MeanRng.new()
	for _i in range(400):
		ExpressionSystemScript.step(low_expression, genome, low_pools, config.tick_dt_min, low_rng, config)
		ExpressionSystemScript.step(high_expression, genome, high_pools, config.tick_dt_min, high_rng, config)

	var target_gene = genome.get_gene_by_locus(2)
	var low_protein: float = ExpressionSystemScript.current_gene_protein(low_expression, target_gene)
	var high_protein: float = ExpressionSystemScript.current_gene_protein(high_expression, target_gene)
	_assert_true(high_protein < low_protein, "ligand-conditioned regulation propagates into different biomass-catalyst protein abundance")

	var assay_low: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	var assay_high: Dictionary = assay_low.duplicate(true)
	for pools in [assay_low, assay_high]:
		pools["ATP"] = 20.0
		pools["ADP"] = 10.0
		pools["AA"] = 10.0
		pools["LIP"] = 5.0
		pools["NUC"] = 10.0
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var low_flux: Dictionary = MetabolicSolverScript.step(assay_low, genome, low_expression, reactions, 0.25, 1.0, config)
	var high_flux: Dictionary = MetabolicSolverScript.step(assay_high, genome, high_expression, reactions, 0.25, 1.0, config)
	_assert_true(float(high_flux["R12"]) < float(low_flux["R12"]), "X information propagates through protein abundance into structural-biomass reaction flux")
	_assert_true(float(assay_high["BIO"]) < float(assay_low["BIO"]), "same genotype can produce environment-conditioned BIO outcome without a behavior flag")

func _test_historical_protein_cohort_preserves_regulatory_effect_after_dna_change() -> void:
	var config = _base_config()
	config.allostery_enabled = false
	var ancestral_regulator = GeneScript.new(1, 6200, 0x1357, 1, 0x0000)
	var target = GeneScript.new(2, 6000, 0xAF14, 2, 0x1357)
	var old_genome = GenomeScript.new([ancestral_regulator, target])
	var inherited_state: Dictionary = ExpressionSystemScript.create_equilibrium_state(old_genome, config)

	var mutated_genome = old_genome.deep_copy()
	mutated_genome.get_gene_by_locus(1).protein_signature = 0x2468
	var inherited_factor: float = float(_regulation_summary(mutated_genome, inherited_state, _rich_pools(config), config)[2])

	var rewritten_control: Dictionary = inherited_state.duplicate(true)
	var old_amount: float = float(rewritten_control[1]["protein"].get(0x1357, 0.0))
	rewritten_control[1]["protein"].erase(0x1357)
	rewritten_control[1]["protein"][0x2468] = old_amount
	var rewritten_factor: float = float(_regulation_summary(mutated_genome, rewritten_control, _rich_pools(config), config)[2])

	_assert_true(inherited_factor > 1.0, "old-sequence protein inherited after DNA mutation still regulates through its physical sequence")
	_assert_close(rewritten_factor, 1.0, 1e-12, "regulatory effect disappears if the historical protein cohort is physically replaced by the new distant sequence")

func _test_active_regulatory_network_is_order_independent_and_replayable() -> void:
	var config = _base_config()
	var genome_a = _create_x_regulated_biomass_genome()
	var genome_b = genome_a.deep_copy()
	genome_b.genes.reverse()
	var state_a: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome_a, config)
	var state_b: Dictionary = state_a.duplicate(true)
	var pools_a: Dictionary = _rich_pools(config)
	var pools_b: Dictionary = pools_a.duplicate(true)
	pools_a["X"] = 3.0
	pools_b["X"] = 3.0
	var rng_a = DeterministicRngScript.new(88221)
	var rng_b = DeterministicRngScript.new(88221)
	for _i in range(150):
		ExpressionSystemScript.step(state_a, genome_a, pools_a, config.tick_dt_min, rng_a, config)
		ExpressionSystemScript.step(state_b, genome_b, pools_b, config.tick_dt_min, rng_b, config)
	_assert_close(ExpressionSystemScript.checksum(state_a), ExpressionSystemScript.checksum(state_b), 1e-12, "active M5-B network is independent of genome-array ordering")
	_assert_true(pools_a == pools_b, "active regulation/allostery preserves exact resource replay under gene-order reversal")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
