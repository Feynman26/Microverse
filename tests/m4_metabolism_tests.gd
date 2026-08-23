extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
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
	_test_reaction_order_does_not_change_state()
	_test_metabolism_conserves_structural_units_and_currency_pairs()
	_test_no_protein_has_zero_catalytic_flux()
	_test_oxygen_changes_metabolic_strategy_without_named_behavior()
	_test_one_coding_mutation_and_matching_protein_activate_dormant_waste_recovery()
	_test_biomass_is_assembled_catalytically()

	if failures == 0:
		print("PASS: %d M4 metabolic-integration tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M4 metabolic-integration tests failed" % [failures, tests_run])
		quit(1)

func _rich_pools() -> Dictionary:
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["G"] = 6.0
	pools["O2"] = 6.0
	pools["NH4"] = 5.0
	pools["P"] = 5.0
	pools["ATP"] = 2.0
	pools["ADP"] = 12.0
	pools["NAD"] = 8.0
	pools["NADH"] = 2.0
	return pools

func _equilibrium_expression(genome) -> Dictionary:
	return ExpressionSystemScript.create_equilibrium_state(genome, config)

func _test_reaction_order_does_not_change_state() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = _equilibrium_expression(genome)
	var forward: Dictionary = _rich_pools()
	var reverse: Dictionary = _rich_pools()
	var reversed_reactions: Array = reactions.duplicate()
	reversed_reactions.reverse()
	MetabolicSolverScript.step(forward, genome, expression, reactions, 0.5, 1.0, config)
	MetabolicSolverScript.step(reverse, genome, expression, reversed_reactions, 0.5, 1.0, config)
	var equal: bool = true
	for metabolite_id in forward.keys():
		if absf(float(forward[metabolite_id]) - float(reverse[metabolite_id])) > 1e-10:
			equal = false
			break
	_assert_true(equal, "reversing reaction array does not change metabolic state")

func _test_metabolism_conserves_structural_units_and_currency_pairs() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = _equilibrium_expression(genome)
	var pools: Dictionary = _rich_pools()
	var structural_before: Dictionary = MetabolicSolverScript.structural_totals(pools)
	var adenylate_before: float = float(pools["ATP"]) + float(pools["ADP"])
	var redox_before: float = float(pools["NAD"]) + float(pools["NADH"])
	MetabolicSolverScript.step(pools, genome, expression, reactions, 2.0, 1.0, config)
	var structural_after: Dictionary = MetabolicSolverScript.structural_totals(pools)
	_assert_close(float(structural_after["C"]), float(structural_before["C"]), 1e-9, "intracellular reactions conserve structural carbon")
	_assert_close(float(structural_after["N"]), float(structural_before["N"]), 1e-9, "intracellular reactions conserve structural nitrogen")
	_assert_close(float(structural_after["P"]), float(structural_before["P"]), 1e-9, "intracellular reactions conserve structural phosphorus")
	_assert_close(float(pools["ATP"]) + float(pools["ADP"]), adenylate_before, 1e-9, "reaction network conserves ATP+ADP currency pool")
	_assert_close(float(pools["NAD"]) + float(pools["NADH"]), redox_before, 1e-9, "reaction network conserves NAD+NADH redox pool")

func _test_no_protein_has_zero_catalytic_flux() -> void:
	var genome = GenomeScript.create_ancestor()
	var silent: Dictionary = _equilibrium_expression(genome)
	for locus_id in silent.keys():
		silent[locus_id]["protein"] = {}
	var pools: Dictionary = _rich_pools()
	var before: Dictionary = pools.duplicate(true)
	var fluxes: Dictionary = MetabolicSolverScript.step(pools, genome, silent, reactions, 1.0, 1.0, config)
	var total_flux: float = 0.0
	for value in fluxes.values():
		total_flux += float(value)
	_assert_close(total_flux, 0.0, 1e-12, "no protein molecules means zero catalytic flux")
	_assert_true(pools == before, "zero catalytic flux leaves metabolite pools unchanged")

func _test_oxygen_changes_metabolic_strategy_without_named_behavior() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = _equilibrium_expression(genome)
	var aerobic: Dictionary = _rich_pools()
	var hypoxic: Dictionary = _rich_pools()
	hypoxic["O2"] = 0.0
	var aerobic_flux: Dictionary = MetabolicSolverScript.step(aerobic, genome, expression, reactions, 2.0, 1.0, config)
	var hypoxic_flux: Dictionary = MetabolicSolverScript.step(hypoxic, genome, expression, reactions, 2.0, 1.0, config)
	_assert_true(float(aerobic_flux["R03"]) > 0.0, "oxygen permits oxidative phosphorylation through ordinary substrate availability")
	_assert_close(float(hypoxic_flux["R03"]), 0.0, 1e-12, "oxidative phosphorylation stops without oxygen")
	_assert_true(float(hypoxic_flux["R04"]) > 0.0, "ancestral weak fermentation remains available in hypoxia")
	_assert_true(float(aerobic["ATP"]) > float(hypoxic["ATP"]), "oxygen-rich chemistry yields a larger ATP pool than hypoxic chemistry")

func _test_one_coding_mutation_and_matching_protein_activate_dormant_waste_recovery() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var mutant = ancestor.deep_copy()
	var recovery = ReactionCatalogScript.by_id(reactions, "R05")
	var gene4 = mutant.get_gene_by_locus(4)
	var difference: int = (int(gene4.protein_signature) ^ int(recovery.signature)) & 0xFFFF
	var chosen_bit: int = -1
	for bit_index in range(16):
		if (difference & (1 << bit_index)) != 0:
			chosen_bit = bit_index
			break
	assert(chosen_bit >= 0)
	gene4.protein_signature = int(gene4.protein_signature) ^ (1 << chosen_bit)
	_assert_true(CatalyticLandscapeScript.hamming_distance(gene4.protein_signature, recovery.signature) == 4, "single coding mutation moves W1 recovery into active radius")

	var ancestral_pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	var mutant_pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	for pools in [ancestral_pools, mutant_pools]:
		pools["W1"] = 3.0
		pools["NAD"] = 4.0
		pools["ADP"] = 6.0
		pools["ATP"] = 0.0
		pools["NADH"] = 0.0
	var ancestor_expression: Dictionary = _equilibrium_expression(ancestor)
	var mutant_expression: Dictionary = _equilibrium_expression(mutant)
	var ancestor_flux: Dictionary = MetabolicSolverScript.step(ancestral_pools, ancestor, ancestor_expression, reactions, 1.0, 1.0, config)
	var mutant_flux: Dictionary = MetabolicSolverScript.step(mutant_pools, mutant, mutant_expression, reactions, 1.0, 1.0, config)
	_assert_close(float(ancestor_flux["R05"]), 0.0, 1e-12, "ancestral proteome cannot recover W1")
	_assert_true(float(mutant_flux["R05"]) > 0.0, "mutant protein cohort creates measurable W1 recovery flux")
	_assert_true(float(mutant_pools["C2"]) > float(ancestral_pools["C2"]), "new catalytic activity creates downstream carbon intermediate")

func _test_biomass_is_assembled_catalytically() -> void:
	var genome = GenomeScript.create_ancestor()
	var expression: Dictionary = _equilibrium_expression(genome)
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["AA"] = 4.0
	pools["LIP"] = 2.0
	pools["NUC"] = 4.0
	pools["ATP"] = 6.0
	pools["ADP"] = 4.0
	var before_bio: float = float(pools["BIO"])
	var fluxes: Dictionary = MetabolicSolverScript.step(pools, genome, expression, reactions, 1.0, 1.0, config)
	_assert_true(float(fluxes["R12"]) > 0.0, "ancestral biomass-assembly protein produces catalytic flux")
	_assert_true(float(pools["BIO"]) > before_bio, "BIO increases only through precursor-consuming catalytic assembly")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
