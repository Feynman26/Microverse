extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_initial_proteome_is_bounded_without_changing_relative_allocation()
	_test_initial_allocation_is_genome_order_independent()
	_test_overexpression_displaces_other_proteins_under_shared_budget()
	_test_runtime_overflow_recycles_material_but_not_expression_atp()

	if failures == 0:
		print("PASS: %d M5 finite-proteome tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5 finite-proteome tests failed" % [failures, tests_run])
		quit(1)

func _new_cell(genome, config):
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2.ZERO, config.ancestor_volume)
	cell.genome = genome.deep_copy()
	cell.initialize_molecular_state(config)
	return cell

func _combined_structural(cell, config) -> Dictionary:
	var chemical: Dictionary = MetabolicSolverScript.structural_totals(cell.metabolites)
	var stored: Dictionary = ExpressionSystemScript.structural_storage_totals(cell.expression_state, config)
	return {
		"C": float(chemical["C"]) + float(stored["C"]),
		"N": float(chemical["N"]) + float(stored["N"]),
		"P": float(chemical["P"]) + float(stored["P"])
	}

func _test_initial_proteome_is_bounded_without_changing_relative_allocation() -> void:
	var config = SimConfigScript.new()
	var genome = GenomeScript.create_ancestor()
	var unconstrained: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var desired_total: float = ExpressionSystemScript.total_protein(unconstrained)
	var cell = _new_cell(genome, config)
	var capacity: float = cell.proteome_capacity(config)
	var locus_1 = genome.get_gene_by_locus(1)
	var locus_2 = genome.get_gene_by_locus(2)
	var p1: float = ExpressionSystemScript.current_gene_protein(cell.expression_state, locus_1)
	var p2: float = ExpressionSystemScript.current_gene_protein(cell.expression_state, locus_2)
	var expected_ratio: float = float(locus_1.promoter_strength()) / float(locus_2.promoter_strength())
	_assert_true(desired_total > capacity, "ancestral unconstrained expression exceeds the finite default proteome budget")
	_assert_close(cell.total_protein(), capacity, 1e-9, "initial cell proteome is normalized exactly to the shared capacity")
	_assert_close(p1 / p2, expected_ratio, 1e-12, "capacity normalization preserves relative promoter-driven protein allocation")

func _test_initial_allocation_is_genome_order_independent() -> void:
	var config = SimConfigScript.new()
	var genome_a = GenomeScript.create_ancestor()
	var genome_b = genome_a.deep_copy()
	genome_b.genes.reverse()
	var cell_a = _new_cell(genome_a, config)
	var cell_b = _new_cell(genome_b, config)
	_assert_close(ExpressionSystemScript.checksum(cell_a.expression_state), ExpressionSystemScript.checksum(cell_b.expression_state), 1e-12, "finite proteome allocation is independent of genome array ordering")
	_assert_close(cell_a.total_protein(), cell_b.total_protein(), 1e-12, "gene ordering cannot change total allocated proteome")

func _test_overexpression_displaces_other_proteins_under_shared_budget() -> void:
	var config = SimConfigScript.new()
	config.proteome_capacity_reference_units = 0.75
	var baseline = GenomeScript.new([
		GeneScript.new(1, 5000, 0x1357, 1, 0xF139),
		GeneScript.new(2, 5000, 0x2468, 2, 0xF139)
	])
	var overexpressor = GenomeScript.new([
		GeneScript.new(1, 9000, 0x1357, 1, 0xF139),
		GeneScript.new(2, 5000, 0x2468, 2, 0xF139)
	])
	var baseline_cell = _new_cell(baseline, config)
	var over_cell = _new_cell(overexpressor, config)
	var baseline_gene2: float = ExpressionSystemScript.current_gene_protein(baseline_cell.expression_state, baseline.get_gene_by_locus(2))
	var over_gene2: float = ExpressionSystemScript.current_gene_protein(over_cell.expression_state, overexpressor.get_gene_by_locus(2))
	_assert_close(baseline_cell.total_protein(), baseline_cell.proteome_capacity(config), 1e-9, "baseline two-gene system fills the shared proteome budget")
	_assert_close(over_cell.total_protein(), over_cell.proteome_capacity(config), 1e-9, "overexpressing system cannot exceed the same shared proteome budget")
	_assert_true(over_gene2 < baseline_gene2, "overexpressing one locus displaces another protein instead of receiving free cellular capacity")

func _test_runtime_overflow_recycles_material_but_not_expression_atp() -> void:
	var config = SimConfigScript.new()
	config.proteome_capacity_reference_units = 1.5
	var genome = GenomeScript.new([
		GeneScript.new(1, 9000, 0x1357, 1, 0xF139),
		GeneScript.new(2, 9000, 0x2468, 2, 0xF139)
	])
	var cell = _new_cell(genome, config)
	cell.metabolites["ATP"] = 100.0
	cell.metabolites["ADP"] = 10.0
	cell.metabolites["AA"] = 20.0
	cell.metabolites["NUC"] = 20.0
	cell.metabolites["G"] = 0.0
	cell.metabolites["O2"] = 0.0
	cell.metabolites["NH4"] = 0.0
	cell.metabolites["P"] = 0.0
	var before: Dictionary = _combined_structural(cell, config)
	var rng = DeterministicRngScript.new(606060)
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var saw_overflow: bool = false
	var spent_atp: float = 0.0
	for _i in range(80):
		cell.step_intracellular(config.tick_dt_min, config, reactions, rng)
		spent_atp += float(cell.last_expression_summary.get("atp_spent", 0.0))
		if float(cell.last_expression_summary.get("proteome_removed", 0.0)) > 0.0:
			saw_overflow = true
		_assert_true(cell.total_protein() <= cell.proteome_capacity(config) + 1e-8, "runtime proteome never exceeds shared capacity after enforcement")
	var after: Dictionary = _combined_structural(cell, config)
	_assert_true(saw_overflow, "stochastic expression can hit the finite proteome ceiling")
	_assert_true(spent_atp > 0.0, "protein synthesis still spends ATP when excess protein is subsequently removed")
	_assert_close(float(after["C"]), float(before["C"]), 1e-8, "proteome turnover recycles modeled structural carbon")
	_assert_close(float(after["N"]), float(before["N"]), 1e-8, "proteome turnover recycles modeled structural nitrogen")
	_assert_close(float(after["P"]), float(before["P"]), 1e-8, "proteome turnover preserves modeled structural phosphorus")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])