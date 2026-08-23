extends SceneTree

const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const ExperimentScript = preload("res://src/experiments/m5c_regulatory_selection.gd")

class MeanRng:
	func poisson(lambda_value: float) -> float:
		return maxf(0.0, lambda_value)

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_oxygen_schedules_match_empty_cell_mean_transport_opportunity()
	_test_competitors_differ_only_in_target_regulatory_motif()
	_test_oxygen_regulator_has_no_hidden_m4_catalytic_activity()
	_test_responsive_r03_is_inducible_but_constitutive_r03_is_not()
	_test_fluctuating_environment_selects_regulatory_architecture_differently()

	if failures == 0:
		print("PASS: %d M5-C regulatory-selection tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5-C tests failed" % [failures, tests_run])
		quit(1)

func _test_oxygen_schedules_match_empty_cell_mean_transport_opportunity() -> void:
	var config = ExperimentScript.create_config(1)
	var means: Dictionary = ExperimentScript.schedule_transport_means(config)
	_assert_close(float(means["stable"]), float(means["fluctuating"]), 1e-12, "stable and fluctuating O2 schedules match mean empty-cell transport opportunity")

func _test_competitors_differ_only_in_target_regulatory_motif() -> void:
	var genomes: Dictionary = ExperimentScript.create_competitor_genomes()
	var constitutive = genomes["constitutive"]
	var responsive = genomes["responsive"]
	var differences: int = 0
	var valid: bool = constitutive.gene_count() == responsive.gene_count()
	for constitutive_gene in constitutive.genes:
		var responsive_gene = responsive.get_gene_by_locus(int(constitutive_gene.locus_id))
		valid = valid and responsive_gene != null
		valid = valid and int(constitutive_gene.promoter_code) == int(responsive_gene.promoter_code)
		valid = valid and int(constitutive_gene.protein_signature) == int(responsive_gene.protein_signature)
		valid = valid and int(constitutive_gene.neutral_marker) == int(responsive_gene.neutral_marker)
		if int(constitutive_gene.regulatory_signature) != int(responsive_gene.regulatory_signature):
			differences += 1
			valid = valid and int(constitutive_gene.locus_id) == 3
	_assert_true(valid and differences == 1, "M5-C competitors differ only in the R03 promoter regulatory motif")

func _test_oxygen_regulator_has_no_hidden_m4_catalytic_activity() -> void:
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var total_activity: float = 0.0
	for reaction in reactions:
		total_activity += CatalyticLandscapeScript.affinity(ExperimentScript.OXYGEN_LIGAND_SIGNATURE, int(reaction.signature))
	_assert_close(total_activity, 0.0, 1e-12, "O2-compatible regulator is outside every M4 catalytic activity radius")

func _regulation_factor(genome, internal_oxygen: float, config) -> float:
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 50.0
	pools["ADP"] = 10.0
	pools["AA"] = 20.0
	pools["NUC"] = 20.0
	pools["O2"] = internal_oxygen
	var summary: Dictionary = ExpressionSystemScript.step(expression, genome, pools, 0.01, MeanRng.new(), config)
	return float(summary["regulation"][3])

func _test_responsive_r03_is_inducible_but_constitutive_r03_is_not() -> void:
	var config = ExperimentScript.create_config(2)
	var genomes: Dictionary = ExperimentScript.create_competitor_genomes()
	var constitutive = genomes["constitutive"]
	var responsive = genomes["responsive"]
	var constitutive_low: float = _regulation_factor(constitutive, 0.0, config)
	var constitutive_high: float = _regulation_factor(constitutive, 6.0, config)
	var responsive_low: float = _regulation_factor(responsive, 0.0, config)
	var responsive_high: float = _regulation_factor(responsive, 6.0, config)
	_assert_close(constitutive_low, 1.0, 1e-12, "constitutive R03 promoter has no regulatory multiplier in low O2")
	_assert_close(constitutive_high, 1.0, 1e-12, "constitutive R03 promoter remains unregulated in high O2")
	_assert_true(responsive_low < responsive_high, "responsive R03 is more repressed without O2 than with O2")
	_assert_true(responsive_high > 0.90, "high O2 releases most R03 repression in the responsive architecture")

func _test_fluctuating_environment_selects_regulatory_architecture_differently() -> void:
	var seeds: Array = [1011, 2022, 3033, 4044]
	var result: Dictionary = ExperimentScript.run_paired_replicates(seeds, 3600)
	var stable: Array = result["stable"]
	var fluctuating: Array = result["fluctuating"]
	var paired: Array = result["paired_differences"]
	var positive_pairs: int = 0
	var viable: bool = true
	for i in range(seeds.size()):
		var s: Dictionary = stable[i]
		var f: Dictionary = fluctuating[i]
		print("M5-C seed=%d stable R=%d C=%d log_ratio=%.6f gen=%d | fluctuating R=%d C=%d log_ratio=%.6f gen=%d | delta=%.6f" % [
			int(seeds[i]), int(s["responsive"]), int(s["constitutive"]), float(s["log_ratio"]), int(s["max_generation"]),
			int(f["responsive"]), int(f["constitutive"]), float(f["log_ratio"]), int(f["max_generation"]), float(paired[i])
		])
		if float(paired[i]) > 0.0:
			positive_pairs += 1
		viable = viable and int(s["population"]) > ExperimentScript.FOUNDERS_PER_GENOTYPE * 2
		viable = viable and int(f["population"]) > ExperimentScript.FOUNDERS_PER_GENOTYPE * 2
		viable = viable and int(s["max_generation"]) >= 1 and int(f["max_generation"]) >= 1
	print("M5-C mean stable log_ratio=%.6f fluctuating=%.6f paired_delta=%.6f positive_pairs=%d/%d" % [
		float(result["mean_stable_log_ratio"]), float(result["mean_fluctuating_log_ratio"]),
		float(result["mean_paired_difference"]), positive_pairs, seeds.size()
	])
	_assert_true(viable, "both M5-C treatments support reproduction rather than comparing extinction artifacts")
	_assert_true(positive_pairs >= 3, "fluctuation increases responsive-vs-constitutive selection in at least three of four paired seeds")
	_assert_true(float(result["mean_paired_difference"]) > 0.05, "mean responsive selection is materially stronger under fluctuating than transport-matched stable O2")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
