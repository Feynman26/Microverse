extends SceneTree

const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
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
	_test_oxygen_is_the_unique_compatible_ligand_for_assay_regulator()
	_test_responsive_r03_is_inducible_but_constitutive_r03_is_not()
	_test_environment_phase_matches_slow_protein_response_timescale()
	_test_r03_proteome_lags_after_low_to_high_oxygen_transition()
	_test_lineage_selection_depends_reproducibly_on_environment_dynamics()

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

func _test_oxygen_is_the_unique_compatible_ligand_for_assay_regulator() -> void:
	var config = ExperimentScript.create_config(1)
	var valid: bool = true
	for metabolite_id in MetaboliteCatalogScript.ids():
		var distance: int = _hamming_distance(
			ExperimentScript.OXYGEN_LIGAND_SIGNATURE,
			MetaboliteCatalogScript.ligand_signature(metabolite_id)
		)
		if String(metabolite_id) == "O2":
			valid = valid and distance == 0
		else:
			valid = valid and distance > int(config.allosteric_max_distance)
	_assert_true(valid, "O2 is the only modeled ligand inside the M5-C regulator allosteric radius")

func _regulation_factor(genome, internal_oxygen: float, config) -> float:
	var expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config)
	var pools: Dictionary = _rich_expression_pools(config, internal_oxygen)
	var summary: Dictionary = ExpressionSystemScript.step(expression, genome, pools, 0.01, MeanRng.new(), config)
	return float(summary["regulation"][3])

func _rich_expression_pools(config, internal_oxygen: float) -> Dictionary:
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 1000.0
	pools["ADP"] = 100.0
	pools["AA"] = 1000.0
	pools["NUC"] = 1000.0
	pools["O2"] = internal_oxygen
	return pools

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

func _test_environment_phase_matches_slow_protein_response_timescale() -> void:
	var config = ExperimentScript.create_config(3)
	var phase_minutes: float = float(ExperimentScript.PHASE_TICKS) * float(config.tick_dt_min)
	var mrna_tau_minutes: float = 1.0 / float(config.mrna_decay_rate_per_min)
	var protein_tau_minutes: float = 1.0 / float(config.protein_decay_rate_per_min)
	_assert_close(phase_minutes, protein_tau_minutes, 1e-12, "M5-C O2 phase duration equals one protein turnover time constant")
	_assert_true(mrna_tau_minutes < phase_minutes, "mRNA can respond substantially faster than the encoded protein pool")

func _test_r03_proteome_lags_after_low_to_high_oxygen_transition() -> void:
	var config = ExperimentScript.create_config(4)
	var genomes: Dictionary = ExperimentScript.create_competitor_genomes()
	var constitutive = genomes["constitutive"]
	var responsive = genomes["responsive"]
	var constitutive_expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(constitutive, config)
	var responsive_expression: Dictionary = ExpressionSystemScript.create_equilibrium_state(responsive, config)
	var constitutive_pools: Dictionary = _rich_expression_pools(config, 0.0)
	var responsive_pools: Dictionary = _rich_expression_pools(config, 0.0)
	var mean_rng = MeanRng.new()
	for _i in range(ExperimentScript.PHASE_TICKS):
		constitutive_pools["O2"] = 0.0
		responsive_pools["O2"] = 0.0
		ExpressionSystemScript.step(constitutive_expression, constitutive, constitutive_pools, config.tick_dt_min, mean_rng, config)
		ExpressionSystemScript.step(responsive_expression, responsive, responsive_pools, config.tick_dt_min, mean_rng, config)
	var constitutive_gene = constitutive.get_gene_by_locus(3)
	var responsive_gene = responsive.get_gene_by_locus(3)
	var low_constitutive: float = ExpressionSystemScript.current_gene_protein(constitutive_expression, constitutive_gene)
	var low_responsive: float = ExpressionSystemScript.current_gene_protein(responsive_expression, responsive_gene)
	_assert_true(low_responsive < low_constitutive, "one anoxic phase leaves less R03 protein in the responsive architecture")

	for _i in range(20):
		constitutive_pools["O2"] = 6.0
		responsive_pools["O2"] = 6.0
		ExpressionSystemScript.step(constitutive_expression, constitutive, constitutive_pools, config.tick_dt_min, mean_rng, config)
		ExpressionSystemScript.step(responsive_expression, responsive, responsive_pools, config.tick_dt_min, mean_rng, config)
	var early_high_constitutive: float = ExpressionSystemScript.current_gene_protein(constitutive_expression, constitutive_gene)
	var early_high_responsive: float = ExpressionSystemScript.current_gene_protein(responsive_expression, responsive_gene)
	_assert_true(early_high_responsive < early_high_constitutive, "R03 protein remains lagged during the early high-O2 opportunity despite immediate sensing")

func _test_lineage_selection_depends_reproducibly_on_environment_dynamics() -> void:
	# Population discovery and first population confirmation have already been
	# observed and are excluded from this gate. These six sequential seeds were
	# declared before the lineage panel was executed. Four equally spaced phase
	# offsets remove the arbitrary choice of where a fluctuating square wave starts.
	# The criterion is intentionally two-sided: M5-C asks whether temporal
	# environment changes selection on regulatory architecture, not which genotype
	# must win.
	var seeds: Array = [13001, 13002, 13003, 13004, 13005, 13006]
	var offsets: Array = [0, 100, 200, 300]
	var panel: Dictionary = ExperimentScript.run_lineage_selection_panel(
		seeds,
		offsets,
		ExperimentScript.DEFAULT_MAX_TICKS
	)
	var seed_results: Array = panel["seeds"]
	var differentials: Array = panel["normalized_differentials"]
	var all_divided: bool = true
	for seed_result in seed_results:
		var stable_c: Dictionary = seed_result["stable_constitutive"]
		var stable_r: Dictionary = seed_result["stable_responsive"]
		all_divided = all_divided and bool(stable_c["reached_division"]) and bool(stable_r["reached_division"])
		for fluct in seed_result["fluctuating_runs"]:
			all_divided = all_divided and bool(fluct["constitutive"]["reached_division"])
			all_divided = all_divided and bool(fluct["responsive"]["reached_division"])
		print("M5-C LINEAGE seed=%d stable_adv=%.8f fluct_adv=%.8f normalized_delta=%.6f" % [
			int(seed_result["seed"]),
			float(seed_result["stable_advantage"]),
			float(seed_result["fluctuating_advantage"]),
			float(seed_result["normalized_differential"])
		])

	var mean_differential: float = float(panel["mean_normalized_differential"])
	var direction: float = 1.0 if mean_differential > 0.0 else -1.0
	var same_direction: int = 0
	for value_variant in differentials:
		var value: float = float(value_variant)
		if value * direction > 0.0:
			same_direction += 1
	print("M5-C LINEAGE mean_normalized_delta=%.6f same_direction=%d/%d" % [
		mean_differential, same_direction, differentials.size()
	])

	_assert_true(all_divided, "all paired M5-C lineage trajectories reach ordinary division criteria before timeout")
	_assert_true(same_direction >= 5, "environment-selection differential has the same direction in at least five of six unseen lineage seeds")
	_assert_true(absf(mean_differential) >= 0.02, "stable-vs-fluctuating dynamics shift the regulatory reproductive selection gradient by at least two percent of baseline growth rate")

func _hamming_distance(first_signature: int, second_signature: int) -> int:
	var value: int = (first_signature ^ second_signature) & 0xFFFF
	var distance: int = 0
	while value != 0:
		distance += value & 1
		value >>= 1
	return distance

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
