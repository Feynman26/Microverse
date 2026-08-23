extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const Circuit = preload("res://src/experiments/m5c_ros_circuit.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_genomes_differ_only_at_ros_control_motif()
	_test_oxygen_regulator_is_not_a_hidden_catalyst()
	_test_r03_generates_ros_and_r09_is_exact_ros_control_catalyst()
	_test_oxygen_relief_of_repression_is_specific_to_responsive_r09()

	if failures == 0:
		print("PASS: %d M5-C ROS-circuit mechanics tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5-C ROS-circuit mechanics tests failed" % [failures, tests_run])
		quit(1)

func _test_genomes_differ_only_at_ros_control_motif() -> void:
	var genomes: Dictionary = Circuit.create_genomes()
	var constitutive = genomes["constitutive"]
	var responsive = genomes["responsive"]
	_assert_true(constitutive.gene_count() == responsive.gene_count(), "ROS validation genomes contain the same loci")
	var differences: int = 0
	for c_gene in constitutive.genes:
		var r_gene = responsive.get_gene_by_locus(int(c_gene.locus_id))
		_assert_true(r_gene != null, "responsive genome contains matching locus %d" % int(c_gene.locus_id))
		if r_gene == null:
			continue
		if int(c_gene.promoter_code) != int(r_gene.promoter_code): differences += 1
		if int(c_gene.protein_signature) != int(r_gene.protein_signature): differences += 1
		if int(c_gene.neutral_marker) != int(r_gene.neutral_marker): differences += 1
		if int(c_gene.regulatory_signature) != int(r_gene.regulatory_signature):
			differences += 1
			_assert_true(int(c_gene.locus_id) == Circuit.ROS_CONTROL_LOCUS, "the sole inherited difference is the R09 regulatory motif")
	_assert_true(differences == 1, "constitutive and responsive ROS circuits differ in exactly one inherited field")

func _test_oxygen_regulator_is_not_a_hidden_catalyst() -> void:
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var active_count: int = 0
	for reaction in reactions:
		if CatalyticLandscapeScript.affinity(Circuit.OXYGEN_REPRESSOR_SIGNATURE, int(reaction.signature)) > 0.0:
			active_count += 1
	_assert_true(active_count == 0, "O2-compatible regulator has no catalytic activity anywhere in the M4 reaction catalog")

func _test_r03_generates_ros_and_r09_is_exact_ros_control_catalyst() -> void:
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var r03 = ReactionCatalogScript.by_id(reactions, "R03")
	var r09 = ReactionCatalogScript.by_id(reactions, "R09")
	_assert_true(r03 != null and r09 != null, "R03 and R09 exist in the production reaction catalog")
	if r03 != null:
		_assert_true(float(r03.products.get("ROS", 0.0)) > 0.0 and float(r03.substrates.get("O2", 0.0)) > 0.0, "ordinary oxidative phosphorylation couples O2 use to ROS production")
	if r09 != null:
		_assert_true(int(r09.signature) == Circuit.ROS_CONTROL_PROTEIN_SIGNATURE, "R09 validation protein is the exact production ROS-control catalyst")
		_assert_true(float(r09.substrates.get("ROS", 0.0)) > 0.0 and float(r09.substrates.get("ATP", 0.0)) > 0.0, "R09 removes ROS through the existing ATP-consuming chemistry")

func _test_oxygen_relief_of_repression_is_specific_to_responsive_r09() -> void:
	var config = SimConfigScript.new()
	var genomes: Dictionary = Circuit.create_genomes()
	var responsive = genomes["responsive"]
	var constitutive = genomes["constitutive"]
	var responsive_state: Dictionary = ExpressionSystemScript.create_equilibrium_state(responsive, config)
	var constitutive_state: Dictionary = ExpressionSystemScript.create_equilibrium_state(constitutive, config)
	var low_pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	var high_pools: Dictionary = low_pools.duplicate(true)
	low_pools["O2"] = 0.0
	high_pools["O2"] = 6.0
	low_pools["ATP"] = 100.0
	high_pools["ATP"] = 100.0
	low_pools["AA"] = 100.0
	high_pools["AA"] = 100.0
	low_pools["NUC"] = 100.0
	high_pools["NUC"] = 100.0

	var resp_low_summary: Dictionary = ExpressionSystemScript.step(responsive_state.duplicate(true), responsive, low_pools.duplicate(true), config.tick_dt_min, DeterministicRngScript.new(230001), config)
	var resp_high_summary: Dictionary = ExpressionSystemScript.step(responsive_state.duplicate(true), responsive, high_pools.duplicate(true), config.tick_dt_min, DeterministicRngScript.new(230001), config)
	var const_low_summary: Dictionary = ExpressionSystemScript.step(constitutive_state.duplicate(true), constitutive, low_pools.duplicate(true), config.tick_dt_min, DeterministicRngScript.new(230001), config)
	var const_high_summary: Dictionary = ExpressionSystemScript.step(constitutive_state.duplicate(true), constitutive, high_pools.duplicate(true), config.tick_dt_min, DeterministicRngScript.new(230001), config)

	var responsive_low_factor: float = float(resp_low_summary["regulation"][Circuit.ROS_CONTROL_LOCUS])
	var responsive_high_factor: float = float(resp_high_summary["regulation"][Circuit.ROS_CONTROL_LOCUS])
	var constitutive_low_factor: float = float(const_low_summary["regulation"][Circuit.ROS_CONTROL_LOCUS])
	var constitutive_high_factor: float = float(const_high_summary["regulation"][Circuit.ROS_CONTROL_LOCUS])
	_assert_true(responsive_high_factor > responsive_low_factor, "high intracellular O2 relieves repression of responsive R09 expression")
	_assert_close(constitutive_low_factor, 1.0, 1e-12, "constitutive R09 has no low-O2 regulatory edge")
	_assert_close(constitutive_high_factor, 1.0, 1e-12, "constitutive R09 remains unregulated at high O2")
	_assert_true(absf(responsive_high_factor - responsive_low_factor) >= 0.05, "generic allostery creates a material O2-conditioned R09 transcription difference")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
