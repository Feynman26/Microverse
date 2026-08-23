extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ExpressionSolverScript = preload("res://src/expression/expression_solver.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_ancestor_is_outside_allosteric_radius()
	_test_semantically_neutral_x_can_control_biomass_flux()
	if failures == 0:
		print("PASS: %d M5 sensing/information tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5 sensing tests failed" % [failures, tests_run])
		quit(1)

func _test_ancestor_is_outside_allosteric_radius() -> void:
	var enabled = SimConfigScript.new()
	enabled.expression_noise_fraction = 0.0
	var disabled = SimConfigScript.new()
	disabled.expression_noise_fraction = 0.0
	disabled.allostery_enabled = false
	var genome = GenomeScript.create_ancestor()
	var a = ExpressionSolverScript.initialize(genome, enabled)
	var b = ExpressionSolverScript.initialize(genome, disabled)
	var pools_a: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, enabled)
	var pools_b: Dictionary = pools_a.duplicate(true)
	# Populate every ordinary intracellular pool. If any ancestral protein were
	# accidentally within the active ligand radius, regulation would diverge.
	for key in pools_a.keys():
		if String(key) != "BIO":
			pools_a[key] = maxf(1.0, float(pools_a[key]))
			pools_b[key] = pools_a[key]
	var stats_a: Dictionary = ExpressionSolverScript.step(a, genome, pools_a, 0.1, DeterministicRngScript.new(1), enabled)
	var stats_b: Dictionary = ExpressionSolverScript.step(b, genome, pools_b, 0.1, DeterministicRngScript.new(1), disabled)
	_assert_true(stats_a["regulation"] == stats_b["regulation"], "new ligand physics does not gift the ancestor an allosteric response")

func _test_semantically_neutral_x_can_control_biomass_flux() -> void:
	var config = SimConfigScript.new()
	config.expression_noise_fraction = 0.0
	config.regulatory_gain = 1.0
	config.allosteric_gain = 1.0
	var genome = GenomeScript.create_ancestor()
	# Turn locus 1 into an X-compatible generic regulator. X has signature 0x9696.
	# High bit is set -> repressor; bit 14 is clear -> X binding potentiates it.
	genome.get_gene_by_locus(1).protein_signature = 0x9696
	# Locus 10 is the exact R12 structural-biomass catalyst (0xAF14). Make its
	# promoter susceptible to the same regulator. Nothing in the engine knows
	# that R12 means growth or that X is a signal.
	genome.get_gene_by_locus(10).regulatory_signature = 0x9696

	var low_state = ExpressionSolverScript.initialize(genome, config)
	var high_state = low_state.deep_copy()
	var low_pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	var high_pools: Dictionary = low_pools.duplicate(true)
	for pools in [low_pools, high_pools]:
		pools["ATP"] = 40.0
		pools["ADP"] = 20.0
		pools["AA"] = 20.0
		pools["LIP"] = 10.0
		pools["NUC"] = 20.0
	low_pools["X"] = 0.0
	high_pools["X"] = 10.0

	var low_rng = DeterministicRngScript.new(2026)
	var high_rng = DeterministicRngScript.new(2026)
	for _i in range(150):
		ExpressionSolverScript.step(low_state, genome, low_pools, 0.1, low_rng, config)
		ExpressionSolverScript.step(high_state, genome, high_pools, 0.1, high_rng, config)

	_assert_true(high_state.protein_for(10) < low_state.protein_for(10), "neutral X can become information by altering evolved regulator abundance effects")

	# Compare catalytic consequence from identical precursor pools after the
	# regulatory phase. Explicit proteomes are the only difference.
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var assay_low: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	var assay_high: Dictionary = assay_low.duplicate(true)
	for pools in [assay_low, assay_high]:
		pools["ATP"] = 20.0
		pools["ADP"] = 10.0
		pools["AA"] = 10.0
		pools["LIP"] = 5.0
		pools["NUC"] = 10.0
	var low_flux: Dictionary = MetabolicSolverScript.step(assay_low, genome, reactions, 0.5, 1.0, config, low_state.proteins)
	var high_flux: Dictionary = MetabolicSolverScript.step(assay_high, genome, reactions, 0.5, 1.0, config, high_state.proteins)
	_assert_true(float(high_flux["R12"]) < float(low_flux["R12"]), "information use propagates from ligand -> regulation -> protein -> biomass reaction flux")
	_assert_true(float(assay_high["BIO"]) < float(assay_low["BIO"]), "the same genotype can produce environment-conditioned biomass outcome through generic molecular rules")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)
