extends SceneTree

const Assay = preload("res://src/experiments/m5c_burden_circuit.gd")
const BaseExperiment = preload("res://src/experiments/m5c_regulatory_selection.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	var genomes: Dictionary = Assay.create_genomes()
	var c = genomes["constitutive"]
	var r = genomes["responsive"]
	_assert_true(c.gene_count() == r.gene_count(), "burden architectures contain the same loci")
	var differences: int = 0
	for gene in c.genes:
		var other = r.get_gene_by_locus(int(gene.locus_id))
		_assert_true(other != null, "responsive burden architecture contains matching locus %d" % int(gene.locus_id))
		if gene.canonical_key() != other.canonical_key():
			differences += 1
			_assert_true(Assay.BURDEN_LOCI.has(int(gene.locus_id)), "only neutral burden promoter motifs differ")
			_assert_true(int(gene.protein_signature) == int(other.protein_signature), "burden coding sequence is identical across architectures")
	_assert_true(differences == Assay.BURDEN_LOCI.size(), "architectures differ only at the three burden regulatory motifs")

	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var catalytic_total: float = 0.0
	for reaction in reactions:
		catalytic_total += CatalyticLandscapeScript.affinity(Assay.NEUTRAL_BURDEN_SIGNATURE, int(reaction.signature))
	_assert_close(catalytic_total, 0.0, 1e-12, "neutral burden protein has zero catalytic activity across M4")

	var config = BaseExperiment.create_config(99173)
	var state: Dictionary = ExpressionSystemScript.create_equilibrium_state(r, config)
	var low_pools: Dictionary = {"O2": 0.0}
	var high_pools: Dictionary = {"O2": 6.0}
	var target = r.get_gene_by_locus(11)
	var low_factor: float = ExpressionSystemScript._regulation_factor(target, state, low_pools, config)
	var high_factor: float = ExpressionSystemScript._regulation_factor(target, state, high_pools, config)
	_assert_true(low_factor < 0.75, "anoxia materially represses responsive neutral burden expression")
	_assert_true(high_factor > low_factor + 0.15, "high O2 allosterically releases neutral burden repression")
	var constitutive_target = c.get_gene_by_locus(11)
	var constitutive_low: float = ExpressionSystemScript._regulation_factor(constitutive_target, ExpressionSystemScript.create_equilibrium_state(c, config), low_pools, config)
	_assert_close(constitutive_low, 1.0, 1e-12, "constitutive burden locus has no O2 regulatory edge")

	if failures == 0:
		print("PASS: %d M5-C expression-burden mechanics tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5-C expression-burden mechanics tests failed" % [failures, tests_run])
		quit(1)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
