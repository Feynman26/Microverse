extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_translation_is_capped_by_shared_ribosome_budget()
	_test_rejected_translation_spends_no_atp_or_amino_acid()
	_test_overexpression_uses_capacity_needed_by_another_locus()
	_test_ribosome_allocation_is_genome_order_independent()

	if failures == 0:
		print("PASS: %d M5 shared-translation-capacity tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M5 shared-translation-capacity tests failed" % [failures, tests_run])
		quit(1)

func _config():
	var config = SimConfigScript.new()
	config.regulation_enabled = false
	config.allostery_enabled = false
	config.transcription_max_events_per_min = 0.0
	config.proteome_capacity_reference_units = 1.0
	config.translation_capacity_fraction_of_proteome_per_min = 0.05
	return config

func _genome():
	return GenomeScript.new([
		GeneScript.new(1, 5000, 0x1357, 1, 0xF139),
		GeneScript.new(2, 5000, 0x2468, 2, 0xF139)
	])

func _state(genome, mrna_1: float, mrna_2: float) -> Dictionary:
	var state: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, _config())
	var gene_1 = genome.get_gene_by_locus(1)
	var gene_2 = genome.get_gene_by_locus(2)
	state[1]["mrna"] = {int(gene_1.protein_signature): mrna_1}
	state[2]["mrna"] = {int(gene_2.protein_signature): mrna_2}
	state[1]["protein"] = {}
	state[2]["protein"] = {}
	return state

func _pools() -> Dictionary:
	return {"ATP": 100.0, "ADP": 0.0, "AA": 100.0, "NUC": 100.0}

func _test_translation_is_capped_by_shared_ribosome_budget() -> void:
	var config = _config()
	var genome = _genome()
	var state: Dictionary = _state(genome, 400.0, 400.0)
	var pools: Dictionary = _pools()
	var summary: Dictionary = ExpressionSystemScript.step(state, genome, pools, config.tick_dt_min, DeterministicRngScript.new(7001), config)
	var expected_capacity: float = (
		config.proteome_capacity_reference_units
		* config.expression_reference_protein_count
		* config.translation_capacity_fraction_of_proteome_per_min
		* config.tick_dt_min
	)
	_assert_close(float(summary["translation_capacity"]), expected_capacity, 1e-12, "translation capacity is derived from shared proteome size and elapsed time")
	_assert_close(float(summary["translated"]), expected_capacity, 1e-9, "abundant mRNA saturates but cannot exceed shared translation capacity")
	_assert_true(float(summary["ribosome_scale"]) < 1.0, "ribosome scale reports competition when proposals exceed capacity")

func _test_rejected_translation_spends_no_atp_or_amino_acid() -> void:
	var config = _config()
	var genome = _genome()
	var state: Dictionary = _state(genome, 400.0, 400.0)
	var pools: Dictionary = _pools()
	var atp_before: float = float(pools["ATP"])
	var aa_before: float = float(pools["AA"])
	var summary: Dictionary = ExpressionSystemScript.step(state, genome, pools, config.tick_dt_min, DeterministicRngScript.new(7002), config)
	var accepted: float = float(summary["translated"])
	_assert_close(atp_before - float(pools["ATP"]), accepted * config.translation_atp_cost_per_event, 1e-9, "only accepted ribosomal events spend translation ATP")
	_assert_close(aa_before - float(pools["AA"]), accepted * config.translation_aa_cost_per_event, 1e-9, "only accepted ribosomal events consume amino-acid material")

func _test_overexpression_uses_capacity_needed_by_another_locus() -> void:
	var config = _config()
	var baseline = _genome()
	var overexpressor = _genome()
	var baseline_state: Dictionary = _state(baseline, 100.0, 100.0)
	var over_state: Dictionary = _state(overexpressor, 400.0, 100.0)
	var baseline_pools: Dictionary = _pools()
	var over_pools: Dictionary = _pools()
	ExpressionSystemScript.step(baseline_state, baseline, baseline_pools, config.tick_dt_min, DeterministicRngScript.new(7003), config)
	ExpressionSystemScript.step(over_state, overexpressor, over_pools, config.tick_dt_min, DeterministicRngScript.new(7003), config)
	var gene_2 = baseline.get_gene_by_locus(2)
	var baseline_gene_2: float = ExpressionSystemScript.current_gene_protein(baseline_state, gene_2)
	var over_gene_2: float = ExpressionSystemScript.current_gene_protein(over_state, overexpressor.get_gene_by_locus(2))
	_assert_true(over_gene_2 < baseline_gene_2, "overexpressing one locus consumes shared translation capacity and reduces another locus output")

func _test_ribosome_allocation_is_genome_order_independent() -> void:
	var config = _config()
	var genome_a = _genome()
	var genome_b = genome_a.deep_copy()
	genome_b.genes.reverse()
	var state_a: Dictionary = _state(genome_a, 250.0, 150.0)
	var state_b: Dictionary = _state(genome_b, 250.0, 150.0)
	var pools_a: Dictionary = _pools()
	var pools_b: Dictionary = _pools()
	var summary_a: Dictionary = ExpressionSystemScript.step(state_a, genome_a, pools_a, config.tick_dt_min, DeterministicRngScript.new(7004), config)
	var summary_b: Dictionary = ExpressionSystemScript.step(state_b, genome_b, pools_b, config.tick_dt_min, DeterministicRngScript.new(7004), config)
	_assert_close(ExpressionSystemScript.checksum(state_a), ExpressionSystemScript.checksum(state_b), 1e-12, "shared ribosome allocation is invariant to genome array order")
	_assert_close(float(summary_a["translated"]), float(summary_b["translated"]), 1e-12, "genome order cannot change accepted translation throughput")
	_assert_close(float(pools_a["ATP"]), float(pools_b["ATP"]), 1e-12, "genome order cannot change translation ATP spending")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])