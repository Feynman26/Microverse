extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const MutationDynamicsAnalyticsScript = preload("res://src/experiments/mutation_dynamics_analytics.gd")
const ExperimentAnalyticsScript = preload("res://src/experiments/experiment_analytics.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_population_sampling_is_observational()
	_test_event_ledger_reconstructs_realized_mutations_per_birth()
	_test_functional_identity_separates_neutral_history()
	_test_ancestral_departure_proxy_is_structural_not_fitness()
	_test_environment_shift_association_tracks_mutator_abundance_descriptively()
	_test_runner_detectors_use_m10_mechanistic_metrics()
	if failures == 0:
		print("PASS: %d M10 mutation-analytics tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 mutation-analytics tests failed" % [failures, tests_run])
		quit(1)

func _test_population_sampling_is_observational() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	config.mutation_enabled = false
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor()
	var checksum_before: float = sim.checksum()
	var rng_before: int = int(sim.rng.get_state())
	var sampled: Dictionary = MutationDynamicsAnalyticsScript.sample_population(sim)
	_assert_close(sim.checksum(), checksum_before, 1e-12, "mutation analytics cannot alter authoritative simulation checksum")
	_assert_true(int(sim.rng.get_state()) == rng_before, "mutation analytics consume no RNG draws")
	_assert_true(int(sampled["population"]) == 1 and sampled["lineages"].size() == 1, "population sample reports one ancestral lineage")
	_assert_true(int(sampled["genetic_genotype_count"]) == 1 and int(sampled["functional_genotype_count"]) == 1, "ancestral population has one genetic and one functional genotype")
	_assert_close(float(sampled["point_error_rate_per_gene"]["mean"]), config.baseline_point_error_rate_per_gene, 1e-12, "repair-naive ancestor reports baseline copy-error probability")
	_assert_close(float(sampled["gene_count"]["mean"]), 12.0, 1e-12, "population sample reports physical genome size")
	_assert_close(float(sampled["ancestral_departure"]["mean"]), 0.0, 1e-12, "canonical ancestor has zero molecular-departure proxy")

func _test_event_ledger_reconstructs_realized_mutations_per_birth() -> void:
	var events: Array = [
		{"kind": "birth", "cell_id": 1, "parent_id": -1},
		{
			"kind": "division", "parent_id": 1, "daughter_ids": [2, 3],
			"replication_profile": {
				"point_error_rate_per_gene": 0.002,
				"structural_error_rate_per_genome": 0.001,
				"mean_repair_activity": 0.0
			}
		},
		{"kind": "birth", "cell_id": 2, "parent_id": 1},
		{"kind": "birth", "cell_id": 3, "parent_id": 1},
		{"kind": "mutation", "cell_id": 2, "mutation_type": "protein_signature_bit_flip", "parent_genome_size": 12, "resulting_genome_size": 12},
		{"kind": "mutation", "cell_id": 2, "mutation_type": "gene_duplication", "parent_genome_size": 12, "resulting_genome_size": 13},
		{"kind": "mutation", "cell_id": 3, "mutation_type": "gene_duplication_blocked", "reason": "insufficient_nucleotide_material", "parent_genome_size": 12, "resulting_genome_size": 12}
	]
	var source_before: Array = events.duplicate(true)
	var summary: Dictionary = MutationDynamicsAnalyticsScript.summarize_event_log(events)
	_assert_true(events == source_before, "event analytics leave input ledger byte-for-byte semantically unchanged")
	_assert_true(int(summary["births"]) == 2, "analytics count only descendant births")
	_assert_true(int(summary["mutations"]) == 2, "analytics count only realized DNA changes as mutations")
	_assert_true(int(summary["mutation_attempts"]) == 3, "analytics preserve all sampled mutation attempts for auditability")
	_assert_true(int(summary["blocked_mutation_attempts"]) == 1, "physically blocked structural event is separated from realized mutation count")
	_assert_close(float(summary["realized_mutations_per_birth"]["mean"]), 1.0, 1e-12, "realized mutations per birth include zero-mutation sisters")
	_assert_close(float(summary["attempted_mutations_per_birth"]["mean"]), 1.5, 1e-12, "attempted mutation burden retains blocked structural draws")
	_assert_true(int(summary["mutation_types"]["protein_signature_bit_flip"]) == 1 and int(summary["mutation_types"]["gene_duplication"]) == 1, "realized mutation-type ledger excludes blocked attempts")
	_assert_true(int(summary["blocked_mutation_types"]["gene_duplication_blocked"]) == 1, "blocked structural attempt remains inspectable by exact type")
	_assert_close(float(summary["expected_point_error_rate_per_gene"]["mean"]), 0.002, 1e-12, "births inherit the replication-derived expected point-error context")
	_assert_true(int(summary["births_with_replication_profile"]) == 2, "both daughters resolve to their maternal replication profile")

func _test_functional_identity_separates_neutral_history() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var ancestral_function: String = MutationDynamicsAnalyticsScript.functional_key(ancestor)

	var neutral = ancestor.deep_copy()
	neutral.get_gene_by_locus(1).neutral_marker += 1
	_assert_true(neutral.canonical_key() != ancestor.canonical_key(), "neutral marker changes exact genetic identity")
	_assert_true(MutationDynamicsAnalyticsScript.functional_key(neutral) == ancestral_function, "neutral marker does not manufacture functional diversity")

	var reordered = ancestor.deep_copy()
	reordered.genes.reverse()
	_assert_true(reordered.canonical_key() != ancestor.canonical_key(), "gene-order rearrangement changes exact genome history")
	_assert_true(MutationDynamicsAnalyticsScript.functional_key(reordered) == ancestral_function, "gene order is not counted as functional diversity while order has no execution semantics")

	var coding = ancestor.deep_copy()
	coding.get_gene_by_locus(1).protein_signature ^= 1
	_assert_true(MutationDynamicsAnalyticsScript.functional_key(coding) != ancestral_function, "coding-sequence change produces distinct functional identity")

	var duplicated = ancestor.deep_copy()
	var mutator = MutationEngineScript.new()
	mutator.apply_structural_mutation(duplicated, "gene_duplication", DeterministicRngScript.new(7721))
	_assert_true(MutationDynamicsAnalyticsScript.functional_key(duplicated) != ancestral_function, "extra expressed locus changes functional inventory even when duplicated sequence is identical")

func _test_ancestral_departure_proxy_is_structural_not_fitness() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var baseline: Dictionary = MutationDynamicsAnalyticsScript.genome_departure_from_ancestor(ancestor)
	_assert_true(int(baseline["departure_count"]) == 0, "ancestral genome has no departure from itself")

	var changed = ancestor.deep_copy()
	changed.get_gene_by_locus(1).protein_signature ^= 1
	changed.get_gene_by_locus(2).promoter_copy_number = 2
	var mutator = MutationEngineScript.new()
	var duplication: Dictionary = mutator.apply_structural_mutation(changed, "gene_duplication", DeterministicRngScript.new(99117))
	var departure: Dictionary = MutationDynamicsAnalyticsScript.genome_departure_from_ancestor(changed)
	_assert_true(int(departure["departure_count"]) >= 3, "coding, cis and genome-expansion changes all contribute to molecular-departure proxy")
	_assert_true(departure["coding_signature_changes"].has(1), "proxy reports altered ancestral coding locus")
	_assert_true(departure["cis_copy_changes"].has(2), "proxy reports altered ancestral cis copy number")
	_assert_true(departure["extra_loci"].has(int(duplication["new_locus_id"])), "proxy reports duplicated non-ancestral locus")
	_assert_true(String(departure["note"]).contains("no beneficial/deleterious"), "departure proxy explicitly refuses a fitness interpretation")

func _test_environment_shift_association_tracks_mutator_abundance_descriptively() -> void:
	var before: Dictionary = {
		"tick": 100,
		"lineages": {
			"stable": {"count": 9, "mean_point_error_rate_per_gene": 0.001},
			"higher_error": {"count": 1, "mean_point_error_rate_per_gene": 0.002}
		}
	}
	var after: Dictionary = {
		"tick": 500,
		"lineages": {
			"stable": {"count": 1, "mean_point_error_rate_per_gene": 0.001},
			"higher_error": {"count": 9, "mean_point_error_rate_per_gene": 0.002}
		}
	}
	var association: Dictionary = MutationDynamicsAnalyticsScript.environment_shift_association(before, after)
	_assert_true(float(association["weighted_point_error_rate_change"]) > 0.0, "weighted mutation propensity rises when higher-error lineage expands after shift")
	_assert_close(float(association["above_baseline_error_abundance_before"]), 0.1, 1e-12, "pre-shift higher-error lineage abundance is measured")
	_assert_close(float(association["above_baseline_error_abundance_after"]), 0.9, 1e-12, "post-shift higher-error lineage expansion is measured")
	_assert_true(String(association["interpretation"]).contains("descriptive"), "environment-shift metric explicitly avoids causal overclaim")
	_assert_true(String(association["interpretation"]).contains("paired fork"), "analytics direct causal claims to M9 counterfactual machinery")

func _test_runner_detectors_use_m10_mechanistic_metrics() -> void:
	var trajectory: Array = [
		{
			"population": 2, "genotype_count": 1, "mutation_events": 0, "resources": {},
			"mutation_dynamics": {
				"point_error_rate_per_gene": {"count": 2, "mean": 0.001},
				"gene_count": {"count": 2, "mean": 12.0},
				"replication_units": {"count": 2, "mean": 15.0}
			}
		},
		{
			"population": 3, "genotype_count": 2, "mutation_events": 1, "resources": {},
			"mutation_dynamics": {
				"point_error_rate_per_gene": {"count": 3, "mean": 0.0015},
				"gene_count": {"count": 3, "mean": 13.0},
				"replication_units": {"count": 3, "mean": 16.25}
			}
		}
	]
	var run: Dictionary = {
		"termination_reason": "horizon",
		"environment": {"mode": "constant"},
		"trajectory": trajectory
	}
	var before: Dictionary = run.duplicate(true)
	var detected: Dictionary = ExperimentAnalyticsScript.detect_candidates(run)
	_assert_true(run == before, "M10-integrated experiment detectors remain observational")
	_assert_true(bool(detected["genome_expansion_or_reduction"]["evidence_available"]), "genome-size detector now has sampled M10 evidence")
	_assert_true(bool(detected["genome_expansion_or_reduction"]["candidate"]), "12-to-13 mean gene expansion is detected")
	_assert_close(float(detected["genome_expansion_or_reduction"]["gene_count_change"]), 1.0, 1e-12, "genome detector reports quantitative mean gene-count change")
	_assert_true(bool(detected["mutation_rate_shift"]["evidence_available"]), "mutation-rate detector now uses mechanistic replication-derived evidence")
	_assert_true(bool(detected["mutation_rate_shift"]["candidate"]), "change in expected copy-error propensity is detected")
	_assert_close(float(detected["mutation_rate_shift"]["change"]), 0.0005, 1e-12, "mutation-rate detector reports expected-rate change")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
