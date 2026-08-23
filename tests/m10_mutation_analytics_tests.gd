extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const MutationDynamicsAnalyticsScript = preload("res://src/experiments/mutation_dynamics_analytics.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_population_sampling_is_observational()
	_test_event_ledger_reconstructs_realized_mutations_per_birth()
	_test_ancestral_departure_proxy_is_structural_not_fitness()
	_test_environment_shift_association_tracks_mutator_abundance_descriptively()
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
		{"kind": "mutation", "cell_id": 2, "mutation_type": "gene_duplication", "parent_genome_size": 12, "resulting_genome_size": 13}
	]
	var source_before: Array = events.duplicate(true)
	var summary: Dictionary = MutationDynamicsAnalyticsScript.summarize_event_log(events)
	_assert_true(events == source_before, "event analytics leave input ledger byte-for-byte semantically unchanged")
	_assert_true(int(summary["births"]) == 2, "analytics count only descendant births")
	_assert_true(int(summary["mutations"]) == 2, "analytics count realized mutation events")
	_assert_close(float(summary["realized_mutations_per_birth"]["mean"]), 1.0, 1e-12, "realized mutations per birth include zero-mutation sisters")
	_assert_true(int(summary["mutation_types"]["protein_signature_bit_flip"]) == 1 and int(summary["mutation_types"]["gene_duplication"]) == 1, "mutation-type ledger remains inspectable without value labels")
	_assert_close(float(summary["expected_point_error_rate_per_gene"]["mean"]), 0.002, 1e-12, "births inherit the replication-derived expected point-error context")
	_assert_true(int(summary["births_with_replication_profile"]) == 2, "both daughters resolve to their maternal replication profile")

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

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
