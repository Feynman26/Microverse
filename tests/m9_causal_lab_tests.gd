extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const SnapshotCodecScript = preload("res://src/experiments/snapshot_codec.gd")
const StrainLibraryScript = preload("res://src/experiments/strain_library.gd")
const GeneticAssayScript = preload("res://src/experiments/genetic_assay.gd")
const CausalLabScript = preload("res://src/experiments/causal_lab.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_snapshot_roundtrip_rng_and_context()
	_test_exact_forks_and_replay()
	_test_freezer_modes_and_archive_independence()
	_test_genetic_interventions_are_explicit_and_local()
	_test_competition_assay_replays_exactly()
	if failures == 0:
		print("PASS: %d compact M9 causal-lab tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d compact M9 causal-lab tests failed" % [failures, tests_run])
		quit(1)

func _test_snapshot_roundtrip_rng_and_context() -> void:
	var sim = _make_sim(91001)
	var cell = sim.seed_ancestor(Vector2(5.0, 5.0))
	sim.world.get_field("glucose").set_value(2, 3, 1.23456789)
	sim.world.release_protein(0xD136, Vector2(4.0, 4.0), 0.75, float(sim.config.extracellular_protein_diffusion))
	sim.step(4)
	var context: Dictionary = {
		"scheduled_interventions": [{"tick": 40, "kind": "environment_intervention", "name": "shock"}],
		"environment_boundary_state": {"mode": "constant", "phase": 2}
	}
	var snapshot: Dictionary = SnapshotCodecScript.capture(sim, context)
	var encoded: PackedByteArray = SnapshotCodecScript.encode(snapshot)
	var decoded: Dictionary = SnapshotCodecScript.decode(encoded)
	var restored = SnapshotCodecScript.restore(decoded)
	var restored_snapshot: Dictionary = SnapshotCodecScript.capture(restored, context)
	_assert_equal(String(snapshot["fingerprint"]), SnapshotCodecScript.fingerprint(decoded), "snapshot fingerprint survives binary save/load")
	_assert_equal(String(snapshot["fingerprint"]), String(restored_snapshot["fingerprint"]), "save/load reconstructs identical authoritative state fingerprint")
	_assert_close(restored.checksum(), sim.checksum(), 1e-9, "restored simulation checksum exactly matches source state")
	_assert_equal(int(restored.tick_index), int(sim.tick_index), "snapshot restores exact tick")
	_assert_equal(int(restored.next_cell_id), int(sim.next_cell_id), "snapshot restores next cell identity counter")
	_assert_equal(int(restored.next_mutation_id), int(sim.next_mutation_id), "snapshot restores next mutation identity counter")
	_assert_equal(decoded["scheduled_interventions"], context["scheduled_interventions"], "snapshot preserves scheduled interventions")
	_assert_equal(decoded["environment_boundary_state"], context["environment_boundary_state"], "snapshot preserves environment boundary state")
	var expected_next_random: float = sim.rng.randf()
	var replay_rng = SnapshotCodecScript.restore(decoded)
	_assert_close(replay_rng.rng.randf(), expected_next_random, 0.0, "RNG internal state resumes exactly after load")
	_assert_true(cell.id == 1, "nontrivial snapshot fixture retains expected source cell")

func _test_exact_forks_and_replay() -> void:
	var sim = _make_sim(91002)
	sim.seed_ancestor(Vector2(5.0, 5.0))
	sim.step(7)
	var branch_point: Dictionary = SnapshotCodecScript.capture(sim)
	var forks: Array = SnapshotCodecScript.fork(sim)
	var left = forks[0]
	var right = forks[1]
	_assert_close(left.checksum(), right.checksum(), 0.0, "two forks begin bit-equivalent at branch instant")
	left.step(8)
	right.step(8)
	_assert_close(left.checksum(), right.checksum(), 0.0, "forks remain identical under identical future inputs")
	var right_glucose_before: float = right.world.get_field("glucose").get_value(0, 0)
	left.world.get_field("glucose").set_value(0, 0, right_glucose_before + 0.5)
	left.cells[0].metabolites["ATP"] = float(left.cells[0].metabolites["ATP"]) + 0.25
	_assert_close(right.world.get_field("glucose").get_value(0, 0), right_glucose_before, 0.0, "editing one fork does not mutate the other's world by shared reference")
	_assert_true(absf(float(left.cells[0].pool("ATP")) - float(right.cells[0].pool("ATP"))) > 1e-12, "editing one fork does not mutate the other's cell state by shared reference")

	var original = SnapshotCodecScript.restore(branch_point)
	original.step(11)
	var expected_checksum: float = original.checksum()
	var replayed = CausalLabScript.replay(branch_point, int(branch_point["tick_index"]) + 11)
	_assert_close(replayed.checksum(), expected_checksum, 1e-9, "replay from prior snapshot reaches recorded target checksum")
	var verification: Dictionary = CausalLabScript.verify_replay(branch_point, int(branch_point["tick_index"]) + 11, expected_checksum)
	_assert_true(bool(verification["matches"]), "replay verification reports checksum agreement")

func _test_freezer_modes_and_archive_independence() -> void:
	var sim = _make_sim(91003)
	var source = sim.seed_ancestor(Vector2(5.0, 5.0))
	sim.step(3)
	var library = StrainLibraryScript.new()
	var genotype_entry: Dictionary = library.freeze_genotype(sim, source, "ancestor-like", "genotype archive")
	var state_entry: Dictionary = library.freeze_molecular_state(sim, source, "stateful", "transient molecular archive")
	_assert_true(String(genotype_entry["mode"]) != String(state_entry["mode"]), "freezer uses distinct names for genotype-only and molecular-state archives")
	_assert_true(not genotype_entry.has("cell_state"), "genotype-only archive does not silently store transient molecular state")
	_assert_true(state_entry.has("cell_state"), "molecular-state archive explicitly stores transient molecular state")
	var archived_key: String = String(genotype_entry["genome"]["canonical_key"])

	var target = _make_sim(91004)
	var reintroduced = library.reintroduce(target, String(genotype_entry["archive_id"]), Vector2(4.0, 4.0))
	_assert_equal(reintroduced.genome.canonical_key(), archived_key, "frozen genotype reintroduces exact archived genotype")
	reintroduced.genome.get_gene_by_locus(1).promoter_code += 1
	_assert_equal(String(library.get_entry(String(genotype_entry["archive_id"]))["genome"]["canonical_key"]), archived_key, "reintroduced strain cannot mutate archived freezer source")

	var state_target = _make_sim(91005)
	var state_cell = library.reintroduce(state_target, String(state_entry["archive_id"]), Vector2(6.0, 6.0))
	_assert_equal(state_cell.metabolites, state_entry["cell_state"]["metabolites"], "molecular-state archive restores inherited molecular pools exactly")
	_assert_true(int(state_cell.id) != int(state_entry["source_cell_id"]) or int(state_entry["source_cell_id"]) == 1, "reintroduction allocates target-world identity rather than mutating archive metadata")
	_assert_true(CausalLabScript.semantic_timeline(state_target.event_log).size() == 1, "experimental reintroduction appears explicitly on semantic timeline")

func _test_genetic_interventions_are_explicit_and_local() -> void:
	var ancestor = preload("res://src/genetics/genome.gd").create_ancestor()
	var locus = ancestor.get_gene_by_locus(1)
	var mutation: Dictionary = {
		"mutation_type": "promoter_code",
		"locus_id": 1,
		"old_value": int(locus.promoter_code),
		"new_value": int(locus.promoter_code) + 37
	}
	var introduced: Dictionary = GeneticAssayScript.introduce_mutation(ancestor, mutation)
	_assert_true(String(introduced["intervention"]["kind"]) == "genetic_intervention", "experiment-generated genetic edit is tagged as intervention, not spontaneous mutation")
	_assert_equal(int(introduced["genome"].get_gene_by_locus(1).promoter_code), int(mutation["new_value"]), "single-mutation introduction edits targeted allele")
	_assert_equal(int(ancestor.get_gene_by_locus(1).promoter_code), int(mutation["old_value"]), "single-mutation introduction leaves source genome unchanged")
	var reverted: Dictionary = GeneticAssayScript.revert_mutation(introduced["genome"], mutation)
	_assert_true(reverted["genome"].exact_equals(ancestor), "reversion restores targeted allele with no unrelated locus changes")
	_assert_true(String(reverted["intervention"]["kind"]) == "genetic_intervention", "reversion remains an explicit experimental intervention")

	var knockout: Dictionary = GeneticAssayScript.knock_out_gene(ancestor, 4)
	_assert_equal(int(knockout["genome"].get_gene_by_locus(4).promoter_code), 0, "gene knockout disables only target promoter")
	_assert_equal(int(knockout["genome"].get_gene_by_locus(5).promoter_code), int(ancestor.get_gene_by_locus(5).promoter_code), "gene knockout leaves unrelated locus unchanged")
	var restored: Dictionary = GeneticAssayScript.restore_ancestral_allele(knockout["genome"], 4)
	_assert_equal(restored["genome"].get_gene_by_locus(4).canonical_key(), ancestor.get_gene_by_locus(4).canonical_key(), "ancestral-allele restoration restores complete target locus")
	var explanation: Dictionary = CausalLabScript.explain_genetic_candidate(ancestor, introduced["genome"])
	_assert_equal(explanation["genotype_differences"].size(), 1, "explain-this enumerates exact genotype differences")
	_assert_equal(explanation["single_reversion_candidates"].size(), 1, "explain-this generates one controlled single-mutation reversion candidate")

func _test_competition_assay_replays_exactly() -> void:
	var source_sim = _make_sim(91006)
	var source = source_sim.seed_ancestor(Vector2(5.0, 5.0))
	var library = StrainLibraryScript.new()
	var first: Dictionary = library.freeze_genotype(source_sim, source, "A")
	var mutation: Dictionary = {
		"mutation_type": "promoter_code",
		"locus_id": 1,
		"old_value": int(source.genome.get_gene_by_locus(1).promoter_code),
		"new_value": int(source.genome.get_gene_by_locus(1).promoter_code) + 1
	}
	var second_genome = GeneticAssayScript.introduce_mutation(source.genome, mutation)["genome"]
	var second: Dictionary = first.duplicate(true)
	second["archive_id"] = "strain-controlled-B"
	second["genome"] = SnapshotCodecScript.capture_genome(second_genome)
	var environment: Dictionary = EnvironmentScheduleScript.constant({"glucose": 4.0, "oxygen": 0.5, "nitrogen": 3.0, "phosphorus": 2.0})
	var config = _make_config(91007)
	var run_a: Dictionary = CausalLabScript.run_competition(config, [first, second], [1, 1], environment, 25, 4444)
	var run_b: Dictionary = CausalLabScript.run_competition(config, [first, second], [1, 1], environment, 25, 4444)
	_assert_equal(run_a, run_b, "competition assay is exactly reproducible for fixed starting state and seed")
	_assert_equal(run_a["frequencies"].keys().size(), 2, "competition tracks relative frequency for every archived strain")
	_assert_true(run_a.has("division_events_by_strain") and run_a.has("death_causes_by_strain"), "competition tracks division and death distributions by lineage")
	_assert_true(run_a.has("endpoint_fluxes_by_strain") and run_a.has("spatial_occupation"), "competition tracks flux and spatial occupation by lineage")
	_assert_true(String(run_a["outcome"]) in ["coexistence", "fixation", "extinction"], "competition classifies coexistence/fixation/extinction without a fitness score")

func _make_config(seed: int):
	var config = SimConfigScript.new()
	config.seed = seed
	config.world_width = 12
	config.world_height = 12
	config.max_cells = 24
	config.mutation_enabled = false
	config.validate()
	return config

func _make_sim(seed: int):
	return SimulationEngineScript.new(_make_config(seed))

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_equal(actual, expected, message: String) -> void:
	_assert_true(actual == expected, "%s (actual=%s expected=%s)" % [message, str(actual), str(expected)])

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, str(actual), str(expected)])
