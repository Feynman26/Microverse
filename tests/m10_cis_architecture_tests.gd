extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const SnapshotCodecScript = preload("res://src/experiments/snapshot_codec.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_single_copy_preserves_canon()
	_test_promoter_copy_changes_transcription_opportunity()
	_test_regulatory_copy_changes_binding_opportunity()
	_test_region_structure_changes_replication_units()
	_test_region_architecture_survives_snapshot()
	if failures == 0:
		print("PASS: %d M10 cis-architecture tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 cis-architecture tests failed" % [failures, tests_run])
		quit(1)

func _test_single_copy_preserves_canon() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var first = ancestor.get_gene_by_locus(1)
	_assert_true(first.canonical_key() == "1:6200:4951:101:4952", "single-copy cis architecture preserves historical gene canonical key")
	_assert_true(not first.canonical_key().contains("PC") and not first.canonical_key().contains("RC"), "default cis-copy metadata stays out of ancestral representation")
	var copied = ancestor.deep_copy()
	_assert_true(copied.canonical_key() == ancestor.canonical_key(), "deep copy preserves historical genome representation")
	_assert_true(copied.fingerprint() == ancestor.fingerprint(), "default cis representation preserves historical fingerprint")

func _test_promoter_copy_changes_transcription_opportunity() -> void:
	var mutator = MutationEngineScript.new()
	var duplicated = GenomeScript.create_ancestor()
	var before_units: float = duplicated.replication_unit_count()
	var event: Dictionary = mutator.apply_structural_mutation(duplicated, "promoter_region_duplication", DeterministicRngScript.new(44001))
	var gene = duplicated.get_gene_by_locus(int(event["locus_id"]))
	_assert_true(int(gene.promoter_copy_number) == 2, "promoter-region duplication creates a heritable second copy")
	var single_strength: float = float(gene.promoter_code) / float(GeneScript.PROMOTER_CODE_MAX)
	_assert_close(gene.promoter_strength(), 2.0 * single_strength, 1e-12, "two identical promoter copies create two transcription-initiation opportunities")
	_assert_close(duplicated.replication_unit_count() - before_units, GenomeScript.CIS_REGION_REPLICATION_UNITS, 1e-12, "duplicated promoter adds explicit DNA copying burden")
	_assert_true(gene.canonical_key().contains("PC2"), "promoter structural duplication is inherited in canonical genome")

	var deleted = GenomeScript.create_ancestor()
	var deletion: Dictionary = mutator.apply_structural_mutation(deleted, "promoter_region_deletion", DeterministicRngScript.new(44001))
	var deleted_gene = deleted.get_gene_by_locus(int(deletion["locus_id"]))
	_assert_true(int(deleted_gene.promoter_copy_number) == 0, "promoter-region deletion removes the cis promoter while retaining coding DNA")
	_assert_close(deleted_gene.promoter_strength(), 0.0, 1e-12, "zero promoter copies eliminate new transcription opportunity")

func _test_regulatory_copy_changes_binding_opportunity() -> void:
	var config = SimConfigScript.new()
	config.allostery_enabled = false
	var factor_one: float = _target_regulation_factor(_regulatory_fixture_genome(1), config)
	var factor_two: float = _target_regulation_factor(_regulatory_fixture_genome(2), config)
	var factor_zero: float = _target_regulation_factor(_regulatory_fixture_genome(0), config)
	_assert_true(factor_one > 1.0, "one compatible cis-regulatory site responds to an activator protein")
	_assert_true(factor_two > factor_one, "duplicating the same regulatory site strengthens occupancy")
	_assert_close(factor_zero, 1.0, 1e-12, "deleting all regulatory copies removes the cis edge while retaining basal promoter activity")

func _test_region_structure_changes_replication_units() -> void:
	var mutator = MutationEngineScript.new()
	for pair in [
		["promoter_region_duplication", 1.0],
		["promoter_region_deletion", -1.0],
		["regulatory_region_duplication", 1.0],
		["regulatory_region_deletion", -1.0]
	]:
		var genome = GenomeScript.create_ancestor()
		var before: float = genome.replication_unit_count()
		var event: Dictionary = mutator.apply_structural_mutation(genome, String(pair[0]), DeterministicRngScript.new(55002))
		var expected_delta: float = float(pair[1]) * GenomeScript.CIS_REGION_REPLICATION_UNITS
		_assert_close(genome.replication_unit_count() - before, expected_delta, 1e-12, "%s changes DNA-copying target by exact cis-region size" % pair[0])
		_assert_close(float(event["replication_units_delta"]), expected_delta, 1e-12, "%s event ledger reports exact genome-size delta" % pair[0])

func _test_region_architecture_survives_snapshot() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	config.mutation_enabled = false
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor()
	var mutator = MutationEngineScript.new()
	mutator.apply_structural_mutation(cell.genome, "promoter_region_duplication", DeterministicRngScript.new(66003))
	mutator.apply_structural_mutation(cell.genome, "regulatory_region_deletion", DeterministicRngScript.new(66004))
	var key_before: String = cell.genome.canonical_key()
	var units_before: float = cell.genome.replication_unit_count()
	var checksum_before: float = sim.checksum()
	var restored = SnapshotCodecScript.restore(SnapshotCodecScript.capture(sim))
	_assert_true(restored.cells[0].genome.canonical_key() == key_before, "snapshot preserves promoter/regulatory copy architecture exactly")
	_assert_close(restored.cells[0].genome.replication_unit_count(), units_before, 1e-12, "snapshot preserves cis-DNA contribution to copying target")
	_assert_close(restored.checksum(), checksum_before, 1e-9, "cis architecture participates in exact replay checksum")

func _regulatory_fixture_genome(copy_number: int):
	return GenomeScript.new([
		GeneScript.new(1, 6000, 0x1234, 1, 0xEEEE),
		GeneScript.new(2, 5000, 0x5678, 2, 0x1234, 1, copy_number)
	])

func _target_regulation_factor(genome, config) -> float:
	var state: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 1.0)
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(1.0, config)
	pools["ATP"] = 100.0
	pools["NUC"] = 100.0
	pools["AA"] = 100.0
	var summary: Dictionary = ExpressionSystemScript.step(
		state, genome, pools, config.tick_dt_min, DeterministicRngScript.new(77005), config, 1.0
	)
	return float(summary["regulation"][2])

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])