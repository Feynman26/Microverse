extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const DNAReplicationScript = preload("res://src/genetics/dna_replication.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const SnapshotCodecScript = preload("res://src/experiments/snapshot_codec.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_repair_is_sequence_derived_and_not_ancestrally_free()
	_test_repair_lowers_bounded_replication_error()
	_test_repair_has_explicit_atp_cost_and_can_slow_copying()
	_test_genome_size_has_material_and_time_cost()
	_test_division_waits_for_completed_genome_copy()
	_test_structural_operators_are_valid_and_deterministic()
	_test_architecture_change_does_not_teleport_molecules()
	_test_replication_state_survives_exact_snapshot()

	if failures == 0:
		print("PASS: %d M10 replication/fidelity tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 replication/fidelity tests failed" % [failures, tests_run])
		quit(1)

func _test_repair_is_sequence_derived_and_not_ancestrally_free() -> void:
	var config = SimConfigScript.new()
	var ancestor = _replicating_cell(1, GenomeScript.create_ancestor(), config)
	var ancestral_activity: float = DNAReplicationScript.repair_activity(ancestor.expression_state, ancestor.volume, config)
	_assert_close(ancestral_activity, 0.0, 1e-12, "canonical ancestor has no free DNA-repair activity")

	var ancestral_locus = ancestor.genome.get_gene_by_locus(12)
	_assert_true(ancestral_locus != null, "controlled repair accessibility uses canonical ancestral locus 12")
	_assert_true(DNAReplicationScript.hamming_distance(int(ancestral_locus.protein_signature), DNAReplicationScript.REPAIR_TARGET_SIGNATURE) == 5, "repair target begins exactly one coding mutation outside active radius")
	var accessible_signature: int = int(ancestral_locus.protein_signature) ^ (1 << 2)
	_assert_true(DNAReplicationScript.hamming_distance(accessible_signature, DNAReplicationScript.REPAIR_TARGET_SIGNATURE) == 4, "one ordinary bit flip enters generic repair-affinity radius")
	ancestor.expression_state[12]["protein"][accessible_signature] = float(config.expression_reference_protein_count) * ancestor.volume
	var mutant_activity: float = DNAReplicationScript.repair_activity(ancestor.expression_state, ancestor.volume, config)
	_assert_true(mutant_activity > 0.0, "repair activity appears only when a compatible realized protein cohort physically exists")

func _test_repair_lowers_bounded_replication_error() -> void:
	var config = SimConfigScript.new()
	var baseline = _replicating_cell(1, GenomeScript.create_ancestor(), config)
	var repaired = _replicating_cell(2, GenomeScript.create_ancestor(), config)
	_install_exact_repair_cohort(repaired, config, 1.0)
	_complete_replication_record(baseline, 0.0)
	_complete_replication_record(repaired, 1.0)
	var baseline_profile: Dictionary = DNAReplicationScript.mutation_profile(baseline, config)
	var repaired_profile: Dictionary = DNAReplicationScript.mutation_profile(repaired, config)
	_assert_true(float(repaired_profile["point_error_rate_per_gene"]) < float(baseline_profile["point_error_rate_per_gene"]), "greater functional repair investment lowers point-copy error probability")
	_assert_true(float(repaired_profile["structural_error_rate_per_genome"]) < float(baseline_profile["structural_error_rate_per_genome"]), "same repair chemistry lowers structural-copy error probability")
	_assert_true(float(repaired_profile["point_error_rate_per_gene"]) >= float(config.minimum_point_error_rate_per_gene), "repair cannot drive point error below configured chemical floor")
	_assert_true(float(repaired_profile["point_error_rate_per_gene"]) <= float(config.baseline_point_error_rate_per_gene), "point error remains bounded by unrepaired replication baseline")
	_assert_true(float(repaired_profile["structural_error_rate_per_genome"]) >= float(config.minimum_structural_error_rate_per_genome), "structural error cannot become negative or cross its chemical floor")

func _test_repair_has_explicit_atp_cost_and_can_slow_copying() -> void:
	var config = SimConfigScript.new()
	var baseline = _replicating_cell(1, GenomeScript.create_ancestor(), config)
	var repaired = _replicating_cell(2, GenomeScript.create_ancestor(), config)
	_install_exact_repair_cohort(repaired, config, 1.0)
	baseline.set_pool("ATP", 100.0)
	baseline.set_pool("NUC", 100.0)
	repaired.set_pool("ATP", 100.0)
	repaired.set_pool("NUC", 100.0)
	var plain_summary: Dictionary = DNAReplicationScript.step(baseline, 1.0, config)
	var repair_summary: Dictionary = DNAReplicationScript.step(repaired, 1.0, config)
	_assert_close(float(repair_summary["copied_this_tick"]), float(plain_summary["copied_this_tick"]), 1e-12, "abundant resources allow repaired and unrepaired forks to copy equal DNA amount")
	_assert_true(float(repair_summary["atp_spent_this_tick"]) > float(plain_summary["atp_spent_this_tick"]), "functional repair pays explicit extra ATP per copied gene")

	var starved_plain = _replicating_cell(3, GenomeScript.create_ancestor(), config)
	var starved_repair = _replicating_cell(4, GenomeScript.create_ancestor(), config)
	_install_exact_repair_cohort(starved_repair, config, 1.0)
	starved_plain.set_pool("ATP", 0.021)
	starved_repair.set_pool("ATP", 0.021)
	starved_plain.set_pool("NUC", 100.0)
	starved_repair.set_pool("NUC", 100.0)
	var starved_plain_summary: Dictionary = DNAReplicationScript.step(starved_plain, 1.0, config)
	var starved_repair_summary: Dictionary = DNAReplicationScript.step(starved_repair, 1.0, config)
	_assert_true(float(starved_repair_summary["copied_this_tick"]) < float(starved_plain_summary["copied_this_tick"]), "under equal ATP scarcity high-fidelity repair slows physical genome copying")
	_assert_true(starved_repair.pool("ATP") >= -1e-12, "repair cost cannot overdraw ATP")

func _test_genome_size_has_material_and_time_cost() -> void:
	var config = SimConfigScript.new()
	var small_genome = GenomeScript.create_ancestor()
	var large_genome = _double_genome(small_genome)
	var small = _replicating_cell(1, small_genome, config)
	var large = _replicating_cell(2, large_genome, config)
	for cell in [small, large]:
		cell.set_pool("ATP", 1000.0)
		cell.set_pool("NUC", 1000.0)
	var small_steps: int = _run_replication_to_completion(small, config)
	var large_steps: int = _run_replication_to_completion(large, config)
	_assert_true(small.genome.gene_count() == 12 and large.genome.gene_count() == 24, "controlled genome-size pair differs exactly twofold in copied gene count")
	_assert_true(large_steps == 2 * small_steps, "twice-larger genome requires twice the replication time at equal copying machinery")
	_assert_close(large.replication_atp_spent, 2.0 * small.replication_atp_spent, 1e-10, "twice-larger genome consumes twice the baseline replication ATP")
	_assert_close(large.replication_nuc_spent, 2.0 * small.replication_nuc_spent, 1e-10, "twice-larger genome consumes twice the nucleotide replication material")

func _test_division_waits_for_completed_genome_copy() -> void:
	var config = SimConfigScript.new()
	var cell = _replicating_cell(1, GenomeScript.create_ancestor(), config)
	cell.set_pool("BIO", float(config.division_volume) * float(config.biomass_units_per_volume))
	cell.volume = float(config.division_volume)
	cell.set_pool("ATP", 10.0)
	cell.replication_progress = 0.99
	cell.replication_gene_equivalents_copied = 0.99 * float(cell.genome.gene_count())
	_assert_true(not cell.ready_to_divide(config), "division-sized ATP-rich cell still waits for complete genome copy")
	cell.replication_progress = 1.0
	cell.replication_gene_equivalents_copied = float(cell.genome.gene_count())
	_assert_true(cell.ready_to_divide(config), "completed physical genome copy releases ordinary size/ATP division gate")

func _test_structural_operators_are_valid_and_deterministic() -> void:
	var mutator = MutationEngineScript.new()
	for operation in MutationEngineScript.STRUCTURAL_OPERATIONS:
		var first = GenomeScript.create_ancestor()
		var second = GenomeScript.create_ancestor()
		var event_a: Dictionary = mutator.apply_structural_mutation(first, String(operation), DeterministicRngScript.new(33001))
		var event_b: Dictionary = mutator.apply_structural_mutation(second, String(operation), DeterministicRngScript.new(33001))
		first.validate()
		second.validate()
		_assert_true(first.canonical_key() == second.canonical_key(), "%s yields deterministic genome representation for same state/seed" % operation)
		_assert_true(event_a == event_b, "%s yields deterministic inspectable event metadata" % operation)
		match String(operation):
			"gene_duplication": _assert_true(first.gene_count() == 13, "gene duplication adds exactly one parseable locus")
			"gene_deletion": _assert_true(first.gene_count() == 11, "gene deletion removes exactly one locus without invalidating genome")
			"segment_inversion", "local_rearrangement": _assert_true(first.gene_count() == 12, "%s preserves gene count while changing order" % operation)

func _test_architecture_change_does_not_teleport_molecules() -> void:
	var config = SimConfigScript.new()
	var mutator = MutationEngineScript.new()
	var duplicated = _replicating_cell(1, GenomeScript.create_ancestor(), config)
	var before_loci: Dictionary = duplicated.expression_state.duplicate(true)
	var duplication: Dictionary = mutator.apply_structural_mutation(duplicated.genome, "gene_duplication", DeterministicRngScript.new(8801))
	var new_locus: int = int(duplication["new_locus_id"])
	ExpressionSystemScript.reconcile_state_with_genome(duplicated.expression_state, duplicated.genome)
	_assert_true(duplicated.expression_state.has(new_locus), "newly duplicated DNA receives a molecular container")
	_assert_close(_molecular_locus_total(duplicated.expression_state[new_locus]), 0.0, 1e-12, "gene duplication does not instantly manufacture transcript or protein")
	for old_locus_variant in before_loci.keys():
		var old_locus: int = int(old_locus_variant)
		_assert_close(_molecular_locus_total(duplicated.expression_state[old_locus]), _molecular_locus_total(before_loci[old_locus]), 1e-12, "duplication leaves pre-existing locus-%d molecules untouched" % old_locus)

	var deleted = _replicating_cell(2, GenomeScript.create_ancestor(), config)
	var deletion: Dictionary = mutator.apply_structural_mutation(deleted.genome, "gene_deletion", DeterministicRngScript.new(9901))
	var deleted_locus: int = int(deletion["locus_id"])
	var inherited_amount: float = _molecular_locus_total(deleted.expression_state[deleted_locus])
	ExpressionSystemScript.reconcile_state_with_genome(deleted.expression_state, deleted.genome)
	_assert_true(deleted.expression_state.has(deleted_locus), "deleted DNA does not erase already-realized molecular cohort")
	_assert_close(_molecular_locus_total(deleted.expression_state[deleted_locus]), inherited_amount, 1e-12, "deleted-locus molecules persist until ordinary decay rather than teleporting away")

func _test_replication_state_survives_exact_snapshot() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	config.mutation_enabled = false
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor()
	cell.set_pool("BIO", 1.5)
	cell.volume = 1.5
	cell.set_pool("ATP", 50.0)
	cell.set_pool("NUC", 50.0)
	DNAReplicationScript.step(cell, 2.5, config)
	var before_checksum: float = sim.checksum()
	var snapshot: Dictionary = SnapshotCodecScript.capture(sim)
	var restored = SnapshotCodecScript.restore(snapshot)
	var restored_cell = restored.cells[0]
	_assert_close(restored_cell.replication_progress, cell.replication_progress, 1e-12, "snapshot restores exact partial DNA-copy progress")
	_assert_close(restored_cell.replication_atp_spent, cell.replication_atp_spent, 1e-12, "snapshot restores exact accumulated replication ATP cost")
	_assert_close(restored_cell.replication_nuc_spent, cell.replication_nuc_spent, 1e-12, "snapshot restores exact accumulated replication nucleotide cost")
	_assert_close(restored.checksum(), before_checksum, 1e-9, "snapshot/fork checksum includes authoritative M10 replication state exactly")

func _replicating_cell(cell_id: int, genome, config):
	var cell = CellStateScript.new(cell_id, -1, 0, 0, Vector2(3.0, 3.0), 1.5)
	cell.genome = genome.deep_copy()
	cell.initialize_molecular_state(config)
	cell.set_pool("BIO", 1.5 * float(config.biomass_units_per_volume))
	cell.volume = 1.5
	cell.replication_progress = 0.0
	cell.replication_gene_equivalents_copied = 0.0
	cell.replication_atp_spent = 0.0
	cell.replication_nuc_spent = 0.0
	cell.replication_repair_activity_integral = 0.0
	return cell

func _install_exact_repair_cohort(cell, config, concentration_scale: float) -> void:
	var amount: float = float(config.expression_reference_protein_count) * cell.volume * concentration_scale
	cell.expression_state[12]["protein"][DNAReplicationScript.REPAIR_TARGET_SIGNATURE] = amount

func _complete_replication_record(cell, mean_repair: float) -> void:
	var copied: float = float(cell.genome.gene_count())
	cell.replication_gene_equivalents_copied = copied
	cell.replication_progress = 1.0
	cell.replication_repair_activity_integral = copied * mean_repair

func _double_genome(source):
	var genes: Array = []
	for gene in source.genes:
		genes.append(gene.deep_copy())
	for gene in source.genes:
		var copy = gene.deep_copy()
		copy.locus_id = int(gene.locus_id) + 1000
		genes.append(copy)
	return GenomeScript.new(genes)

func _run_replication_to_completion(cell, config) -> int:
	var steps: int = 0
	while not DNAReplicationScript.replication_complete(cell):
		DNAReplicationScript.step(cell, 1.0, config)
		steps += 1
		assert(steps < 1000)
	return steps

func _molecular_locus_total(locus_state: Dictionary) -> float:
	var result: float = 0.0
	for species_name in ["mrna", "protein"]:
		for amount in locus_state[species_name].values():
			result += float(amount)
	return result

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])