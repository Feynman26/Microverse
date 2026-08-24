extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const DNAReplicationScript = preload("res://src/genetics/dna_replication.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	var costly = SimConfigScript.new()
	var free = SimConfigScript.new()
	free.dna_repair_atp_cost_per_gene_activity = 0.0

	var costly_plain = _cell(1, costly, false)
	var costly_repair = _cell(2, costly, true)
	var free_plain = _cell(3, free, false)
	var free_repair = _cell(4, free, true)
	for cell in [costly_plain, costly_repair, free_plain, free_repair]:
		cell.set_pool("ATP", 0.021)
		cell.set_pool("NUC", 100.0)

	var costly_plain_summary: Dictionary = DNAReplicationScript.step(costly_plain, 1.0, costly)
	var costly_repair_summary: Dictionary = DNAReplicationScript.step(costly_repair, 1.0, costly)
	var free_plain_summary: Dictionary = DNAReplicationScript.step(free_plain, 1.0, free)
	var free_repair_summary: Dictionary = DNAReplicationScript.step(free_repair, 1.0, free)

	_assert_true(float(costly_repair_summary["copied_this_tick"]) < float(costly_plain_summary["copied_this_tick"]), "with physical repair ATP cost, higher fidelity can slow DNA copying under equal energy scarcity")
	_assert_close(float(free_repair_summary["copied_this_tick"]), float(free_plain_summary["copied_this_tick"]), 1e-12, "removing repair ATP cost makes repaired and unrepaired replication equally fast")

	_complete_record(free_plain, 0.0)
	_complete_record(free_repair, 1.0)
	var plain_profile: Dictionary = DNAReplicationScript.mutation_profile(free_plain, free)
	var repair_profile: Dictionary = DNAReplicationScript.mutation_profile(free_repair, free)
	_assert_true(float(repair_profile["point_error_rate_per_gene"]) < float(plain_profile["point_error_rate_per_gene"]), "zero-cost control still receives higher fidelity")
	_assert_true(float(free_repair_summary["atp_spent_this_tick"]) == float(free_plain_summary["atp_spent_this_tick"]), "zero-cost control obtains fidelity without additional replication ATP burden")
	_assert_true(float(costly_repair_summary["atp_spent_this_tick"]) >= float(costly_plain_summary["atp_spent_this_tick"]), "production control retains a real energetic price for repair investment")

	if failures == 0:
		print("PASS: %d M10 repair-cost control tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 repair-cost control tests failed" % [failures, tests_run])
		quit(1)

func _cell(cell_id: int, config, repaired: bool):
	var cell = CellStateScript.new(cell_id, -1, 0, 0, Vector2(3.0, 3.0), 1.5)
	cell.genome = GenomeScript.create_ancestor()
	cell.initialize_molecular_state(config)
	cell.set_pool("BIO", 1.5 * float(config.biomass_units_per_volume))
	cell.volume = 1.5
	cell.replication_progress = 0.0
	cell.replication_gene_equivalents_copied = 0.0
	cell.replication_atp_spent = 0.0
	cell.replication_nuc_spent = 0.0
	cell.replication_repair_activity_integral = 0.0
	if repaired:
		cell.expression_state[12]["protein"][DNAReplicationScript.REPAIR_TARGET_SIGNATURE] = float(config.expression_reference_protein_count) * cell.volume
	return cell

func _complete_record(cell, mean_repair: float) -> void:
	var copied: float = float(cell.genome.replication_unit_count())
	cell.replication_gene_equivalents_copied = copied
	cell.replication_progress = 1.0
	cell.replication_repair_activity_integral = copied * mean_repair

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
