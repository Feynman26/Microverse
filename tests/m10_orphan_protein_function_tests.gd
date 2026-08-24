extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_deleted_gene_does_not_switch_off_existing_enzyme()
	_test_genome_reordering_cannot_change_realized_catalysis()
	if failures == 0:
		print("PASS: %d M10 orphan-protein function tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 orphan-protein function tests failed" % [failures, tests_run])
		quit(1)

func _test_deleted_gene_does_not_switch_off_existing_enzyme() -> void:
	var config = SimConfigScript.new()
	var cell = _ancestor_cell(config)
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var reaction = ReactionCatalogScript.by_id(reactions, "R01")
	var before: float = CatalyticLandscapeScript.proteome_activity(cell.genome, cell.expression_state, reaction, config)
	_assert_true(before > 0.0, "ancestral realized proteome catalyses R01 before deletion")

	var deleted_locus: int = int(cell.genome.genes[0].locus_id)
	_assert_true(deleted_locus == 1, "controlled deletion fixture targets ancestral locus 1")
	cell.genome.genes.remove_at(0)
	cell.genome.validate()
	_assert_true(cell.expression_state.has(deleted_locus), "gene deletion leaves its existing molecular cohort physically present")
	var after_dna_deletion: float = CatalyticLandscapeScript.proteome_activity(cell.genome, cell.expression_state, reaction, config)
	_assert_close(after_dna_deletion, before, 0.0, "deleting DNA does not instantaneously switch off already-existing enzyme molecules")

	cell.expression_state[deleted_locus]["protein"].clear()
	var after_protein_removal: float = CatalyticLandscapeScript.proteome_activity(cell.genome, cell.expression_state, reaction, config)
	_assert_true(after_protein_removal < after_dna_deletion, "catalytic function declines only when the physical protein cohort is actually removed/decays")

func _test_genome_reordering_cannot_change_realized_catalysis() -> void:
	var config = SimConfigScript.new()
	var cell = _ancestor_cell(config)
	var reaction = ReactionCatalogScript.by_id(ReactionCatalogScript.create_m4_candidate(), "R03")
	var before: float = CatalyticLandscapeScript.proteome_activity(cell.genome, cell.expression_state, reaction, config)
	cell.genome.genes.reverse()
	cell.genome.validate()
	var after: float = CatalyticLandscapeScript.proteome_activity(cell.genome, cell.expression_state, reaction, config)
	_assert_close(after, before, 0.0, "genome array rearrangement cannot change chemistry through floating-point accumulation order")

func _ancestor_cell(config):
	var cell = CellStateScript.new(1, -1, 0, 0, Vector2(4.0, 4.0), config.ancestor_volume)
	cell.genome = GenomeScript.create_ancestor()
	cell.initialize_molecular_state(config)
	return cell

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
