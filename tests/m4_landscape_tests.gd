extends SceneTree

const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const CatalyticLandscapeScript = preload("res://src/chemistry/catalytic_landscape.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	ReactionCatalogScript.validate_unique(reactions)
	_test_catalog_scope(reactions)
	_test_all_reactions_close_structural_balance(reactions)
	_test_ancestral_accessibility_pattern(reactions)
	_test_dormant_recovery_is_one_mutation_from_activity(reactions)
	_test_single_bit_neighbors_preserve_local_continuity(reactions)
	_test_random_landscape_density_is_sparse_but_connected(reactions)

	if failures == 0:
		print("PASS: %d M4 catalytic-landscape tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M4 catalytic-landscape tests failed" % [failures, tests_run])
		quit(1)

func _test_catalog_scope(reactions: Array) -> void:
	_assert_true(MetaboliteCatalogScript.ids().size() == 18, "candidate chemistry defines eighteen primitive metabolites")
	_assert_true(reactions.size() == 11, "candidate chemistry defines eleven catalytic reactions")
	_assert_true(ReactionCatalogScript.by_id(reactions, "R01") != null, "reaction IDs are addressable")

func _test_all_reactions_close_structural_balance(reactions: Array) -> void:
	var balanced: bool = true
	for reaction in reactions:
		if not reaction.is_structurally_balanced():
			balanced = false
			break
	_assert_true(balanced, "all M4 candidate reactions conserve digital C/N/P structural units")

func _test_ancestral_accessibility_pattern(reactions: Array) -> void:
	var ancestor = GenomeScript.create_ancestor()
	var exact_expectations: Dictionary = {
		"R01": 1, "R02": 2, "R03": 3,
		"R06": 5, "R07": 6, "R08": 7, "R09": 8
	}
	var exact_ok: bool = true
	for reaction_id in exact_expectations.keys():
		var reaction = ReactionCatalogScript.by_id(reactions, String(reaction_id))
		var locus: int = int(exact_expectations[reaction_id])
		var gene = ancestor.get_gene_by_locus(locus)
		if CatalyticLandscapeScript.hamming_distance(int(gene.protein_signature), int(reaction.signature)) != 0:
			exact_ok = false
			break
	_assert_true(exact_ok, "ancestral core reactions have exact catalytic matches")

	var fermentation = ReactionCatalogScript.by_id(reactions, "R04")
	var overflow = ReactionCatalogScript.by_id(reactions, "R10")
	var gene4 = ancestor.get_gene_by_locus(4)
	var gene9 = ancestor.get_gene_by_locus(9)
	_assert_true(CatalyticLandscapeScript.hamming_distance(gene4.protein_signature, fermentation.signature) == 3, "ancestral fermentation is a weak but active side route")
	_assert_true(CatalyticLandscapeScript.hamming_distance(gene9.protein_signature, overflow.signature) == 3, "ancestral W2 overflow is a weak but active side route")

	var recovery1 = ReactionCatalogScript.by_id(reactions, "R05")
	var recovery2 = ReactionCatalogScript.by_id(reactions, "R11")
	_assert_true(CatalyticLandscapeScript.hamming_distance(gene4.protein_signature, recovery1.signature) == 5, "W1 recovery begins just outside active radius")
	_assert_true(CatalyticLandscapeScript.hamming_distance(gene9.protein_signature, recovery2.signature) == 5, "W2 recovery begins just outside active radius")
	_assert_close(CatalyticLandscapeScript.affinity(gene4.protein_signature, recovery1.signature), 0.0, 1e-12, "dormant W1 recovery has zero ancestral activity")
	_assert_close(CatalyticLandscapeScript.affinity(gene9.protein_signature, recovery2.signature), 0.0, 1e-12, "dormant W2 recovery has zero ancestral activity")

func _test_dormant_recovery_is_one_mutation_from_activity(reactions: Array) -> void:
	var ancestor = GenomeScript.create_ancestor()
	var pairs: Array = [
		[ancestor.get_gene_by_locus(4).protein_signature, ReactionCatalogScript.by_id(reactions, "R05").signature],
		[ancestor.get_gene_by_locus(9).protein_signature, ReactionCatalogScript.by_id(reactions, "R11").signature]
	]
	var accessible: bool = true
	for pair in pairs:
		var protein_signature: int = int(pair[0])
		var reaction_signature: int = int(pair[1])
		var difference: int = (protein_signature ^ reaction_signature) & 0xFFFF
		var chosen_bit: int = -1
		for bit_index in range(16):
			if (difference & (1 << bit_index)) != 0:
				chosen_bit = bit_index
				break
		if chosen_bit < 0:
			accessible = false
			continue
		var mutant_signature: int = protein_signature ^ (1 << chosen_bit)
		if CatalyticLandscapeScript.hamming_distance(mutant_signature, reaction_signature) != 4:
			accessible = false
		if CatalyticLandscapeScript.affinity(mutant_signature, reaction_signature) <= 0.0:
			accessible = false
	_assert_true(accessible, "each dormant waste-recovery route is reachable by one coding bit flip")

func _test_single_bit_neighbors_preserve_local_continuity(reactions: Array) -> void:
	var reaction = ReactionCatalogScript.by_id(reactions, "R01")
	var exact_signature: int = int(reaction.signature)
	var expected_neighbor_affinity: float = exp(-CatalyticLandscapeScript.DISTANCE_DECAY)
	var all_neighbors_active: bool = true
	for bit_index in range(16):
		var mutant: int = exact_signature ^ (1 << bit_index)
		var affinity: float = CatalyticLandscapeScript.affinity(mutant, exact_signature)
		if affinity <= 0.0 or absf(affinity - expected_neighbor_affinity) > 1e-12:
			all_neighbors_active = false
			break
	_assert_true(all_neighbors_active, "all one-bit neighbors of an exact enzyme retain partial catalytic activity")

func _test_random_landscape_density_is_sparse_but_connected(reactions: Array) -> void:
	var rng = DeterministicRngScript.new(4042026)
	var samples: int = 4096
	var active_pairs: int = 0
	var total_pairs: int = samples * reactions.size()
	var sum_active_reactions: int = 0
	var universal_proteins: int = 0
	for _i in range(samples):
		var signature: int = int(rng.randi_range(0, 0xFFFF))
		var active_count: int = CatalyticLandscapeScript.active_reaction_count_for_signature(signature, reactions)
		sum_active_reactions += active_count
		active_pairs += active_count
		if active_count == reactions.size():
			universal_proteins += 1
	var active_fraction: float = float(active_pairs) / float(total_pairs)
	var mean_active_reactions: float = float(sum_active_reactions) / float(samples)
	print("M4 landscape active_fraction=%.6f mean_active_reactions=%.6f universal=%d" % [active_fraction, mean_active_reactions, universal_proteins])
	_assert_true(active_fraction >= 0.025 and active_fraction <= 0.055, "random protein/reaction activity density remains sparse but nonzero")
	_assert_true(mean_active_reactions >= 0.20 and mean_active_reactions <= 0.80, "random proteins have limited average catalytic promiscuity")
	_assert_true(universal_proteins == 0, "sampled random proteins are not universal catalysts")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
