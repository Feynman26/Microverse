extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const DNAReplicationScript = preload("res://src/genetics/dna_replication.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const MetaboliteCatalogScript = preload("res://src/chemistry/metabolite_catalog.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_resident_and_partial_dna_return_on_lysis()
	_test_structural_expansion_requires_nucleotide_material()
	_test_structural_reduction_releases_nucleotide_material()
	if failures == 0:
		print("PASS: %d M10 DNA-material tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 DNA-material tests failed" % [failures, tests_run])
		quit(1)

func _test_resident_and_partial_dna_return_on_lysis() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	var sim = SimulationEngineScript.new(config)
	var cell = sim.seed_ancestor(Vector2(4.0, 4.0))
	var nuc_field: String = MetaboliteCatalogScript.extracellular_field("NUC")
	var resident: float = DNAReplicationScript.genome_nuc_material(cell.genome, config)
	var baseline_release: float = float(cell.releasable_pools(config)[nuc_field])

	cell.replication_nuc_spent = 0.0375
	var partial_release: float = float(cell.releasable_pools(config)[nuc_field])
	_assert_close(partial_release - baseline_release, 0.0375, 1e-12, "partially synthesized DNA returns its exact consumed nucleotide material on lysis")
	_assert_true(resident > 0.0, "resident chromosome carries explicit derivable nucleotide material")
	_assert_close(
		DNAReplicationScript.total_cell_dna_nuc_material(cell, config),
		resident + 0.0375,
		1e-12,
		"cell DNA inventory equals resident chromosome plus nascent copy material"
	)

	var before_world: float = sim.world.get_field(nuc_field).total_amount()
	var expected_release: float = float(cell.releasable_pools(config)[nuc_field])
	cell.alive = false
	cell.death_reason = "dna_material_fixture"
	sim._process_deaths()
	_assert_close(
		sim.world.get_field(nuc_field).total_amount() - before_world,
		expected_release,
		1e-9,
		"cell death returns resident plus partial DNA nucleotide material to extracellular NUC"
	)

func _test_structural_expansion_requires_nucleotide_material() -> void:
	var config = SimConfigScript.new()
	var profile: Dictionary = {
		"point_error_rate_per_gene": 0.0,
		"structural_error_rate_per_genome": 1.0
	}
	var mutator = MutationEngineScript.new()
	var expansion_seed: int = -1
	var blocked: Dictionary = {}
	for seed in range(1, 512):
		var candidate: Dictionary = mutator.mutate_replicated_copy(
			GenomeScript.create_ancestor(), DeterministicRngScript.new(seed), config, profile, 0.0
		)
		var events: Array = candidate["events"]
		if events.is_empty():
			continue
		var event: Dictionary = events[-1]
		if String(event.get("reason", "")) == "insufficient_nucleotide_material":
			expansion_seed = seed
			blocked = candidate
			break
	_assert_true(expansion_seed > 0, "deterministic structural-error search finds a DNA-expanding event")
	if expansion_seed <= 0:
		return
	_assert_close(float(blocked["dna_nuc_material_delta"]), 0.0, 1e-12, "zero free NUC blocks structural DNA expansion without manufacturing material")
	_assert_true(int(blocked["genome"].gene_count()) == 12, "blocked structural expansion leaves ancestral coding-locus count intact")

	var funded: Dictionary = mutator.mutate_replicated_copy(
		GenomeScript.create_ancestor(), DeterministicRngScript.new(expansion_seed), config, profile, 100.0
	)
	_assert_true(float(funded["dna_nuc_material_delta"]) > 0.0, "same sampled structural expansion proceeds when nucleotide material is available")
	_assert_true(
		float(funded["genome"].replication_unit_count()) > float(GenomeScript.create_ancestor().replication_unit_count()),
		"funded structural expansion increases physical DNA replication units"
	)

func _test_structural_reduction_releases_nucleotide_material() -> void:
	var config = SimConfigScript.new()
	var profile: Dictionary = {
		"point_error_rate_per_gene": 0.0,
		"structural_error_rate_per_genome": 1.0
	}
	var mutator = MutationEngineScript.new()
	var reduction: Dictionary = {}
	for seed in range(1, 512):
		var candidate: Dictionary = mutator.mutate_replicated_copy(
			GenomeScript.create_ancestor(), DeterministicRngScript.new(seed), config, profile, 100.0
		)
		if float(candidate.get("dna_nuc_material_delta", 0.0)) < 0.0:
			reduction = candidate
			break
	_assert_true(not reduction.is_empty(), "deterministic structural-error search finds a DNA-reducing event")
	if reduction.is_empty():
		return
	var delta: float = float(reduction["dna_nuc_material_delta"])
	_assert_true(delta < 0.0, "DNA-reducing mutation reports negative material delta for nucleotide recycling")
	_assert_close(
		-delta,
		(float(GenomeScript.create_ancestor().replication_unit_count()) - float(reduction["genome"].replication_unit_count())) * float(config.genome_replication_nuc_cost_per_gene),
		1e-12,
		"released nucleotide material equals exact lost DNA replication units"
	)

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])