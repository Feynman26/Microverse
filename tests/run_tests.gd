extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const ChemicalFieldScript = preload("res://src/world/chemical_field.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const MutationEngineScript = preload("res://src/genetics/mutation_engine.gd")
const ReactionCatalogScript = preload("res://src/chemistry/reaction_catalog.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const MainUiScript = preload("res://src/ui/main.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_assert_true(MainUiScript != null, "main UI script parses")
	_test_diffusion_conserves_mass_and_nonnegativity()
	_test_competing_identical_cells_are_transport_order_fair()
	_test_division_conserves_metabolites_and_expression_amounts()
	_test_ancestor_genome_is_stable_and_copyable()
	_test_mutation_disabled_inheritance_is_deep_copy()
	_test_forced_mutations_change_only_child_genotype()
	_test_same_seed_reproduces_exact_mutation_sequence()
	_test_neutral_marker_remains_physically_neutral()
	_test_mutation_frequency_matches_configured_probability()
	_test_simulation_records_mutation_ancestry()
	_test_cell_grows_and_divides_with_resources()
	_test_cell_dies_without_energy_source()
	_test_same_seed_reproduces_same_history_including_expression_noise()

	if failures == 0:
		print("PASS: %d Microverse tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d Microverse tests failed" % [failures, tests_run])
		quit(1)

func _test_diffusion_conserves_mass_and_nonnegativity() -> void:
	var field = ChemicalFieldScript.new(9, 9, 1.0, 1.0, 0.0)
	field.set_value(4, 4, 10.0)
	var before: float = field.total_amount()
	for _i in range(100):
		field.step_diffusion(0.1)
	_assert_close(field.total_amount(), before, 1e-9, "diffusion conserves closed-chamber mass")
	_assert_true(field.minimum_value() >= 0.0, "diffusion remains nonnegative")

func _test_competing_identical_cells_are_transport_order_fair() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	config.initial_glucose = 0.001
	config.glucose_diffusion = 0.0
	config.oxygen_diffusion = 0.0
	config.nitrogen_diffusion = 0.0
	config.phosphorus_diffusion = 0.0
	config.glucose_transport_vmax = 10.0
	config.metabolic_rate_scale = 1e-12
	var sim = SimulationEngineScript.new(config)
	var position: Vector2 = Vector2(4.0, 4.0)
	var first = sim.seed_ancestor(position)
	var second = sim.seed_ancestor(position)
	sim.step(1)
	_assert_close(first.pool("G"), second.pool("G"), 1e-10, "identical competitors receive equal scarce carbon before material metabolic consumption")
	_assert_true(sim.world.get_field("glucose").get_value(4, 4) <= 1e-12, "scarce local glucose is exhausted by proportional allocation")

func _test_division_conserves_metabolites_and_expression_amounts() -> void:
	var config = SimConfigScript.new()
	var rng = DeterministicRngScript.new(12345)
	var world = WorldStateScript.new(8, 8, 1.0)
	var parent = CellStateScript.new(1, -1, 0, 0, Vector2(4.0, 4.0), config.division_volume)
	parent.genome = GenomeScript.create_ancestor()
	parent.initialize_molecular_state(config)
	parent.set_pool("BIO", config.division_volume * config.biomass_units_per_volume)
	parent.volume = config.division_volume
	parent.set_pool("G", 1.2)
	parent.set_pool("O2", 2.4)
	parent.set_pool("ATP", 5.0)
	parent.set_pool("ADP", 5.0)
	parent.set_pool("NAD", 3.0)
	parent.set_pool("NADH", 1.0)
	parent.damage = 0.2
	parent.energy_debt = 0.1
	var metabolite_before: Dictionary = parent.metabolites.duplicate(true)
	var mrna_before: float = parent.total_mrna()
	var protein_before: float = parent.total_protein()
	var daughters: Array = parent.create_daughters(2, 3, 7, rng, world, config)
	var a = daughters[0]
	var b = daughters[1]
	_assert_close(a.volume + b.volume, config.division_volume, 1e-12, "division conserves structural cell volume")
	for metabolite_id in metabolite_before.keys():
		var expected: float = float(metabolite_before[metabolite_id])
		if String(metabolite_id) == "ATP":
			expected -= config.division_atp_cost
		elif String(metabolite_id) == "ADP":
			expected += config.division_atp_cost
		_assert_close(a.pool(String(metabolite_id)) + b.pool(String(metabolite_id)), expected, 1e-10, "division conserves %s after explicit ATP->ADP cost" % metabolite_id)
	_assert_close(a.total_mrna() + b.total_mrna(), mrna_before, 1e-9, "stochastic partition conserves total inherited mRNA")
	_assert_close(a.total_protein() + b.total_protein(), protein_before, 1e-9, "stochastic partition conserves total inherited protein")
	_assert_true(a.parent_id == 1 and b.parent_id == 1, "both daughters retain immutable parent identity")
	_assert_true(a.generation == 1 and b.generation == 1, "both daughters advance generation")

func _test_ancestor_genome_is_stable_and_copyable() -> void:
	var ancestor = GenomeScript.create_ancestor()
	var copied = ancestor.deep_copy()
	_assert_true(ancestor.gene_count() == 12, "M3 ancestor contains twelve discrete loci")
	_assert_true(ancestor.exact_equals(copied), "deep-copied genome is exactly equal")
	_assert_true(ancestor.fingerprint() == copied.fingerprint(), "equal genomes share deterministic fingerprint")
	_assert_true(ancestor != copied and ancestor.genes[0] != copied.genes[0], "genome copy has no shared mutable gene objects")

func _test_mutation_disabled_inheritance_is_deep_copy() -> void:
	var config = SimConfigScript.new()
	config.mutation_enabled = false
	var rng = DeterministicRngScript.new(7721)
	var world = WorldStateScript.new(8, 8, 1.0)
	var parent = CellStateScript.new(10, -1, 0, 0, Vector2(4.0, 4.0), config.division_volume)
	parent.genome = GenomeScript.create_ancestor()
	parent.initialize_molecular_state(config)
	parent.set_pool("BIO", config.division_volume * config.biomass_units_per_volume)
	parent.volume = config.division_volume
	parent.set_pool("ATP", 5.0)
	var parent_key: String = parent.genome.canonical_key()
	var daughters: Array = parent.create_daughters(11, 12, 1, rng, world, config)
	var a = daughters[0]
	var b = daughters[1]
	_assert_true(a.genome.canonical_key() == parent_key and b.genome.canonical_key() == parent_key, "mutation-free daughters inherit exact parental genotype")
	_assert_true(a.genome != b.genome and a.genome.genes[0] != b.genome.genes[0], "sister genomes are independent mutable copies")
	_assert_true(a.expression_state != b.expression_state, "sister expression states are independent containers")

func _forced_mutation_config():
	var config = SimConfigScript.new()
	# This is the historical M3 forced-Bernoulli harness. Disable M10 physical
	# replication here explicitly so it continues testing the old mutation API,
	# not the production replication-derived fidelity path.
	config.evolvable_replication_enabled = false
	config.promoter_mutation_rate_per_gene = 1.0
	config.signature_mutation_rate_per_gene = 1.0
	config.neutral_marker_mutation_rate_per_gene = 1.0
	return config

func _test_forced_mutations_change_only_child_genotype() -> void:
	var config = _forced_mutation_config()
	var rng = DeterministicRngScript.new(44321)
	var mutator = MutationEngineScript.new()
	var parent = GenomeScript.create_ancestor()
	var parent_key_before: String = parent.canonical_key()
	var result: Dictionary = mutator.mutate_copy(parent, rng, config)
	var child = result["genome"]
	var events: Array = result["events"]
	_assert_true(parent.canonical_key() == parent_key_before, "mutation engine never mutates parental genome")
	_assert_true(not parent.exact_equals(child), "forced molecular mutations alter child genotype")
	_assert_true(events.size() == parent.gene_count() * 3, "forced mutation emits promoter signature and neutral event per locus")

func _test_same_seed_reproduces_exact_mutation_sequence() -> void:
	var config = _forced_mutation_config()
	var mutator = MutationEngineScript.new()
	var first: Dictionary = mutator.mutate_copy(GenomeScript.create_ancestor(), DeterministicRngScript.new(91731), config)
	var second: Dictionary = mutator.mutate_copy(GenomeScript.create_ancestor(), DeterministicRngScript.new(91731), config)
	_assert_true(first["genome"].canonical_key() == second["genome"].canonical_key(), "same seed reproduces exact mutated genotype")
	_assert_true(first["events"] == second["events"], "same seed reproduces exact mutation event sequence")

func _test_neutral_marker_remains_physically_neutral() -> void:
	var config = SimConfigScript.new()
	var world = WorldStateScript.new(8, 8, 1.0)
	world.register_field("glucose", 0.0, 2.0)
	world.register_field("oxygen", 0.0, 2.0)
	world.register_field("nitrogen", 0.0, 2.0)
	world.register_field("phosphorus", 0.0, 2.0)
	var base_genome = GenomeScript.create_ancestor()
	var neutral_genome = base_genome.deep_copy()
	neutral_genome.genes[0].neutral_marker += 1
	var reactions: Array = ReactionCatalogScript.create_m4_candidate()
	var a = CellStateScript.new(1, -1, 0, 0, Vector2(4.0, 4.0), 1.0)
	var b = CellStateScript.new(2, -1, 0, 0, Vector2(4.0, 4.0), 1.0)
	a.genome = base_genome
	b.genome = neutral_genome
	a.initialize_molecular_state(config)
	b.initialize_molecular_state(config)
	for metabolite_id in ["G", "O2", "NH4", "P", "AA", "NUC"]:
		a.set_pool(metabolite_id, 0.5)
		b.set_pool(metabolite_id, 0.5)
	var request_a: Dictionary = a.transport_requests(config.tick_dt_min, world, config)
	var request_b: Dictionary = b.transport_requests(config.tick_dt_min, world, config)
	_assert_true(request_a == request_b, "neutral marker does not alter membrane requests")
	a.step_intracellular(config.tick_dt_min, config, reactions, DeterministicRngScript.new(811))
	b.step_intracellular(config.tick_dt_min, config, reactions, DeterministicRngScript.new(811))
	_assert_close(a.pool("ATP"), b.pool("ATP"), 1e-12, "neutral marker does not alter ATP physiology under identical molecular noise")
	_assert_close(a.volume, b.volume, 1e-12, "neutral marker does not alter biomass physiology")
	_assert_close(a.damage, b.damage, 1e-12, "neutral marker does not alter damage physiology")
	_assert_close(ExpressionSystemScript.checksum(a.expression_state), ExpressionSystemScript.checksum(b.expression_state), 1e-12, "neutral marker does not alter expression trajectory")

func _test_mutation_frequency_matches_configured_probability() -> void:
	var config = SimConfigScript.new()
	config.promoter_mutation_rate_per_gene = 0.0
	config.signature_mutation_rate_per_gene = 0.0
	config.neutral_marker_mutation_rate_per_gene = 0.05
	var rng = DeterministicRngScript.new(552211)
	var mutator = MutationEngineScript.new()
	var trials: int = 300
	var loci: int = GenomeScript.create_ancestor().gene_count()
	var observed: int = 0
	for _i in range(trials):
		var result: Dictionary = mutator.mutate_copy(GenomeScript.create_ancestor(), rng, config)
		observed += int(result["events"].size())
	var expected: float = float(trials * loci) * config.neutral_marker_mutation_rate_per_gene
	var sigma: float = sqrt(float(trials * loci) * 0.05 * 0.95)
	_assert_true(absf(float(observed) - expected) <= 5.0 * sigma, "realized neutral mutation count is within five sigma of configured Bernoulli rate")

func _test_simulation_records_mutation_ancestry() -> void:
	var config = _forced_mutation_config()
	config.world_width = 16
	config.world_height = 16
	config.max_cells = 4
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor()
	sim.step(1600)
	_assert_true(sim.mutation_event_count() > 0, "simulation records molecular mutation events")
	var checked: bool = false
	for event in sim.event_log:
		if event["kind"] == "mutation":
			checked = true
			_assert_true(event.has("mutation_id") and event.has("cell_id") and event.has("parent_id"), "mutation event contains stable ancestry identifiers")
			_assert_true(event.has("parent_genotype_fingerprint") and event.has("resulting_genotype_fingerprint"), "mutation event links parental and resulting genotype")
			break
	_assert_true(checked, "at least one mutation event is inspectable")

func _test_cell_grows_and_divides_with_resources() -> void:
	var config = SimConfigScript.new()
	config.world_width = 16
	config.world_height = 16
	config.max_cells = 8
	config.mutation_enabled = false
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor()
	sim.step(1600)
	_assert_true(sim.population_size() > 1, "resource-fed ancestor produces descendants through explicit expression plus chemistry")
	_assert_true(sim.maximum_generation() >= 1, "molecular expression and metabolic biomass accumulation advance generation")
	var division_events: int = 0
	for event in sim.event_log:
		if event["kind"] == "division":
			division_events += 1
	_assert_true(division_events > 0, "expression-supported division event is recorded")

func _test_cell_dies_without_energy_source() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	config.initial_glucose = 0.0
	config.initial_oxygen = 0.0
	config.initial_nitrogen = 0.0
	config.initial_phosphorus = 0.0
	config.glucose_diffusion = 0.0
	config.oxygen_diffusion = 0.0
	config.nitrogen_diffusion = 0.0
	config.phosphorus_diffusion = 0.0
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor()
	sim.step(1000)
	_assert_true(sim.population_size() == 0, "cell eventually dies when ATP maintenance/expression cannot be replenished")

func _test_same_seed_reproduces_same_history_including_expression_noise() -> void:
	var config_a = SimConfigScript.new()
	config_a.world_width = 16
	config_a.world_height = 16
	config_a.max_cells = 12
	config_a.seed = 99173
	var config_b = SimConfigScript.new()
	config_b.world_width = 16
	config_b.world_height = 16
	config_b.max_cells = 12
	config_b.seed = 99173
	var first = SimulationEngineScript.new(config_a)
	var second = SimulationEngineScript.new(config_b)
	first.seed_ancestor()
	second.seed_ancestor()
	first.step(900)
	second.step(900)
	_assert_true(first.population_size() == second.population_size(), "same seed reproduces population size with expression noise")
	_assert_true(first.event_log == second.event_log, "same seed reproduces exact semantic event history")
	_assert_close(first.checksum(), second.checksum(), 1e-9, "same seed reproduces chemistry genetics and expression state")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])