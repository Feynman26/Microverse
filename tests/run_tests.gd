extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const ChemicalFieldScript = preload("res://src/world/chemical_field.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const MainUiScript = preload("res://src/ui/main.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	# Referencing the preloaded UI script makes parse failures part of this suite.
	_assert_true(MainUiScript != null, "main UI script parses")
	_test_diffusion_conserves_mass_and_nonnegativity()
	_test_competing_identical_cells_are_order_fair()
	_test_division_conserves_partitioned_pools()
	_test_cell_grows_and_divides_with_resources()
	_test_cell_dies_without_energy_source()
	_test_same_seed_reproduces_same_history()

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

func _test_competing_identical_cells_are_order_fair() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	config.initial_glucose = 0.001
	config.initial_oxygen = 5.0
	config.glucose_diffusion = 0.0
	config.oxygen_diffusion = 0.0
	config.glucose_transport_vmax = 10.0
	var sim = SimulationEngineScript.new(config)
	var position: Vector2 = Vector2(4.0, 4.0)
	var first = sim.seed_ancestor(position)
	var second = sim.seed_ancestor(position)
	sim.step(1)
	_assert_close(first.internal_glucose, second.internal_glucose, 1e-12, "identical competitors receive equal scarce resource")
	_assert_true(sim.world.get_field("glucose").get_value(4, 4) <= 1e-12, "scarce local glucose is exhausted by proportional allocation")

func _test_division_conserves_partitioned_pools() -> void:
	var config = SimConfigScript.new()
	var rng = DeterministicRngScript.new(12345)
	var world = WorldStateScript.new(8, 8, 1.0)
	var parent = CellStateScript.new(1, -1, 0, 0, Vector2(4.0, 4.0), config.division_volume)
	parent.internal_glucose = 1.2
	parent.internal_oxygen = 2.4
	parent.atp = 5.0
	parent.precursor = 0.8
	parent.ros = 0.3
	parent.damage = 0.2
	parent.energy_debt = 0.1

	var expected_atp: float = float(parent.atp) - float(config.division_atp_cost)
	var daughters: Array = parent.create_daughters(2, 3, 7, rng, world, config)
	var a = daughters[0]
	var b = daughters[1]
	_assert_close(a.volume + b.volume, config.division_volume, 1e-12, "division conserves cell volume")
	_assert_close(a.internal_glucose + b.internal_glucose, 1.2, 1e-12, "division conserves intracellular glucose")
	_assert_close(a.internal_oxygen + b.internal_oxygen, 2.4, 1e-12, "division conserves intracellular oxygen")
	_assert_close(a.atp + b.atp, expected_atp, 1e-12, "division conserves ATP after explicit division cost")
	_assert_close(a.precursor + b.precursor, 0.8, 1e-12, "division conserves precursor")
	_assert_true(a.parent_id == 1 and b.parent_id == 1, "both daughters retain immutable parent identity")
	_assert_true(a.generation == 1 and b.generation == 1, "both daughters advance generation")

func _test_cell_grows_and_divides_with_resources() -> void:
	var config = SimConfigScript.new()
	config.world_width = 16
	config.world_height = 16
	config.max_cells = 8
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor()
	sim.step(300)
	_assert_true(sim.population_size() > 1, "resource-fed ancestor produces descendants")
	_assert_true(sim.maximum_generation() >= 1, "division advances generation")
	var division_events: int = 0
	var birth_events: int = 0
	for event in sim.event_log:
		if event["kind"] == "division":
			division_events += 1
		elif event["kind"] == "birth":
			birth_events += 1
	_assert_true(division_events > 0, "division event is recorded")
	_assert_true(birth_events >= 3, "daughter births are recorded explicitly")

func _test_cell_dies_without_energy_source() -> void:
	var config = SimConfigScript.new()
	config.world_width = 8
	config.world_height = 8
	config.initial_glucose = 0.0
	config.initial_oxygen = 0.0
	config.glucose_diffusion = 0.0
	config.oxygen_diffusion = 0.0
	var sim = SimulationEngineScript.new(config)
	sim.seed_ancestor()
	sim.step(1200)
	_assert_true(sim.population_size() == 0, "cell eventually dies when maintenance cannot be paid")

func _test_same_seed_reproduces_same_history() -> void:
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
	first.step(500)
	second.step(500)

	_assert_true(first.population_size() == second.population_size(), "same seed reproduces population size")
	_assert_true(first.event_log.size() == second.event_log.size(), "same seed reproduces event count")
	_assert_close(first.checksum(), second.checksum(), 1e-9, "same seed reproduces complete state checksum")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
