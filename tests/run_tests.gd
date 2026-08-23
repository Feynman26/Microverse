extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const ChemicalFieldScript = preload("res://src/world/chemical_field.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_diffusion_conserves_mass_and_nonnegativity()
	_test_competing_identical_cells_are_order_fair()
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
	# Deliberately make demand exceed the finite local glucose pool so this test
	# exercises proportional allocation rather than merely equal requests.
	config.glucose_transport_vmax = 10.0
	var sim = SimulationEngineScript.new(config)
	var position := Vector2(4.0, 4.0)
	var first = sim.seed_ancestor(position)
	var second = sim.seed_ancestor(position)
	sim.step(1)
	_assert_close(first.internal_glucose, second.internal_glucose, 1e-12, "identical competitors receive equal scarce resource")
	_assert_true(sim.world.get_field("glucose").get_value(4, 4) <= 1e-12, "scarce local glucose is exhausted by proportional allocation")

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
	var saw_division := false
	for event in sim.event_log:
		if event["kind"] == "division":
			saw_division = true
			break
	_assert_true(saw_division, "division event is recorded")

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
