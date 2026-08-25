extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const SnapshotCodecScript = preload("res://src/experiments/snapshot_codec.gd")
const CellMechanicsScript = preload("res://src/physics/cell_mechanics.gd")

const RECEPTOR_SIGNATURE: int = 0xE000
const MOTOR_SIGNATURE: int = 0xF000

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_ancestor_does_not_consume_motility_rng()
	_test_extracellular_turnover_is_live_in_engine()
	_test_constructed_receptor_motor_moves_and_pays_atp()
	_test_active_movement_respects_physical_wall()
	_test_m11_snapshot_roundtrip_preserves_dynamic_state_and_rng()
	if failures == 0:
		print("PASS: %d M11 authoritative integration tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M11 authoritative integration tests failed" % [failures, tests_run])
		quit(1)

func _config() -> Object:
	var config = SimConfigScript.new()
	config.world_width = 12
	config.world_height = 12
	config.max_cells = 16
	config.mutation_enabled = false
	return config

func _construct_molecular_motor_cell(sim, position: Vector2):
	var cell = sim.seed_ancestor(position)
	cell.genome.get_gene_by_locus(1).protein_signature = RECEPTOR_SIGNATURE
	cell.genome.get_gene_by_locus(2).protein_signature = MOTOR_SIGNATURE
	cell.initialize_molecular_state(sim.config)
	cell.set_pool("ATP", 20.0)
	cell.motor_heading = Vector2.RIGHT
	return cell

func _test_ancestor_does_not_consume_motility_rng() -> void:
	var sim = SimulationEngineScript.new(_config())
	sim.seed_ancestor()
	var before: int = int(sim.motility_rng.get_state())
	sim.step(5)
	_assert_true(int(sim.motility_rng.get_state()) == before, "non-motile ancestor consumes no motility RNG draws")
	_assert_close(float(sim.last_motility_summary.get("total_distance", -1.0)), 0.0, 1e-15, "non-motile ancestor has exactly zero active displacement")

func _test_extracellular_turnover_is_live_in_engine() -> void:
	var config = _config()
	config.extracellular_protein_decay_rate_per_min = log(2.0) / float(config.tick_dt_min)
	var sim = SimulationEngineScript.new(config)
	sim.world.release_protein(0xD136, Vector2(6, 6), 20.0, 0.0)
	var protein_before: float = sim.world.total_extracellular_protein()
	var aa_before: float = sim.world.get_field("amino_acids").total_amount()
	sim.step(1)
	var protein_after: float = sim.world.total_extracellular_protein()
	var aa_after: float = sim.world.get_field("amino_acids").total_amount()
	_assert_close(protein_before, 20.0, 1e-12, "engine turnover fixture starts with exact protein material")
	_assert_close(protein_after, 10.0, 1e-9, "engine tick applies configured extracellular protein half-life")
	_assert_close(aa_after - aa_before, 0.01, 1e-10, "engine tick recycles hydrolyzed protein material to extracellular AA")

func _test_constructed_receptor_motor_moves_and_pays_atp() -> void:
	var config = _config()
	config.motor_baseline_turn_rate_per_min = 0.0
	config.motor_speed_grid_per_min_per_activity = 1.0
	config.motor_atp_cost_per_grid_distance = 0.5
	var sim = SimulationEngineScript.new(config)
	var cell = _construct_molecular_motor_cell(sim, Vector2(4, 6))
	var x_before: float = cell.position.x
	var atp_before: float = cell.pool("ATP")
	sim.step(1)
	var moved: float = cell.position.x - x_before
	var motility_spent: float = float(cell.last_motility_summary.get("atp_spent", 0.0))
	_assert_true(moved > 0.0, "constructed realized motor produces active displacement in its physical heading")
	_assert_true(motility_spent > 0.0, "active displacement has explicit ATP charge in production simulation")
	_assert_true(cell.pool("ATP") < atp_before, "motor ATP cost is actually debited from authoritative cell pool")
	_assert_close(motility_spent, float(cell.last_motility_summary["actual_distance"]) * config.motor_atp_cost_per_grid_distance, 1e-10, "production motility ledger matches realized distance cost")

func _test_active_movement_respects_physical_wall() -> void:
	var config = _config()
	config.motor_baseline_turn_rate_per_min = 0.0
	config.motor_speed_grid_per_min_per_activity = 1000.0
	config.motor_atp_cost_per_grid_distance = 0.0
	var sim = SimulationEngineScript.new(config)
	var cell = _construct_molecular_motor_cell(sim, Vector2(10.0, 6.0))
	cell.motor_heading = Vector2.RIGHT
	sim.step(1)
	var radius: float = CellMechanicsScript.radius_for_cell(cell, config)
	_assert_true(cell.position.x <= float(config.world_width - 1) - radius + 1e-10, "active motor cannot push cell center through reflecting physical wall")

func _test_m11_snapshot_roundtrip_preserves_dynamic_state_and_rng() -> void:
	var config = _config()
	config.motor_baseline_turn_rate_per_min = 2.0
	var sim = SimulationEngineScript.new(config)
	var cell = _construct_molecular_motor_cell(sim, Vector2(5, 6))
	cell.signalling_state = {RECEPTOR_SIGNATURE: 17.25}
	cell.motor_heading = Vector2(0.6, 0.8)
	# Advance only the dedicated motor stream so restoration must preserve it.
	sim.motility_rng.randf()
	var snapshot: Dictionary = SnapshotCodecScript.capture(sim)
	var restored = SnapshotCodecScript.restore(SnapshotCodecScript.decode(SnapshotCodecScript.encode(snapshot)))
	_assert_true(int(snapshot["schema_version"]) == 3 and String(snapshot["model_identifier"]) == "microverse-m11", "M11 snapshot self-identifies new authoritative model schema")
	_assert_true(int(restored.motility_rng.get_state()) == int(sim.motility_rng.get_state()), "snapshot preserves independent motility RNG state exactly")
	_assert_true(restored.cells[0].signalling_state == cell.signalling_state, "snapshot preserves reversible molecular activation state")
	_assert_true(restored.cells[0].motor_heading == cell.motor_heading, "snapshot preserves physical motor orientation")
	_assert_close(restored.checksum(), sim.checksum(), 1e-12, "M11 snapshot roundtrip restores exact authoritative checksum")
	sim.step(10)
	restored.step(10)
	_assert_close(restored.checksum(), sim.checksum(), 1e-12, "M11 restored state replays exact future trajectory including stochastic motility")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
