extends SceneTree

const WorldStateScript = preload("res://src/world/world_state.gd")
const TurnoverScript = preload("res://src/chemistry/extracellular_protein_turnover.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_half_life_recycles_local_material()
	_test_sequence_identity_does_not_change_decay_law()
	_test_exhausted_signature_field_is_pruned_without_material_loss()
	if failures == 0:
		print("PASS: %d M11 extracellular-protein-turnover tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M11 extracellular-protein-turnover tests failed" % [failures, tests_run])
		quit(1)

func _world():
	var world = WorldStateScript.new(5, 5, 1.0)
	world.register_field("amino_acids", 0.0, 0.0)
	return world

func _test_half_life_recycles_local_material() -> void:
	var world = _world()
	var site := Vector2(2, 2)
	world.release_protein(0xD136, site, 100.0, 0.0)
	var before_protein: float = world.total_extracellular_protein()
	var summary: Dictionary = TurnoverScript.step(world, 1.0, log(2.0), 0.001)
	var remaining: float = world.total_extracellular_protein()
	var aa: float = world.get_field("amino_acids").total_amount()
	_assert_close(before_protein, 100.0, 1e-12, "turnover fixture starts with exact extracellular protein material")
	_assert_close(remaining, 50.0, 1e-10, "one configured half-life leaves half the protein material")
	_assert_close(float(summary["degraded_protein"]), 50.0, 1e-10, "turnover ledger reports physically degraded protein")
	_assert_close(aa, 0.05, 1e-12, "degraded protein recycles exact modeled amino-acid material")
	_assert_close(world.get_field("amino_acids").get_value(2, 2), 0.05, 1e-12, "protein hydrolysis returns AA at the same spatial site")
	_assert_close(world.get_field("amino_acids").get_value(0, 0), 0.0, 1e-12, "turnover does not teleport recycled material across the chamber")

func _test_sequence_identity_does_not_change_decay_law() -> void:
	var world = _world()
	world.release_protein(0x1111, Vector2(1, 2), 40.0, 0.0)
	world.release_protein(0xEEEE, Vector2(3, 2), 40.0, 0.0)
	TurnoverScript.step(world, 2.0, 0.2, 0.001)
	var first: float = world.get_protein_field(0x1111).total_amount()
	var second: float = world.get_protein_field(0xEEEE).total_amount()
	_assert_close(first, second, 1e-12, "generic extracellular turnover assigns no sequence-specific survival bonus")
	_assert_true(first < 40.0 and first > 0.0, "finite environmental turnover removes protein gradually")

func _test_exhausted_signature_field_is_pruned_without_material_loss() -> void:
	var world = _world()
	world.release_protein(0xABCD, Vector2(2, 2), 12.5, 0.0)
	var summary: Dictionary = TurnoverScript.step(world, 1.0, 1000000.0, 0.001)
	_assert_true(not world.has_protein_field(0xABCD), "fully hydrolyzed extracellular sequence field is removed from dynamic world state")
	_assert_true(summary["removed_signatures"].has(0xABCD), "turnover ledger records pruned sequence field")
	_assert_close(float(summary["degraded_protein"]), 12.5, 1e-10, "field pruning accounts for all original protein material")
	_assert_close(world.get_field("amino_acids").total_amount(), 0.0125, 1e-12, "field pruning recycles all protein material to amino acids")

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])
