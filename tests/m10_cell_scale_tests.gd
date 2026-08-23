extends SceneTree

const SimConfigScript = preload("res://src/core/sim_config.gd")
const DeterministicRngScript = preload("res://src/core/deterministic_rng.gd")
const WorldStateScript = preload("res://src/world/world_state.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")
const ExpressionSystemScript = preload("res://src/expression/expression_system.gd")
const MetabolicSolverScript = preload("res://src/chemistry/metabolic_solver.gd")

var failures: int = 0
var tests_run: int = 0

func _init() -> void:
	call_deferred("_run")

func _run() -> void:
	_test_capacity_and_equilibrium_scale_with_biomass()
	_test_sensing_is_concentration_based()
	_test_division_preserves_machinery_concentration()
	if failures == 0:
		print("PASS: %d M10 cell-scale tests" % tests_run)
		quit(0)
	else:
		push_error("FAIL: %d of %d M10 cell-scale tests failed" % [failures, tests_run])
		quit(1)

func _test_capacity_and_equilibrium_scale_with_biomass() -> void:
	var config = SimConfigScript.new()
	var unit = _new_cell(1, 1.0, config)
	var double = _new_cell(2, 2.0, config)

	_assert_close(double.proteome_capacity(config), 2.0 * unit.proteome_capacity(config), 1e-10, "doubling biomass doubles physical proteome capacity")
	_assert_close(double.total_protein(), 2.0 * unit.total_protein(), 1e-10, "equilibrium protein inventory scales with biomass instead of remaining fixed per cell")
	_assert_close(double.total_mrna(), 2.0 * unit.total_mrna(), 1e-10, "equilibrium mRNA inventory scales with biomass")
	_assert_close(double.total_protein() / double.volume, unit.total_protein() / unit.volume, 1e-10, "equal genotype at equal state keeps protein concentration invariant across cell size")

	var genome = GenomeScript.create_ancestor()
	var state_unit: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 1.0)
	var state_double: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 2.0)
	var pools_unit: Dictionary = _abundant_expression_pools(1.0, config)
	var pools_double: Dictionary = _abundant_expression_pools(2.0, config)
	var summary_unit: Dictionary = ExpressionSystemScript.step(state_unit, genome, pools_unit, config.tick_dt_min, DeterministicRngScript.new(701), config, 1.0)
	var summary_double: Dictionary = ExpressionSystemScript.step(state_double, genome, pools_double, config.tick_dt_min, DeterministicRngScript.new(701), config, 2.0)
	_assert_close(float(summary_double["translation_capacity"]), 2.0 * float(summary_unit["translation_capacity"]), 1e-12, "shared translation throughput scales linearly with physical cell size")

func _test_sensing_is_concentration_based() -> void:
	var config = SimConfigScript.new()
	var genome = GenomeScript.create_ancestor()
	var state_unit: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 1.0)
	var state_double: Dictionary = ExpressionSystemScript.create_equilibrium_state(genome, config, 2.0)
	var pools_unit: Dictionary = _abundant_expression_pools(1.0, config)
	var pools_double: Dictionary = _abundant_expression_pools(2.0, config)
	# Give both cells the same ligand concentrations by doubling every molecular
	# amount in the larger cell. Regulation must therefore see the same chemistry.
	var summary_unit: Dictionary = ExpressionSystemScript.step(state_unit, genome, pools_unit, config.tick_dt_min, DeterministicRngScript.new(991), config, 1.0)
	var summary_double: Dictionary = ExpressionSystemScript.step(state_double, genome, pools_double, config.tick_dt_min, DeterministicRngScript.new(991), config, 2.0)
	for locus_variant in summary_unit["regulation"].keys():
		var locus_id: int = int(locus_variant)
		_assert_close(
			float(summary_double["regulation"][locus_id]),
			float(summary_unit["regulation"][locus_id]),
			1e-12,
			"regulatory occupancy for locus %d depends on concentration rather than absolute cell size" % locus_id
		)

func _test_division_preserves_machinery_concentration() -> void:
	var config = SimConfigScript.new()
	config.partition_jitter = 0.0
	config.expression_partition_noise_scale = 0.0
	var mother = _new_cell(1, config.division_volume, config)
	var world = WorldStateScript.new(config.world_width, config.world_height, config.grid_cell_size_um)
	var rng = DeterministicRngScript.new(424242)
	var protein_before: float = mother.total_protein()
	var mrna_before: float = mother.total_mrna()
	var protein_concentration_before: float = protein_before / mother.volume
	var mrna_concentration_before: float = mrna_before / mother.volume

	_assert_true(mother.ready_to_divide(config), "controlled biomass-scaled mother is ordinarily ready to divide")
	var daughters: Array = mother.create_daughters(2, 3, 10, rng, world, config)
	_assert_close(float(daughters[0].volume), 1.0, 1e-10, "symmetric division returns first daughter to reference cell volume")
	_assert_close(float(daughters[1].volume), 1.0, 1e-10, "symmetric division returns second daughter to reference cell volume")
	_assert_close(float(daughters[0].total_protein() + daughters[1].total_protein()), protein_before, 1e-10, "division conserves explicit protein inventory")
	_assert_close(float(daughters[0].total_mrna() + daughters[1].total_mrna()), mrna_before, 1e-10, "division conserves explicit mRNA inventory")
	for daughter in daughters:
		_assert_close(daughter.total_protein() / daughter.volume, protein_concentration_before, 1e-10, "daughter inherits maternal protein concentration without artificial post-division dilution")
		_assert_close(daughter.total_mrna() / daughter.volume, mrna_concentration_before, 1e-10, "daughter inherits maternal mRNA concentration without artificial post-division dilution")
		_assert_true(daughter.total_protein() <= daughter.proteome_capacity(config) + 1e-8, "daughter inherited proteome remains inside its volume-scaled physical capacity")

func _new_cell(cell_id: int, volume: float, config):
	var cell = CellStateScript.new(cell_id, -1, 0, 0, Vector2(4.0, 4.0), volume)
	cell.genome = GenomeScript.create_ancestor()
	cell.initialize_molecular_state(config)
	return cell

func _abundant_expression_pools(volume: float, config) -> Dictionary:
	var pools: Dictionary = MetabolicSolverScript.create_initial_pools(volume, config)
	# Scale all relevant molecular inventories with volume so the two controlled
	# cells begin at identical concentrations. The large values prevent resource
	# scarcity from obscuring the machinery-capacity test.
	pools["ATP"] = 1000.0 * volume
	pools["ADP"] = 1000.0 * volume
	pools["AA"] = 1000.0 * volume
	pools["NUC"] = 1000.0 * volume
	pools["O2"] = 2.0 * volume
	pools["X"] = 1.0 * volume
	return pools

func _assert_true(condition: bool, message: String) -> void:
	tests_run += 1
	if condition:
		print("PASS: %s" % message)
	else:
		failures += 1
		push_error("FAIL: %s" % message)

func _assert_close(actual: float, expected: float, tolerance: float, message: String) -> void:
	_assert_true(absf(actual - expected) <= tolerance, "%s (actual=%s expected=%s)" % [message, actual, expected])