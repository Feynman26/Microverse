extends RefCounted
class_name CausalLab

const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const SnapshotCodecScript = preload("res://src/experiments/snapshot_codec.gd")
const EnvironmentScheduleScript = preload("res://src/experiments/environment_schedule.gd")
const GeneticAssayScript = preload("res://src/experiments/genetic_assay.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const CellMechanicsScript = preload("res://src/physics/cell_mechanics.gd")

static func replay(snapshot: Dictionary, target_tick: int, environment: Dictionary = {}):
	var sim = SnapshotCodecScript.restore(snapshot)
	assert(target_tick >= int(sim.tick_index))
	while int(sim.tick_index) < target_tick:
		if not environment.is_empty():
			EnvironmentScheduleScript.apply(sim, int(sim.tick_index), environment)
		sim.step(1)
	return sim

static func verify_replay(snapshot: Dictionary, target_tick: int, expected_checksum: float, environment: Dictionary = {}) -> Dictionary:
	var sim = replay(snapshot, target_tick, environment)
	return {
		"target_tick": target_tick,
		"expected_checksum": expected_checksum,
		"actual_checksum": float(sim.checksum()),
		"matches": absf(float(sim.checksum()) - expected_checksum) <= 1e-9
	}

# Select the latest persisted state that does not lie after the requested event.
# This is a pure lookup over supplied snapshots; it never mutates or advances them.
static func nearest_prior_snapshot(snapshots: Array, target_tick: int) -> Dictionary:
	assert(target_tick >= 0)
	var found: bool = false
	var best_tick: int = -1
	var best_snapshot: Dictionary = {}
	for snapshot_variant in snapshots:
		var snapshot: Dictionary = snapshot_variant
		assert(snapshot.has("tick_index"))
		var tick: int = int(snapshot["tick_index"])
		if tick <= target_tick and (not found or tick > best_tick):
			found = true
			best_tick = tick
			best_snapshot = snapshot
	assert(found, "No snapshot exists at or before requested target tick")
	return best_snapshot.duplicate(true)

static func semantic_timeline(event_log: Array, from_tick: int = 0, to_tick: int = 2147483647) -> Array:
	var result: Array = []
	for event_variant in event_log:
		var event: Dictionary = event_variant
		var tick: int = int(event.get("tick", -1))
		if tick < from_tick or tick > to_tick:
			continue
		if String(event.get("kind", "")) in ["birth", "division", "death", "mutation", "experimental_inoculation", "genetic_intervention", "environment_intervention", "detected_transition"]:
			result.append(event.duplicate(true))
	return result

static func run_competition(
	base_config,
	strain_entries: Array,
	initial_counts: Array,
	environment: Dictionary,
	horizon_ticks: int,
	seed: int,
	positions: Array = []
) -> Dictionary:
	assert(not strain_entries.is_empty())
	assert(strain_entries.size() == initial_counts.size())
	assert(horizon_ticks > 0)
	var config = SnapshotCodecScript.restore_config(SnapshotCodecScript.capture_config(base_config))
	config.seed = seed
	config.mutation_enabled = false
	var requested_founders: int = 0
	for count_variant in initial_counts:
		assert(int(count_variant) > 0)
		requested_founders += int(count_variant)
	assert(requested_founders <= int(config.max_cells))
	var sim = SimulationEngineScript.new(config)
	var founder_to_strain: Dictionary = {}
	var founders_by_strain: Dictionary = {}
	var descendants_by_founder: Dictionary = {}
	var position_index: int = 0
	for strain_index in range(strain_entries.size()):
		var entry: Dictionary = strain_entries[strain_index]
		var archive_id: String = String(entry["archive_id"])
		founders_by_strain[archive_id] = []
		for _local_index in range(int(initial_counts[strain_index])):
			var position: Vector2 = _competition_position(config, position_index, requested_founders)
			if position_index < positions.size():
				position = positions[position_index]
			var cell = _seed_entry(sim, entry, position)
			founder_to_strain[int(cell.id)] = archive_id
			founders_by_strain[archive_id].append(int(cell.id))
			descendants_by_founder[int(cell.id)] = 0
			position_index += 1
	sim.relax_mechanics()
	var initial_resources: Dictionary = _field_totals(sim)
	for _tick in range(horizon_ticks):
		if not environment.is_empty():
			EnvironmentScheduleScript.apply(sim, int(sim.tick_index), environment)
		sim.step(1)
		if sim.cells.is_empty():
			break

	var parent_map: Dictionary = _parent_map(sim.event_log)
	var descendants: Dictionary = {}
	var division_counts: Dictionary = {}
	var death_causes: Dictionary = {}
	var fluxes: Dictionary = {}
	var spatial: Dictionary = {}
	for entry_variant in strain_entries:
		var archive_id: String = String(entry_variant["archive_id"])
		descendants[archive_id] = 0
		division_counts[archive_id] = 0
		death_causes[archive_id] = {}
		fluxes[archive_id] = {}
		spatial[archive_id] = {"count": 0, "min_x": INF, "max_x": -INF, "min_y": INF, "max_y": -INF}

	for cell in sim.cells:
		var root_id: int = _root_founder(int(cell.id), parent_map)
		var strain_id: String = String(founder_to_strain.get(root_id, ""))
		if strain_id.is_empty():
			continue
		descendants[strain_id] = int(descendants[strain_id]) + 1
		descendants_by_founder[root_id] = int(descendants_by_founder.get(root_id, 0)) + 1
		var box: Dictionary = spatial[strain_id]
		box["count"] = int(box["count"]) + 1
		box["min_x"] = minf(float(box["min_x"]), float(cell.position.x))
		box["max_x"] = maxf(float(box["max_x"]), float(cell.position.x))
		box["min_y"] = minf(float(box["min_y"]), float(cell.position.y))
		box["max_y"] = maxf(float(box["max_y"]), float(cell.position.y))
		spatial[strain_id] = box
		var strain_fluxes: Dictionary = fluxes[strain_id]
		for reaction_id_variant in cell.last_fluxes.keys():
			var reaction_id: String = String(reaction_id_variant)
			strain_fluxes[reaction_id] = float(strain_fluxes.get(reaction_id, 0.0)) + float(cell.last_fluxes[reaction_id_variant])
		fluxes[strain_id] = strain_fluxes

	for event_variant in sim.event_log:
		var event: Dictionary = event_variant
		var kind: String = String(event.get("kind", ""))
		if kind == "division":
			var root_id: int = _root_founder(int(event["parent_id"]), parent_map)
			var strain_id: String = String(founder_to_strain.get(root_id, ""))
			if not strain_id.is_empty():
				division_counts[strain_id] = int(division_counts[strain_id]) + 1
		elif kind == "death":
			var death_root: int = _root_founder(int(event["cell_id"]), parent_map)
			var death_strain: String = String(founder_to_strain.get(death_root, ""))
			if not death_strain.is_empty():
				var causes: Dictionary = death_causes[death_strain]
				var reason: String = String(event.get("reason", "unknown"))
				causes[reason] = int(causes.get(reason, 0)) + 1
				death_causes[death_strain] = causes

	var frequencies: Dictionary = {}
	var active_strains: int = 0
	for archive_id_variant in descendants.keys():
		var archive_id: String = String(archive_id_variant)
		var count: int = int(descendants[archive_id])
		frequencies[archive_id] = 0.0 if sim.population_size() == 0 else float(count) / float(sim.population_size())
		if count > 0:
			active_strains += 1
		var box: Dictionary = spatial[archive_id]
		if int(box["count"]) == 0:
			spatial[archive_id] = {"count": 0, "min_x": null, "max_x": null, "min_y": null, "max_y": null}

	var outcome: String = "coexistence"
	if sim.population_size() == 0:
		outcome = "extinction"
	elif active_strains == 1:
		outcome = "fixation"
	return {
		"seed": seed,
		"planned_horizon_ticks": horizon_ticks,
		"realized_ticks": int(sim.tick_index),
		"final_population": sim.population_size(),
		"outcome": outcome,
		"founders_by_strain": founders_by_strain,
		"descendants_by_founder": descendants_by_founder,
		"descendants_by_strain": descendants,
		"frequencies": frequencies,
		"division_events_by_strain": division_counts,
		"death_causes_by_strain": death_causes,
		"endpoint_fluxes_by_strain": fluxes,
		"spatial_occupation": spatial,
		"initial_resources": initial_resources,
		"final_resources": _field_totals(sim),
		"final_checksum": float(sim.checksum()),
		"timeline": semantic_timeline(sim.event_log)
	}

static func run_competition_batch(
	base_config,
	strain_entries: Array,
	initial_counts: Array,
	environment: Dictionary,
	horizon_ticks: int,
	seeds: Array,
	positions: Array = []
) -> Dictionary:
	assert(not seeds.is_empty())
	var runs: Array = []
	var by_seed: Dictionary = {}
	var outcome_counts: Dictionary = {"coexistence": 0, "fixation": 0, "extinction": 0}
	for seed_variant in seeds:
		var seed: int = int(seed_variant)
		var result: Dictionary = run_competition(
			base_config, strain_entries, initial_counts, environment, horizon_ticks, seed, positions
		)
		runs.append(result)
		by_seed[seed] = result
		var outcome: String = String(result["outcome"])
		outcome_counts[outcome] = int(outcome_counts.get(outcome, 0)) + 1
	return {
		"seeds": seeds.duplicate(),
		"runs": runs,
		"by_seed": by_seed,
		"outcome_counts": outcome_counts
	}

static func explain_genetic_candidate(reference_genome, candidate_genome, reference_phenotype: Dictionary = {}, candidate_phenotype: Dictionary = {}) -> Dictionary:
	var differences: Array = GeneticAssayScript.differences(reference_genome, candidate_genome)
	return {
		"genotype_differences": differences,
		"single_reversion_candidates": GeneticAssayScript.single_reversion_candidates(reference_genome, candidate_genome),
		"phenotype_delta": _phenotype_delta(reference_phenotype, candidate_phenotype),
		"requires_epistasis_test": differences.size() > 1,
		"interpretation": "candidate_only_until_controlled_assays"
	}

static func _seed_entry(sim, entry: Dictionary, position: Vector2):
	var id: int = int(sim.next_cell_id)
	sim.next_cell_id += 1
	var cell
	if String(entry["mode"]) == "molecular_state" and entry.has("cell_state"):
		cell = SnapshotCodecScript.restore_cell(entry["cell_state"])
		cell.id = id
		cell.parent_id = -1
		cell.generation = 0
		cell.birth_tick = int(sim.tick_index)
		cell.alive = true
		cell.death_reason = ""
	else:
		cell = CellStateScript.new(id, -1, 0, int(sim.tick_index), position, float(sim.config.ancestor_volume))
		cell.genome = SnapshotCodecScript.restore_genome(entry["genome"])
		cell.initialize_molecular_state(sim.config)
	cell.position = CellMechanicsScript.clamp_position(
		sim.world.clamp_position(position), CellMechanicsScript.radius_for_cell(cell, sim.config), sim.world
	)
	sim.cells.append(cell)
	sim._record_event("experimental_inoculation", {
		"cell_id": int(cell.id),
		"archive_id": String(entry["archive_id"]),
		"archive_mode": String(entry["mode"]),
		"intervention": true,
		"genotype_fingerprint": int(cell.genome.fingerprint())
	})
	return cell

static func _parent_map(event_log: Array) -> Dictionary:
	var result: Dictionary = {}
	for event_variant in event_log:
		var event: Dictionary = event_variant
		if String(event.get("kind", "")) == "birth":
			result[int(event["cell_id"])] = int(event.get("parent_id", -1))
		elif String(event.get("kind", "")) == "experimental_inoculation":
			result[int(event["cell_id"])] = -1
	return result

static func _root_founder(cell_id: int, parent_map: Dictionary) -> int:
	var current: int = cell_id
	var guard: int = 0
	while parent_map.has(current) and int(parent_map[current]) >= 0:
		current = int(parent_map[current])
		guard += 1
		assert(guard < 100000, "Ancestry loop detected")
	return current

static func _competition_position(config, index: int, total: int) -> Vector2:
	var center := Vector2(float(config.world_width - 1) * 0.5, float(config.world_height - 1) * 0.5)
	var spacing: float = maxf(1.1, float(config.ancestor_radius_grid) * 2.4)
	var offset: float = (float(index) - float(total - 1) * 0.5) * spacing
	return Vector2(center.x + offset, center.y)

static func _field_totals(sim) -> Dictionary:
	var result: Dictionary = {}
	for field_name in sim.world.field_order:
		result[field_name] = float(sim.world.get_field(field_name).total_amount())
	return result

static func _phenotype_delta(reference: Dictionary, candidate: Dictionary) -> Dictionary:
	var result: Dictionary = {}
	for key_variant in candidate.keys():
		var key: String = String(key_variant)
		if not reference.has(key):
			continue
		if (candidate[key_variant] is int or candidate[key_variant] is float) and (reference[key] is int or reference[key] is float):
			result[key] = float(candidate[key_variant]) - float(reference[key])
	return result
