extends RefCounted
class_name StrainLibrary

const SnapshotCodecScript = preload("res://src/experiments/snapshot_codec.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const CellMechanicsScript = preload("res://src/physics/cell_mechanics.gd")

const MODE_GENOTYPE_ONLY: String = "genotype_only"
const MODE_MOLECULAR_STATE: String = "molecular_state"

var entries: Dictionary = {}
var next_archive_id: int = 1

func freeze_genotype(sim, cell, label: String = "", notes: String = "") -> Dictionary:
	return _freeze(sim, cell, MODE_GENOTYPE_ONLY, label, notes)

func freeze_molecular_state(sim, cell, label: String = "", notes: String = "") -> Dictionary:
	return _freeze(sim, cell, MODE_MOLECULAR_STATE, label, notes)

func get_entry(archive_id: String) -> Dictionary:
	assert(entries.has(archive_id), "Unknown frozen strain: %s" % archive_id)
	return entries[archive_id].duplicate(true)

func all_entries() -> Array:
	var ids: Array = entries.keys()
	ids.sort()
	var result: Array = []
	for archive_id_variant in ids:
		result.append(entries[String(archive_id_variant)].duplicate(true))
	return result

func reintroduce(sim, archive_id: String, position: Vector2):
	var entry: Dictionary = get_entry(archive_id)
	var cell
	if String(entry["mode"]) == MODE_GENOTYPE_ONLY:
		cell = CellStateScript.new(
			_allocate_cell_id(sim), -1, 0, int(sim.tick_index), position, float(sim.config.ancestor_volume)
		)
		cell.genome = SnapshotCodecScript.restore_genome(entry["genome"])
		cell.initialize_molecular_state(sim.config)
	else:
		assert(String(entry["mode"]) == MODE_MOLECULAR_STATE)
		cell = SnapshotCodecScript.restore_cell(entry["cell_state"])
		cell.id = _allocate_cell_id(sim)
		cell.parent_id = -1
		cell.generation = 0
		cell.birth_tick = int(sim.tick_index)
		cell.position = position
		cell.alive = true
		cell.death_reason = ""
	cell.position = CellMechanicsScript.clamp_position(
		sim.world.clamp_position(cell.position), CellMechanicsScript.radius_for_cell(cell, sim.config), sim.world
	)
	sim.cells.append(cell)
	sim._record_event("experimental_inoculation", {
		"cell_id": int(cell.id),
		"archive_id": archive_id,
		"archive_mode": String(entry["mode"]),
		"genotype_fingerprint": int(cell.genome.fingerprint()),
		"intervention": true
	})
	return cell

func _freeze(sim, cell, mode: String, label: String, notes: String) -> Dictionary:
	assert(cell != null and cell.alive and cell.genome != null)
	assert(mode == MODE_GENOTYPE_ONLY or mode == MODE_MOLECULAR_STATE)
	var archive_id: String = "strain-%06d" % next_archive_id
	next_archive_id += 1
	var entry: Dictionary = {
		"archive_id": archive_id,
		"mode": mode,
		"label": label,
		"notes": notes,
		"source_seed": int(sim.config.seed),
		"source_tick": int(sim.tick_index),
		"source_time_min": float(sim.simulation_time_min),
		"source_cell_id": int(cell.id),
		"source_parent_id": int(cell.parent_id),
		"source_generation": int(cell.generation),
		"source_world_checksum": float(sim.world.checksum()),
		"genome": SnapshotCodecScript.capture_genome(cell.genome),
		"phenotype_summary": _phenotype_summary(cell)
	}
	if mode == MODE_MOLECULAR_STATE:
		entry["cell_state"] = SnapshotCodecScript.capture_cell(cell)
	entries[archive_id] = entry.duplicate(true)
	return entry.duplicate(true)

func _phenotype_summary(cell) -> Dictionary:
	return {
		"volume": float(cell.volume),
		"damage": float(cell.damage),
		"energy_debt": float(cell.energy_debt),
		"metabolites": cell.metabolites.duplicate(true),
		"last_fluxes": cell.last_fluxes.duplicate(true),
		"total_mrna": float(cell.total_mrna()),
		"total_protein": float(cell.total_protein())
	}

func _allocate_cell_id(sim) -> int:
	var result: int = int(sim.next_cell_id)
	sim.next_cell_id += 1
	return result
