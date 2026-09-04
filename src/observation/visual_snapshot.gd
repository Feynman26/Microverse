extends RefCounted
class_name VisualSnapshot

const DNAReplicationScript = preload("res://src/genetics/dna_replication.gd")

# Produces an immutable-by-contract presentation view. No authoritative object
# reference is published: the UI receives only scalar values and fresh packed
# arrays/dictionaries that the worker never mutates after publication.

static func capture(simulation) -> Dictionary:
	var glucose_field = simulation.world.get_field("glucose")
	var positions := PackedVector2Array()
	var volumes := PackedFloat32Array()
	positions.resize(simulation.cells.size())
	volumes.resize(simulation.cells.size())
	for index in range(simulation.cells.size()):
		var cell = simulation.cells[index]
		positions[index] = cell.position
		volumes[index] = float(cell.volume)

	return {
		"tick": int(simulation.tick_index),
		"simulation_time_min": float(simulation.simulation_time_min),
		"population": int(simulation.population_size()),
		"max_cells": int(simulation.config.max_cells),
		"maximum_generation": int(simulation.maximum_generation()),
		"genotype_count": int(simulation.genotype_count()),
		"mutation_event_count": int(simulation.mutation_event_count()),
		"total_cell_volume": float(simulation.total_cell_volume()),
		"event_count": int(simulation.event_log.size()),
		"seed": int(simulation.config.seed),
		"checksum": float(simulation.checksum()),
		"world_width": int(simulation.config.world_width),
		"world_height": int(simulation.config.world_height),
		"initial_glucose": float(simulation.config.initial_glucose),
		"glucose_values": glucose_field.values.duplicate(),
		"cell_positions": positions,
		"cell_volumes": volumes,
		"environment": _environment(simulation),
		"division": _division_status(simulation),
		"focal_cell": _focal_cell(simulation)
	}

static func _environment(simulation) -> Dictionary:
	return {
		"glucose": float(simulation.world.get_field("glucose").total_amount()),
		"oxygen": float(simulation.world.get_field("oxygen").total_amount()),
		"nitrogen": float(simulation.world.get_field("nitrogen").total_amount()),
		"phosphorus": float(simulation.world.get_field("phosphorus").total_amount())
	}

static func _division_status(simulation) -> Dictionary:
	var ready: int = 0
	var volume_ready: int = 0
	var volume_ready_atp_blocked: int = 0
	var volume_atp_ready_replication_blocked: int = 0
	for cell in simulation.cells:
		if not cell.alive:
			continue
		var has_volume: bool = float(cell.volume) >= float(simulation.config.division_volume)
		var has_atp: bool = float(cell.pool("ATP")) >= float(simulation.config.division_atp_cost)
		var has_replication: bool = (
			not bool(simulation.config.evolvable_replication_enabled)
			or DNAReplicationScript.replication_complete(cell)
		)
		if has_volume:
			volume_ready += 1
			if not has_atp:
				volume_ready_atp_blocked += 1
			elif not has_replication:
				volume_atp_ready_replication_blocked += 1
			else:
				ready += 1
	return {
		"ready": ready,
		"volume_ready": volume_ready,
		"volume_ready_atp_blocked": volume_ready_atp_blocked,
		"volume_atp_ready_replication_blocked": volume_atp_ready_replication_blocked
	}

static func _focal_cell(simulation) -> Dictionary:
	if simulation.cells.is_empty():
		return {}
	var cell = simulation.cells[0]
	var expression: Dictionary = cell.last_expression_summary
	return {
		"id": int(cell.id),
		"generation": int(cell.generation),
		"genotype": int(cell.genome.fingerprint()),
		"volume": float(cell.volume),
		"bio": float(cell.pool("BIO")),
		"atp": float(cell.pool("ATP")),
		"adp": float(cell.pool("ADP")),
		"adenylate": float(cell.total_adenylate()),
		"redox_currency": float(cell.total_redox_currency()),
		"replication_progress": float(cell.replication_progress),
		"g": float(cell.pool("G")),
		"c3": float(cell.pool("C3")),
		"c2": float(cell.pool("C2")),
		"w1": float(cell.pool("W1")),
		"w2": float(cell.pool("W2")),
		"ros": float(cell.pool("ROS")),
		"damage": float(cell.damage),
		"mrna": float(cell.total_mrna()),
		"protein": float(cell.total_protein()),
		"protein_cohorts": _protein_cohort_count(cell),
		"transcribed": float(expression.get("transcribed", 0.0)),
		"translated": float(expression.get("translated", 0.0)),
		"expression_atp": float(expression.get("atp_spent", 0.0)),
		"dominant_flux": _dominant_flux(cell)
	}

static func _protein_cohort_count(cell) -> int:
	var count: int = 0
	for locus_state in cell.expression_state.values():
		count += locus_state["protein"].size()
	return count

static func _dominant_flux(cell) -> Dictionary:
	var best_id: String = "none"
	var best_flux: float = 0.0
	for reaction_id in cell.last_fluxes.keys():
		var flux: float = float(cell.last_fluxes[reaction_id])
		if flux > best_flux:
			best_flux = flux
			best_id = String(reaction_id)
	return {"reaction_id": best_id, "flux": best_flux}
