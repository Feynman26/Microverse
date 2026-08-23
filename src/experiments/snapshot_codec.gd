extends RefCounted
class_name SnapshotCodec

const SimConfigScript = preload("res://src/core/sim_config.gd")
const SimulationEngineScript = preload("res://src/simulation/simulation_engine.gd")
const CellStateScript = preload("res://src/biology/cell_state.gd")
const GeneScript = preload("res://src/genetics/gene.gd")
const GenomeScript = preload("res://src/genetics/genome.gd")

const SNAPSHOT_SCHEMA_VERSION: int = 1
const MODEL_IDENTIFIER: String = "microverse-m9"

const CONFIG_FIELDS: Array[String] = [
	"tick_dt_min", "world_width", "world_height", "grid_cell_size_um", "max_cells", "seed",
	"initial_glucose", "initial_oxygen", "initial_nitrogen", "initial_phosphorus",
	"glucose_diffusion", "oxygen_diffusion", "nitrogen_diffusion", "phosphorus_diffusion",
	"secondary_extracellular_initial", "secondary_extracellular_diffusion",
	"glucose_transport_vmax", "glucose_transport_km", "oxygen_transport_vmax", "oxygen_transport_km",
	"nitrogen_transport_vmax", "nitrogen_transport_km", "phosphorus_transport_vmax", "phosphorus_transport_km",
	"intracellular_pool_capacity_per_volume", "secondary_transport_vmax_per_reference_protein",
	"secondary_transport_gradient_km", "secondary_transport_atp_cost_per_unit",
	"extracellular_protein_diffusion", "extracellular_protein_secretion_fraction_per_min",
	"extracellular_protein_secretion_atp_cost_per_unit", "extracellular_catalysis_rate_scale",
	"extracellular_catalysis_km", "metabolic_substeps_per_tick", "metabolic_km_per_volume",
	"metabolic_rate_scale", "biomass_units_per_volume", "initial_atp_per_volume", "initial_adp_per_volume",
	"initial_nad_per_volume", "initial_nadh_per_volume", "transcription_max_events_per_min",
	"mrna_decay_rate_per_min", "translation_events_per_mrna_per_min", "protein_decay_rate_per_min",
	"expression_reference_protein_count", "transcription_atp_cost_per_event", "transcription_nuc_cost_per_event",
	"translation_atp_cost_per_event", "translation_aa_cost_per_event", "expression_partition_noise_scale",
	"proteome_capacity_reference_units", "translation_capacity_fraction_of_proteome_per_min",
	"regulation_enabled", "regulatory_max_distance", "regulatory_distance_decay", "regulatory_gain",
	"regulatory_min_factor", "regulatory_max_factor", "allostery_enabled", "allosteric_max_distance",
	"allosteric_distance_decay", "allosteric_km", "allosteric_gain", "allosteric_min_factor",
	"allosteric_max_factor", "maintenance_atp_rate_per_volume", "spontaneous_ros_decay_rate",
	"ros_damage_rate", "basal_repair_rate", "repair_atp_cost", "lethal_damage", "lethal_energy_debt",
	"ancestor_volume", "division_volume", "division_atp_cost", "partition_jitter", "daughter_offset_grid",
	"ancestor_radius_grid", "mechanical_relaxation_iterations", "mechanical_relaxation_fraction",
	"mechanical_overlap_tolerance", "mechanical_use_spatial_index", "mechanical_neighbor_bucket_size_grid",
	"mutation_enabled", "promoter_mutation_rate_per_gene", "signature_mutation_rate_per_gene",
	"regulatory_signature_mutation_rate_per_gene", "neutral_marker_mutation_rate_per_gene",
	"promoter_mutation_step_max", "protein_signature_bits"
]

static func capture(sim, experiment_context: Dictionary = {}) -> Dictionary:
	var snapshot: Dictionary = {
		"schema_version": SNAPSHOT_SCHEMA_VERSION,
		"model_identifier": MODEL_IDENTIFIER,
		"config": capture_config(sim.config),
		"tick_index": int(sim.tick_index),
		"simulation_time_min": float(sim.simulation_time_min),
		"rng_state": int(sim.rng.get_state()),
		"world": _capture_world(sim.world),
		"cells": [],
		"next_cell_id": int(sim.next_cell_id),
		"next_mutation_id": int(sim.next_mutation_id),
		"event_log": sim.event_log.duplicate(true),
		"last_mechanics_summary": sim.last_mechanics_summary.duplicate(true),
		"last_secondary_transport_summary": sim.last_secondary_transport_summary.duplicate(true),
		"last_protein_secretion_summary": sim.last_protein_secretion_summary.duplicate(true),
		"last_extracellular_catalysis_summary": sim.last_extracellular_catalysis_summary.duplicate(true),
		"scheduled_interventions": experiment_context.get("scheduled_interventions", []).duplicate(true),
		"environment_boundary_state": _clone_variant(experiment_context.get("environment_boundary_state", {})),
		"source_checksum": float(sim.checksum())
	}
	for cell in sim.cells:
		snapshot["cells"].append(capture_cell(cell))
	snapshot["fingerprint"] = fingerprint(snapshot)
	return snapshot

static func capture_config(config) -> Dictionary:
	var result: Dictionary = {}
	for field_name in CONFIG_FIELDS:
		result[field_name] = _clone_variant(config.get(field_name))
	return result

static func capture_cell(cell) -> Dictionary:
	return {
		"id": int(cell.id),
		"parent_id": int(cell.parent_id),
		"generation": int(cell.generation),
		"birth_tick": int(cell.birth_tick),
		"position": cell.position,
		"alive": bool(cell.alive),
		"death_reason": String(cell.death_reason),
		"genome": capture_genome(cell.genome),
		"metabolites": cell.metabolites.duplicate(true),
		"expression_state": cell.expression_state.duplicate(true),
		"last_fluxes": cell.last_fluxes.duplicate(true),
		"last_expression_summary": cell.last_expression_summary.duplicate(true),
		"volume": float(cell.volume),
		"damage": float(cell.damage),
		"energy_debt": float(cell.energy_debt)
	}

static func capture_genome(genome) -> Dictionary:
	assert(genome != null)
	var genes: Array = []
	for gene in genome.genes:
		genes.append({
			"locus_id": int(gene.locus_id),
			"promoter_code": int(gene.promoter_code),
			"protein_signature": int(gene.protein_signature),
			"neutral_marker": int(gene.neutral_marker),
			"regulatory_signature": int(gene.regulatory_signature)
		})
	return {"genes": genes, "canonical_key": genome.canonical_key(), "fingerprint": int(genome.fingerprint())}

static func restore(snapshot: Dictionary):
	_validate_snapshot(snapshot)
	var config = restore_config(snapshot["config"])
	var sim = SimulationEngineScript.new(config)
	sim.tick_index = int(snapshot["tick_index"])
	sim.simulation_time_min = float(snapshot["simulation_time_min"])
	sim.rng.set_state(int(snapshot["rng_state"]))
	_restore_world(sim.world, snapshot["world"])
	sim.cells.clear()
	for cell_data_variant in snapshot["cells"]:
		sim.cells.append(restore_cell(cell_data_variant))
	sim.next_cell_id = int(snapshot["next_cell_id"])
	sim.next_mutation_id = int(snapshot["next_mutation_id"])
	sim.event_log = snapshot["event_log"].duplicate(true)
	sim.last_mechanics_summary = snapshot["last_mechanics_summary"].duplicate(true)
	sim.last_secondary_transport_summary = snapshot["last_secondary_transport_summary"].duplicate(true)
	sim.last_protein_secretion_summary = snapshot["last_protein_secretion_summary"].duplicate(true)
	sim.last_extracellular_catalysis_summary = snapshot["last_extracellular_catalysis_summary"].duplicate(true)
	return sim

static func restore_config(data: Dictionary):
	var config = SimConfigScript.new()
	for field_name in CONFIG_FIELDS:
		assert(data.has(field_name), "Snapshot config missing field: %s" % field_name)
		config.set(field_name, _clone_variant(data[field_name]))
	config.validate()
	return config

static func restore_cell(data: Dictionary):
	var cell = CellStateScript.new(
		int(data["id"]), int(data["parent_id"]), int(data["generation"]), int(data["birth_tick"]),
		data["position"], float(data["volume"])
	)
	cell.alive = bool(data["alive"])
	cell.death_reason = String(data["death_reason"])
	cell.genome = restore_genome(data["genome"])
	cell.metabolites = data["metabolites"].duplicate(true)
	cell.expression_state = data["expression_state"].duplicate(true)
	cell.last_fluxes = data["last_fluxes"].duplicate(true)
	cell.last_expression_summary = data["last_expression_summary"].duplicate(true)
	cell.volume = float(data["volume"])
	cell.damage = float(data["damage"])
	cell.energy_debt = float(data["energy_debt"])
	return cell

static func restore_genome(data: Dictionary):
	var genes: Array = []
	for gene_data_variant in data["genes"]:
		var gene_data: Dictionary = gene_data_variant
		genes.append(GeneScript.new(
			int(gene_data["locus_id"]), int(gene_data["promoter_code"]), int(gene_data["protein_signature"]),
			int(gene_data["neutral_marker"]), int(gene_data["regulatory_signature"])
		))
	var genome = GenomeScript.new(genes)
	assert(genome.canonical_key() == String(data["canonical_key"]))
	return genome

static func fork(sim, experiment_context: Dictionary = {}) -> Array:
	var snapshot: Dictionary = capture(sim, experiment_context)
	return [restore(snapshot), restore(snapshot)]

static func encode(snapshot: Dictionary) -> PackedByteArray:
	_validate_snapshot(snapshot)
	return var_to_bytes(snapshot)

static func decode(bytes: PackedByteArray) -> Dictionary:
	var decoded = bytes_to_var(bytes)
	assert(decoded is Dictionary)
	_validate_snapshot(decoded)
	return decoded

static func fingerprint(snapshot: Dictionary) -> String:
	var payload: Dictionary = snapshot.duplicate(true)
	payload.erase("fingerprint")
	var context := HashingContext.new()
	var error: Error = context.start(HashingContext.HASH_SHA256)
	assert(error == OK)
	context.update(var_to_bytes(_canonical(payload)))
	return context.finish().hex_encode()

static func _capture_world(world) -> Dictionary:
	var fields: Array = []
	for field_name in world.field_order:
		var field = world.get_field(field_name)
		fields.append({
			"name": field_name,
			"diffusion_coefficient": float(field.diffusion_coefficient),
			"values": field.values.duplicate()
		})
	var protein_fields: Array = []
	for signature_variant in world.protein_signatures():
		var signature: int = int(signature_variant)
		var protein_field = world.get_protein_field(signature)
		protein_fields.append({
			"signature": signature,
			"diffusion_coefficient": float(protein_field.diffusion_coefficient),
			"values": protein_field.values.duplicate()
		})
	return {
		"width": int(world.width),
		"height": int(world.height),
		"cell_size": float(world.cell_size),
		"field_order": world.field_order.duplicate(),
		"fields": fields,
		"protein_fields": protein_fields
	}

static func _restore_world(world, data: Dictionary) -> void:
	assert(int(data["width"]) == int(world.width))
	assert(int(data["height"]) == int(world.height))
	assert(absf(float(data["cell_size"]) - float(world.cell_size)) <= 1e-12)
	assert(data["field_order"] == world.field_order)
	for field_data_variant in data["fields"]:
		var field_data: Dictionary = field_data_variant
		var field = world.get_field(String(field_data["name"]))
		field.diffusion_coefficient = float(field_data["diffusion_coefficient"])
		field.replace_values(field_data["values"])
	world.protein_fields.clear()
	for protein_data_variant in data["protein_fields"]:
		var protein_data: Dictionary = protein_data_variant
		var protein_field = world.ensure_protein_field(int(protein_data["signature"]), float(protein_data["diffusion_coefficient"]))
		protein_field.replace_values(protein_data["values"])
	world.assert_nonnegative()

static func _validate_snapshot(snapshot: Dictionary) -> void:
	assert(int(snapshot.get("schema_version", -1)) == SNAPSHOT_SCHEMA_VERSION)
	assert(String(snapshot.get("model_identifier", "")) == MODEL_IDENTIFIER)
	for required_key in ["config", "tick_index", "simulation_time_min", "rng_state", "world", "cells", "next_cell_id", "next_mutation_id", "event_log"]:
		assert(snapshot.has(required_key), "Snapshot missing authoritative key: %s" % required_key)

static func _clone_variant(value):
	if value is Dictionary or value is Array:
		return value.duplicate(true)
	if value is PackedFloat64Array:
		return value.duplicate()
	return value

static func _canonical(value):
	if value is Dictionary:
		var keys: Array = value.keys()
		keys.sort()
		var result: Dictionary = {}
		for key in keys:
			result[key] = _canonical(value[key])
		return result
	if value is Array:
		var result_array: Array = []
		for item in value:
			result_array.append(_canonical(item))
		return result_array
	if value is PackedFloat64Array:
		var floats: Array = []
		for item in value:
			floats.append(item)
		return floats
	return value
