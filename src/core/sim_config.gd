extends RefCounted
class_name SimConfig

var tick_dt_min: float = 0.10
var world_width: int = 64
var world_height: int = 64
var grid_cell_size_um: float = 1.0
var max_cells: int = 64
var seed: int = 824718

var initial_glucose: float = 4.0
var initial_oxygen: float = 5.0
var initial_nitrogen: float = 3.0
var initial_phosphorus: float = 2.0
var glucose_diffusion: float = 0.80
var oxygen_diffusion: float = 1.00
var nitrogen_diffusion: float = 0.70
var phosphorus_diffusion: float = 0.50

var glucose_transport_vmax: float = 0.60
var glucose_transport_km: float = 0.80
var oxygen_transport_vmax: float = 0.90
var oxygen_transport_km: float = 0.60
var nitrogen_transport_vmax: float = 0.45
var nitrogen_transport_km: float = 0.50
var phosphorus_transport_vmax: float = 0.30
var phosphorus_transport_km: float = 0.40
var intracellular_pool_capacity_per_volume: float = 8.0

var metabolic_substeps_per_tick: int = 6
var metabolic_km_per_volume: float = 0.20
var metabolic_rate_scale: float = 0.85
var biomass_units_per_volume: float = 1.0
var initial_atp_per_volume: float = 2.0
var initial_adp_per_volume: float = 8.0
var initial_nad_per_volume: float = 6.0
var initial_nadh_per_volume: float = 1.0

# M5-A canonical stochastic expression.
var transcription_max_events_per_min: float = 1.0
var mrna_decay_rate_per_min: float = 0.25
var translation_events_per_mrna_per_min: float = 2.0
var protein_decay_rate_per_min: float = 0.05
var expression_reference_protein_count: float = 160.0
var transcription_atp_cost_per_event: float = 0.010
var transcription_nuc_cost_per_event: float = 0.002
var translation_atp_cost_per_event: float = 0.003
var translation_aa_cost_per_event: float = 0.001
var expression_partition_noise_scale: float = 0.12

# M5-C finite proteome budget. The compressed cell can physically maintain only
# this many reference-protein equivalents at once. Expression remains free to
# propose any proteins, but excess realized protein is removed proportionally
# and its amino-acid material is recycled. ATP already spent on unnecessary
# synthesis is not refunded, creating a generic expression opportunity cost.
# This is a global physical constraint, not a genotype-specific fitness trait.
var proteome_capacity_reference_units: float = 5.0

# Shared translation machinery. All mRNA cohorts compete simultaneously for a
# finite ribosomal throughput proportional to the maximum proteome. At the
# default 800-protein capacity, 0.06 corresponds to 48 translation events/min.
# Basal protein turnover is ~40 proteins/min, leaving limited spare throughput
# for environmental remodeling. No locus or protein identity receives priority.
var translation_capacity_fraction_of_proteome_per_min: float = 0.06

# M5-B generic promoter regulation. Any physical protein sequence cohort can
# bind any promoter motif if its 16-bit signature is close enough. No gene is
# tagged as a named transcription factor or behavior controller.
var regulation_enabled: bool = true
var regulatory_max_distance: int = 3
var regulatory_distance_decay: float = 0.80
var regulatory_gain: float = 0.60
var regulatory_min_factor: float = 0.25
var regulatory_max_factor: float = 1.75

# M5-B generic chemical sensing. Metabolites carry ligand signatures; compatible
# protein cohorts have their promoter-binding contribution allosterically
# potentiated/inhibited. Protein bit 14 selects direction; this is a digital
# chemistry convention, not a named sensor class.
var allostery_enabled: bool = true
var allosteric_max_distance: int = 2
var allosteric_distance_decay: float = 0.90
var allosteric_km: float = 0.50
var allosteric_gain: float = 0.75
var allosteric_min_factor: float = 0.25
var allosteric_max_factor: float = 2.00

var maintenance_atp_rate_per_volume: float = 0.08
var spontaneous_ros_decay_rate: float = 0.01
var ros_damage_rate: float = 0.04
var basal_repair_rate: float = 0.02
var repair_atp_cost: float = 0.20
var lethal_damage: float = 3.0
var lethal_energy_debt: float = 4.0

var ancestor_volume: float = 1.0
var division_volume: float = 2.0
var division_atp_cost: float = 1.0
var partition_jitter: float = 0.02
var daughter_offset_grid: float = 0.20

var mutation_enabled: bool = true
var promoter_mutation_rate_per_gene: float = 0.001
var signature_mutation_rate_per_gene: float = 0.001
var regulatory_signature_mutation_rate_per_gene: float = 0.001
var neutral_marker_mutation_rate_per_gene: float = 0.001
var promoter_mutation_step_max: int = 250
var protein_signature_bits: int = 16

func validate() -> void:
	assert(tick_dt_min > 0.0)
	assert(world_width > 2 and world_height > 2)
	assert(grid_cell_size_um > 0.0)
	assert(max_cells >= 1)
	assert(initial_glucose >= 0.0 and initial_oxygen >= 0.0)
	assert(initial_nitrogen >= 0.0 and initial_phosphorus >= 0.0)
	assert(glucose_diffusion >= 0.0 and oxygen_diffusion >= 0.0)
	assert(nitrogen_diffusion >= 0.0 and phosphorus_diffusion >= 0.0)
	var dx2: float = grid_cell_size_um * grid_cell_size_um
	assert(glucose_diffusion * tick_dt_min / dx2 <= 0.25)
	assert(oxygen_diffusion * tick_dt_min / dx2 <= 0.25)
	assert(nitrogen_diffusion * tick_dt_min / dx2 <= 0.25)
	assert(phosphorus_diffusion * tick_dt_min / dx2 <= 0.25)
	assert(glucose_transport_vmax >= 0.0 and oxygen_transport_vmax >= 0.0)
	assert(nitrogen_transport_vmax >= 0.0 and phosphorus_transport_vmax >= 0.0)
	assert(metabolic_substeps_per_tick >= 1)
	assert(metabolic_km_per_volume > 0.0 and metabolic_rate_scale > 0.0)
	assert(biomass_units_per_volume > 0.0)
	assert(initial_atp_per_volume >= 0.0 and initial_adp_per_volume >= 0.0)
	assert(initial_nad_per_volume >= 0.0 and initial_nadh_per_volume >= 0.0)
	assert(transcription_max_events_per_min >= 0.0)
	assert(mrna_decay_rate_per_min > 0.0)
	assert(translation_events_per_mrna_per_min >= 0.0)
	assert(protein_decay_rate_per_min > 0.0)
	assert(expression_reference_protein_count > 0.0)
	assert(transcription_atp_cost_per_event >= 0.0 and transcription_nuc_cost_per_event >= 0.0)
	assert(translation_atp_cost_per_event >= 0.0 and translation_aa_cost_per_event >= 0.0)
	assert(expression_partition_noise_scale >= 0.0 and expression_partition_noise_scale < 0.5)
	assert(proteome_capacity_reference_units > 0.0)
	assert(translation_capacity_fraction_of_proteome_per_min > 0.0 and translation_capacity_fraction_of_proteome_per_min <= 1.0)
	assert(regulatory_max_distance >= 0 and regulatory_max_distance <= 16)
	assert(regulatory_distance_decay >= 0.0 and regulatory_gain >= 0.0)
	assert(regulatory_min_factor >= 0.0 and regulatory_max_factor >= regulatory_min_factor)
	assert(allosteric_max_distance >= 0 and allosteric_max_distance <= 16)
	assert(allosteric_distance_decay >= 0.0 and allosteric_km > 0.0 and allosteric_gain >= 0.0)
	assert(allosteric_min_factor >= 0.0 and allosteric_max_factor >= allosteric_min_factor)
	assert(division_volume > ancestor_volume)
	assert(partition_jitter >= 0.0 and partition_jitter < 0.5)
	assert(promoter_mutation_rate_per_gene >= 0.0 and promoter_mutation_rate_per_gene <= 1.0)
	assert(signature_mutation_rate_per_gene >= 0.0 and signature_mutation_rate_per_gene <= 1.0)
	assert(regulatory_signature_mutation_rate_per_gene >= 0.0 and regulatory_signature_mutation_rate_per_gene <= 1.0)
	assert(neutral_marker_mutation_rate_per_gene >= 0.0 and neutral_marker_mutation_rate_per_gene <= 1.0)
	assert(promoter_mutation_step_max >= 1)
	assert(protein_signature_bits >= 1 and protein_signature_bits <= 16)