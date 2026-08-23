extends RefCounted
class_name SimConfig

# Time and world scale.
var tick_dt_min: float = 0.10
var world_width: int = 64
var world_height: int = 64
var grid_cell_size_um: float = 1.0
var max_cells: int = 64
var seed: int = 824718

# Environment concentrations.
var initial_glucose: float = 4.0
var initial_oxygen: float = 5.0
var initial_nitrogen: float = 3.0
var initial_phosphorus: float = 2.0
var glucose_diffusion: float = 0.80
var oxygen_diffusion: float = 1.00
var nitrogen_diffusion: float = 0.70
var phosphorus_diffusion: float = 0.50

# Transitional basal membrane transport. Intracellular catalytic capacity is
# genome/proteome-derived from M4 onward; membrane transport becomes similarly
# evolvable after the expression/regulation layer is stable.
var glucose_transport_vmax: float = 0.60
var glucose_transport_km: float = 0.80
var oxygen_transport_vmax: float = 0.90
var oxygen_transport_km: float = 0.60
var nitrogen_transport_vmax: float = 0.45
var nitrogen_transport_km: float = 0.50
var phosphorus_transport_vmax: float = 0.30
var phosphorus_transport_km: float = 0.40
var intracellular_pool_capacity_per_volume: float = 8.0

# M4 intracellular chemistry.
var metabolic_substeps_per_tick: int = 6
var metabolic_km_per_volume: float = 0.20
var metabolic_rate_scale: float = 0.85
var biomass_units_per_volume: float = 1.0
var initial_atp_per_volume: float = 2.0
var initial_adp_per_volume: float = 8.0
var initial_nad_per_volume: float = 6.0
var initial_nadh_per_volume: float = 1.0

# M5-A explicit expression kinetics. Values are compressed digital kinetics, not
# quantitative bacterial copy-number claims. At basal steady state the chosen
# rates give approximately `expression_reference_protein_count * promoter_strength`
# proteins per locus, preserving the M4 phenotype scale while replacing the
# shortcut with explicit state and stochastic turnover.
var transcription_max_events_per_min: float = 1.0
var mrna_decay_rate_per_min: float = 0.25
var translation_events_per_mrna_per_min: float = 2.0
var protein_decay_rate_per_min: float = 0.05
var expression_reference_protein_count: float = 160.0

# Expression resource costs per accepted molecular event. Material withdrawn
# from NUC/AA is stored in mRNA/protein state and fully returned on degradation
# in M5-A (compressed recycling). ATP is converted to ADP.
var transcription_atp_cost_per_event: float = 0.010
var transcription_nuc_cost_per_event: float = 0.002
var translation_atp_cost_per_event: float = 0.003
var translation_aa_cost_per_event: float = 0.001
var expression_partition_noise_scale: float = 0.12

# Basal physiology that is not yet part of the evolvable molecular grammar.
var maintenance_atp_rate_per_volume: float = 0.08
var spontaneous_ros_decay_rate: float = 0.01
var ros_damage_rate: float = 0.04
var basal_repair_rate: float = 0.02
var repair_atp_cost: float = 0.20
var lethal_damage: float = 3.0
var lethal_energy_debt: float = 4.0

# Cell cycle. Volume is derived from BIO / biomass_units_per_volume.
var ancestor_volume: float = 1.0
var division_volume: float = 2.0
var division_atp_cost: float = 1.0
var partition_jitter: float = 0.02
var daughter_offset_grid: float = 0.20

# M3 genetics. These rates are per gene per replication event. They remain
# externally configured until M10 derives mutation rate from evolved
# replication/repair machinery.
var mutation_enabled: bool = true
var promoter_mutation_rate_per_gene: float = 0.001
var signature_mutation_rate_per_gene: float = 0.001
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
	var alpha_g: float = glucose_diffusion * tick_dt_min / dx2
	var alpha_o: float = oxygen_diffusion * tick_dt_min / dx2
	var alpha_n: float = nitrogen_diffusion * tick_dt_min / dx2
	var alpha_p: float = phosphorus_diffusion * tick_dt_min / dx2
	assert(alpha_g <= 0.25, "Explicit 2D diffusion unstable for glucose: alpha must be <= 0.25")
	assert(alpha_o <= 0.25, "Explicit 2D diffusion unstable for oxygen: alpha must be <= 0.25")
	assert(alpha_n <= 0.25, "Explicit 2D diffusion unstable for nitrogen: alpha must be <= 0.25")
	assert(alpha_p <= 0.25, "Explicit 2D diffusion unstable for phosphorus: alpha must be <= 0.25")
	assert(glucose_transport_vmax >= 0.0 and oxygen_transport_vmax >= 0.0)
	assert(nitrogen_transport_vmax >= 0.0 and phosphorus_transport_vmax >= 0.0)
	assert(metabolic_substeps_per_tick >= 1)
	assert(metabolic_km_per_volume > 0.0)
	assert(metabolic_rate_scale > 0.0)
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
	assert(division_volume > ancestor_volume)
	assert(partition_jitter >= 0.0 and partition_jitter < 0.5)
	assert(promoter_mutation_rate_per_gene >= 0.0 and promoter_mutation_rate_per_gene <= 1.0)
	assert(signature_mutation_rate_per_gene >= 0.0 and signature_mutation_rate_per_gene <= 1.0)
	assert(neutral_marker_mutation_rate_per_gene >= 0.0 and neutral_marker_mutation_rate_per_gene <= 1.0)
	assert(promoter_mutation_step_max >= 1)
	assert(protein_signature_bits >= 1 and protein_signature_bits <= 16)
