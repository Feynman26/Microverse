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

var transcription_rate: float = 0.08
var mrna_decay_rate: float = 0.20
var translation_rate: float = 0.10
var protein_decay_rate: float = 0.04
var transcription_atp_cost_per_unit: float = 0.20
var translation_atp_cost_per_unit: float = 0.35
var nuc_cost_per_mrna_unit: float = 0.015
var aa_cost_per_protein_unit: float = 0.020
var expression_noise_fraction: float = 0.05
var expression_partition_jitter: float = 0.03

var regulation_enabled: bool = true
var regulatory_max_distance: int = 3
var regulatory_distance_decay: float = 0.80
var regulatory_gain: float = 0.20
var regulatory_min_factor: float = 0.50
var regulatory_max_factor: float = 1.50

# M5 chemical sensing is generic allostery, not a named sensor API. Any protein
# whose signature is sufficiently close to any intracellular metabolite ligand
# signature has its regulatory occupancy modulated by that molecule. Bit 14 of
# the protein signature chooses potentiation versus inhibition. Default ligand
# signatures are deliberately >2 bits from ancestral proteins, so this adds an
# evolvable capability without silently rewriting the ancestral baseline.
var allostery_enabled: bool = true
var allosteric_max_distance: int = 2
var allosteric_distance_decay: float = 0.90
var allosteric_km_per_volume: float = 0.50
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
	assert(transcription_rate >= 0.0 and translation_rate >= 0.0)
	assert(mrna_decay_rate > 0.0 and protein_decay_rate > 0.0)
	assert(transcription_atp_cost_per_unit >= 0.0 and translation_atp_cost_per_unit >= 0.0)
	assert(nuc_cost_per_mrna_unit >= 0.0 and aa_cost_per_protein_unit >= 0.0)
	assert(expression_noise_fraction >= 0.0 and expression_noise_fraction <= 1.0)
	assert(expression_partition_jitter >= 0.0 and expression_partition_jitter < 0.45)
	assert(regulatory_max_distance >= 0 and regulatory_max_distance <= 16)
	assert(regulatory_distance_decay >= 0.0 and regulatory_gain >= 0.0)
	assert(regulatory_min_factor >= 0.0 and regulatory_max_factor >= regulatory_min_factor)
	assert(allosteric_max_distance >= 0 and allosteric_max_distance <= 16)
	assert(allosteric_distance_decay >= 0.0 and allosteric_km_per_volume > 0.0 and allosteric_gain >= 0.0)
	assert(allosteric_min_factor >= 0.0 and allosteric_max_factor >= allosteric_min_factor)
	assert(division_volume > ancestor_volume)
	assert(partition_jitter >= 0.0 and partition_jitter < 0.5)
	assert(promoter_mutation_rate_per_gene >= 0.0 and promoter_mutation_rate_per_gene <= 1.0)
	assert(signature_mutation_rate_per_gene >= 0.0 and signature_mutation_rate_per_gene <= 1.0)
	assert(regulatory_signature_mutation_rate_per_gene >= 0.0 and regulatory_signature_mutation_rate_per_gene <= 1.0)
	assert(neutral_marker_mutation_rate_per_gene >= 0.0 and neutral_marker_mutation_rate_per_gene <= 1.0)
	assert(promoter_mutation_step_max >= 1)
	assert(protein_signature_bits >= 1 and protein_signature_bits <= 16)
