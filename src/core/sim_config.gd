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

# M7 secondary extracellular chemistry. These are generic chamber properties,
# keyed by intracellular metabolite identity so the world/cell mapping remains
# centralized in MetaboliteCatalog. Zero initial amounts keep the M0-M6 basal
# environment unchanged while allowing lysis/secretion to populate the fields.
const SECONDARY_EXTRACELLULAR_IDS: Array[String] = [
	"C2", "C3", "W1", "W2", "CO2", "AA", "LIP", "NUC", "ROS", "X"
]
var secondary_extracellular_initial: Dictionary = {
	"C2": 0.0,
	"C3": 0.0,
	"W1": 0.0,
	"W2": 0.0,
	"CO2": 0.0,
	"AA": 0.0,
	"LIP": 0.0,
	"NUC": 0.0,
	"ROS": 0.0,
	"X": 0.0
}
var secondary_extracellular_diffusion: Dictionary = {
	"C2": 0.65,
	"C3": 0.55,
	"W1": 0.60,
	"W2": 0.50,
	"CO2": 1.20,
	"AA": 0.35,
	"LIP": 0.10,
	"NUC": 0.25,
	"ROS": 1.00,
	"X": 0.75
}

var glucose_transport_vmax: float = 0.60
var glucose_transport_km: float = 0.80
var oxygen_transport_vmax: float = 0.90
var oxygen_transport_km: float = 0.60
var nitrogen_transport_vmax: float = 0.45
var nitrogen_transport_km: float = 0.50
var phosphorus_transport_vmax: float = 0.30
var phosphorus_transport_km: float = 0.40
var intracellular_pool_capacity_per_volume: float = 8.0

# M7-B generic secondary membrane exchange. Activity comes from realized
# proteins in the transport recognition landscape. The gradient determines the
# direction, while every moved unit pays ATP and the transporter already pays
# the finite-expression/proteome opportunity cost inherited from M5.
var secondary_transport_vmax_per_reference_protein: float = 0.40
var secondary_transport_gradient_km: float = 0.50
var secondary_transport_atp_cost_per_unit: float = 0.02

# M7-E generic protein secretion and extracellular catalysis. A protein must
# carry the sequence-level secretion motif before any of its realized molecules
# can leave the cell. Secretion removes those exact protein molecules and spends
# ATP. Extracellular protein diffuses slowly and catalyses only reactions for
# which its sequence has ordinary Hamming-distance affinity.
var extracellular_protein_diffusion: float = 0.08
var extracellular_protein_secretion_fraction_per_min: float = 0.20
var extracellular_protein_secretion_atp_cost_per_unit: float = 0.01
var extracellular_catalysis_rate_scale: float = 1.00
var extracellular_catalysis_km: float = 0.25

var metabolic_substeps_per_tick: int = 6
var metabolic_km_per_volume: float = 0.20
var metabolic_rate_scale: float = 0.85
var biomass_units_per_volume: float = 1.0
var initial_atp_per_volume: float = 2.0
var initial_adp_per_volume: float = 8.0
var initial_nad_per_volume: float = 6.0
var initial_nadh_per_volume: float = 1.0

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
var proteome_capacity_reference_units: float = 5.0
var translation_capacity_fraction_of_proteome_per_min: float = 0.06

var regulation_enabled: bool = true
var regulatory_max_distance: int = 3
var regulatory_distance_decay: float = 0.80
var regulatory_gain: float = 0.60
var regulatory_min_factor: float = 0.25
var regulatory_max_factor: float = 1.75

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

# M10 mechanistic DNA replication. Genome copying is a time-resolved molecular
# process rather than an instantaneous division surcharge. Every copied gene
# consumes nucleotide material and ATP. A fixed genes/min copying throughput
# makes larger genomes take proportionally more time even under abundant
# resources. Sequence-derived repair activity can lower replication errors, but
# increases ATP cost per copied gene; there is no mutator/fidelity phenotype flag.
var evolvable_replication_enabled: bool = true
var genome_replication_initiation_volume: float = 1.20
var genome_replication_gene_copy_rate_per_min: float = 1.00
var genome_replication_nuc_cost_per_gene: float = 0.010
var genome_replication_atp_cost_per_gene: float = 0.020
var dna_repair_max_distance: int = 4
var dna_repair_distance_decay: float = 0.70
var dna_repair_fidelity_gain: float = 8.0
var dna_repair_atp_cost_per_gene_activity: float = 0.010
var baseline_point_error_rate_per_gene: float = 0.002
var minimum_point_error_rate_per_gene: float = 0.00001
var baseline_structural_error_rate_per_genome: float = 0.002
var minimum_structural_error_rate_per_genome: float = 0.00001

# M6 physical-cell mechanics and broad-phase indexing.
var ancestor_radius_grid: float = 0.45
var mechanical_relaxation_iterations: int = 16
var mechanical_relaxation_fraction: float = 0.80
var mechanical_overlap_tolerance: float = 1e-4
var mechanical_use_spatial_index: bool = true
var mechanical_neighbor_bucket_size_grid: float = 2.0

# Legacy M3 point-mutation knobs remain available for direct historical assays.
# The production simulation switches to M10 replication-derived probabilities
# when evolvable_replication_enabled is true.
var mutation_enabled: bool = true
var promoter_mutation_rate_per_gene: float = 0.001
var signature_mutation_rate_per_gene: float = 0.001
var regulatory_signature_mutation_rate_per_gene: float = 0.001
var neutral_marker_mutation_rate_per_gene: float = 0.001
var promoter_mutation_step_max: int = 250
var protein_signature_bits: int = 16

func extracellular_initial_amount(metabolite_id: String) -> float:
	match metabolite_id:
		"G": return initial_glucose
		"O2": return initial_oxygen
		"NH4": return initial_nitrogen
		"P": return initial_phosphorus
		_:
			assert(secondary_extracellular_initial.has(metabolite_id), "Missing extracellular initial amount: %s" % metabolite_id)
			return float(secondary_extracellular_initial[metabolite_id])

func extracellular_diffusion_coefficient(metabolite_id: String) -> float:
	match metabolite_id:
		"G": return glucose_diffusion
		"O2": return oxygen_diffusion
		"NH4": return nitrogen_diffusion
		"P": return phosphorus_diffusion
		_:
			assert(secondary_extracellular_diffusion.has(metabolite_id), "Missing extracellular diffusion coefficient: %s" % metabolite_id)
			return float(secondary_extracellular_diffusion[metabolite_id])

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
	assert(secondary_extracellular_initial.size() == SECONDARY_EXTRACELLULAR_IDS.size())
	assert(secondary_extracellular_diffusion.size() == SECONDARY_EXTRACELLULAR_IDS.size())
	for metabolite_id in SECONDARY_EXTRACELLULAR_IDS:
		assert(secondary_extracellular_initial.has(metabolite_id), "Missing secondary extracellular initial amount: %s" % metabolite_id)
		assert(secondary_extracellular_diffusion.has(metabolite_id), "Missing secondary extracellular diffusion coefficient: %s" % metabolite_id)
		var initial_amount: float = float(secondary_extracellular_initial[metabolite_id])
		var diffusion: float = float(secondary_extracellular_diffusion[metabolite_id])
		assert(initial_amount >= 0.0)
		assert(diffusion >= 0.0)
		assert(diffusion * tick_dt_min / dx2 <= 0.25, "Unstable secondary extracellular diffusion: %s" % metabolite_id)
	assert(glucose_transport_vmax >= 0.0 and oxygen_transport_vmax >= 0.0)
	assert(nitrogen_transport_vmax >= 0.0 and phosphorus_transport_vmax >= 0.0)
	assert(secondary_transport_vmax_per_reference_protein >= 0.0)
	assert(secondary_transport_gradient_km > 0.0)
	assert(secondary_transport_atp_cost_per_unit >= 0.0)
	assert(extracellular_protein_diffusion >= 0.0)
	assert(extracellular_protein_diffusion * tick_dt_min / dx2 <= 0.25, "Unstable extracellular protein diffusion")
	assert(extracellular_protein_secretion_fraction_per_min >= 0.0)
	assert(extracellular_protein_secretion_atp_cost_per_unit >= 0.0)
	assert(extracellular_catalysis_rate_scale >= 0.0)
	assert(extracellular_catalysis_km > 0.0)
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
	assert(genome_replication_initiation_volume >= ancestor_volume and genome_replication_initiation_volume < division_volume)
	assert(genome_replication_gene_copy_rate_per_min > 0.0)
	assert(genome_replication_nuc_cost_per_gene >= 0.0)
	assert(genome_replication_atp_cost_per_gene >= 0.0)
	assert(dna_repair_max_distance >= 0 and dna_repair_max_distance <= 16)
	assert(dna_repair_distance_decay >= 0.0)
	assert(dna_repair_fidelity_gain >= 0.0)
	assert(dna_repair_atp_cost_per_gene_activity >= 0.0)
	assert(baseline_point_error_rate_per_gene >= 0.0 and baseline_point_error_rate_per_gene <= 1.0)
	assert(minimum_point_error_rate_per_gene >= 0.0 and minimum_point_error_rate_per_gene <= baseline_point_error_rate_per_gene)
	assert(baseline_structural_error_rate_per_genome >= 0.0 and baseline_structural_error_rate_per_genome <= 1.0)
	assert(minimum_structural_error_rate_per_genome >= 0.0 and minimum_structural_error_rate_per_genome <= baseline_structural_error_rate_per_genome)
	assert(ancestor_radius_grid > 0.0)
	assert(mechanical_relaxation_iterations >= 1)
	assert(mechanical_relaxation_fraction > 0.0 and mechanical_relaxation_fraction <= 1.0)
	assert(mechanical_overlap_tolerance >= 0.0)
	assert(mechanical_neighbor_bucket_size_grid > 0.0)
	assert(2.0 * ancestor_radius_grid < float(mini(world_width - 1, world_height - 1)))
	assert(promoter_mutation_rate_per_gene >= 0.0 and promoter_mutation_rate_per_gene <= 1.0)
	assert(signature_mutation_rate_per_gene >= 0.0 and signature_mutation_rate_per_gene <= 1.0)
	assert(regulatory_signature_mutation_rate_per_gene >= 0.0 and regulatory_signature_mutation_rate_per_gene <= 1.0)
	assert(neutral_marker_mutation_rate_per_gene >= 0.0 and neutral_marker_mutation_rate_per_gene <= 1.0)
	assert(promoter_mutation_step_max >= 1)
	assert(protein_signature_bits >= 1 and protein_signature_bits <= 16)