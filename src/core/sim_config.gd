extends RefCounted
class_name SimConfig

# Time and world scale. Units are intentionally explicit even though M0-M2
# remain a compressed biochemical model rather than a quantitative bacterium.
var tick_dt_min: float = 0.10
var world_width: int = 64
var world_height: int = 64
var grid_cell_size_um: float = 1.0
var max_cells: int = 64
var seed: int = 824718

# Environment concentrations (arbitrary concentration units for the first
# vertical slice; the model is dimensionless-but-consistent at this stage).
var initial_glucose: float = 4.0
var initial_oxygen: float = 5.0
var glucose_diffusion: float = 0.80
var oxygen_diffusion: float = 1.00

# Membrane transport.
var glucose_transport_vmax: float = 0.60
var glucose_transport_km: float = 0.80
var oxygen_transport_vmax: float = 0.90
var oxygen_transport_km: float = 0.60
var intracellular_pool_capacity_per_volume: float = 8.0

# Compressed respiration: glucose + 2 O2 -> ATP + precursor + ROS.
var respiration_vmax: float = 0.50
var oxygen_per_glucose: float = 2.0
var atp_yield_per_glucose: float = 6.0
var precursor_yield_per_glucose: float = 1.0
var ros_yield_per_glucose: float = 0.04

# Basal physiology.
var maintenance_atp_rate_per_volume: float = 0.08
var growth_vmax: float = 0.20
var growth_atp_per_precursor: float = 1.5
var volume_yield_per_precursor: float = 0.50
var basal_ros_decay_rate: float = 0.08
var ros_damage_rate: float = 0.04
var basal_repair_rate: float = 0.02
var repair_atp_cost: float = 0.20
var lethal_damage: float = 3.0
var lethal_energy_debt: float = 4.0

# Cell cycle.
var ancestor_volume: float = 1.0
var division_volume: float = 2.0
var division_atp_cost: float = 1.0
var partition_jitter: float = 0.02
var daughter_offset_grid: float = 0.20

func validate() -> void:
	assert(tick_dt_min > 0.0)
	assert(world_width > 2 and world_height > 2)
	assert(grid_cell_size_um > 0.0)
	assert(max_cells >= 1)
	assert(initial_glucose >= 0.0 and initial_oxygen >= 0.0)
	assert(glucose_diffusion >= 0.0 and oxygen_diffusion >= 0.0)
	var alpha_g := glucose_diffusion * tick_dt_min / (grid_cell_size_um * grid_cell_size_um)
	var alpha_o := oxygen_diffusion * tick_dt_min / (grid_cell_size_um * grid_cell_size_um)
	assert(alpha_g <= 0.25, "Explicit 2D diffusion unstable for glucose: alpha must be <= 0.25")
	assert(alpha_o <= 0.25, "Explicit 2D diffusion unstable for oxygen: alpha must be <= 0.25")
	assert(division_volume > ancestor_volume)
	assert(partition_jitter >= 0.0 and partition_jitter < 0.5)
