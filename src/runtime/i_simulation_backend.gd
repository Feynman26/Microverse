extends RefCounted
class_name ISimulationBackend

# P2 contract implemented by every authoritative simulation backend. GDScript
# has no nominal interface type, so this base class supplies explicit failing
# stubs and one version number that callers can validate.

const INTERFACE_VERSION: int = 1

func execute(_command: Dictionary) -> Dictionary:
	assert(false, "Simulation backend must implement execute(command)")
	return {}

func capture_visual_snapshot() -> Dictionary:
	assert(false, "Simulation backend must implement capture_visual_snapshot()")
	return {}

func state_summary() -> Dictionary:
	assert(false, "Simulation backend must implement state_summary()")
	return {}

func telemetry_report() -> Dictionary:
	assert(false, "Simulation backend must implement telemetry_report()")
	return {}

func metadata() -> Dictionary:
	assert(false, "Simulation backend must implement metadata()")
	return {}

# Temporary migration seam. P2 routes all mutations through commands, while
# existing M8 analytics may inspect the M10 object through this adapter. The
# returned object is read-only by contract and will be replaced by dense views
# during P3.
func legacy_inspection_state():
	assert(false, "Simulation backend does not expose a legacy inspection state")
	return null
