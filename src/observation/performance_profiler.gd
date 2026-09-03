extends RefCounted
class_name PerformanceProfiler

# P0 observational profiler. It owns wall-clock measurements only and is kept
# outside checksums, snapshots and all biological decisions.

var enabled: bool = false
var _phase_samples_usec: Dictionary = {}

func _init(p_enabled: bool = false) -> void:
	enabled = p_enabled

func begin_phase() -> int:
	return Time.get_ticks_usec() if enabled else 0

func end_phase(phase: StringName, started_usec: int) -> void:
	if not enabled:
		return
	var elapsed: int = maxi(0, Time.get_ticks_usec() - started_usec)
	if not _phase_samples_usec.has(phase):
		_phase_samples_usec[phase] = PackedInt64Array()
	var samples: PackedInt64Array = _phase_samples_usec[phase]
	samples.append(elapsed)
	_phase_samples_usec[phase] = samples

func reset() -> void:
	_phase_samples_usec.clear()

func report() -> Dictionary:
	var result: Dictionary = {}
	var phases: Array = _phase_samples_usec.keys()
	phases.sort()
	for phase_variant in phases:
		var phase: StringName = phase_variant
		result[String(phase)] = _summarize(_phase_samples_usec[phase])
	return result

func _summarize(source: PackedInt64Array) -> Dictionary:
	if source.is_empty():
		return {
			"samples": 0,
			"total_usec": 0,
			"mean_usec": 0.0,
			"min_usec": 0,
			"p50_usec": 0,
			"p95_usec": 0,
			"max_usec": 0
		}
	var ordered: Array[int] = []
	var total: int = 0
	for sample in source:
		var value: int = int(sample)
		ordered.append(value)
		total += value
	ordered.sort()
	return {
		"samples": ordered.size(),
		"total_usec": total,
		"mean_usec": float(total) / float(ordered.size()),
		"min_usec": ordered[0],
		"p50_usec": _percentile(ordered, 0.50),
		"p95_usec": _percentile(ordered, 0.95),
		"max_usec": ordered[ordered.size() - 1]
	}

func _percentile(ordered: Array[int], fraction: float) -> int:
	var index: int = ceili(fraction * float(ordered.size())) - 1
	return ordered[clampi(index, 0, ordered.size() - 1)]
