extends RefCounted
class_name PerformanceProfiler

# P0 observational profiler. It owns wall-clock measurements only and is kept
# outside checksums, snapshots and all biological decisions.

const STRUCTURED_REPORT_SCHEMA_VERSION: int = 1

var enabled: bool = false
var _phase_samples_usec: Dictionary = {}
var _work_counters: Dictionary = {}

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

func add_work(counter: StringName, amount: int = 1) -> void:
	if not enabled:
		return
	assert(amount >= 0)
	_work_counters[counter] = int(_work_counters.get(counter, 0)) + amount

func reset() -> void:
	_phase_samples_usec.clear()
	_work_counters.clear()

func report() -> Dictionary:
	var result: Dictionary = {}
	var phases: Array = _phase_samples_usec.keys()
	phases.sort()
	for phase_variant in phases:
		var phase: StringName = phase_variant
		result[String(phase)] = _summarize(_phase_samples_usec[phase])
	return result

func work_report() -> Dictionary:
	var result: Dictionary = {}
	var counters: Array = _work_counters.keys()
	counters.sort()
	for counter_variant in counters:
		result[String(counter_variant)] = int(_work_counters[counter_variant])
	return result

func structured_report() -> Dictionary:
	return {
		"schema_version": STRUCTURED_REPORT_SCHEMA_VERSION,
		"phase_timers": report(),
		"work_counters": work_report()
	}

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
