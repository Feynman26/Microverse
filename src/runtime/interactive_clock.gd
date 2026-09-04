extends RefCounted
class_name InteractiveClock

# Wall-clock demand controller for interactive execution. It schedules exact
# biological ticks but never changes their dt or skips scientific work. When
# hardware cannot satisfy the requested clock, pending demand is bounded and
# the discarded wall-clock demand remains observable.

const BASE_SIM_MIN_PER_REAL_SEC: float = 1.0
const BACKLOG_HORIZON_SEC: float = 0.50

var tick_dt_min: float
var requested_multiplier: float = 1.0
var backlog_min: float = 0.0
var unserved_requested_min: float = 0.0
var overload_events: int = 0

func _init(p_tick_dt_min: float = 0.10) -> void:
	tick_dt_min = p_tick_dt_min
	assert(tick_dt_min > 0.0)

func set_requested_multiplier(value: float) -> void:
	assert(value > 0.0)
	requested_multiplier = value
	_bound_backlog()

func add_wall_time(delta_sec: float) -> void:
	assert(delta_sec >= 0.0)
	backlog_min += delta_sec * BASE_SIM_MIN_PER_REAL_SEC * requested_multiplier
	_bound_backlog()

func can_consume_tick() -> bool:
	return backlog_min + 1e-12 >= tick_dt_min

func consume_tick() -> void:
	assert(can_consume_tick())
	backlog_min = maxf(0.0, backlog_min - tick_dt_min)

func clear_backlog() -> void:
	backlog_min = 0.0

func backlog_ticks() -> float:
	return backlog_min / tick_dt_min

func requested_ticks_per_second() -> float:
	return BASE_SIM_MIN_PER_REAL_SEC * requested_multiplier / tick_dt_min

func achieved_multiplier(ticks_per_second: float) -> float:
	return ticks_per_second * tick_dt_min / BASE_SIM_MIN_PER_REAL_SEC

func maximum_backlog_min() -> float:
	return maxf(tick_dt_min, BASE_SIM_MIN_PER_REAL_SEC * requested_multiplier * BACKLOG_HORIZON_SEC)

func _bound_backlog() -> void:
	var maximum: float = maximum_backlog_min()
	if backlog_min <= maximum:
		return
	unserved_requested_min += backlog_min - maximum
	backlog_min = maximum
	overload_events += 1
