extends RefCounted
class_name InteractiveSimulationRuntime

const InteractiveClockScript = preload("res://src/runtime/interactive_clock.gd")
const LegacySimulationBackendScript = preload("res://src/runtime/legacy_simulation_backend.gd")
const SimulationCommandScript = preload("res://src/runtime/simulation_command.gd")

const SPEED_MEASURE_USEC: int = 500000
const OVERLOAD_VISIBLE_USEC: int = 2000000
const IDLE_DELAY_USEC: int = 1000
const SNAPSHOT_NORMAL_USEC: int = 100000
const SNAPSHOT_ACCELERATED_USEC: int = 250000
const SNAPSHOT_EXTREME_USEC: int = 500000

var _thread: Thread = null
var _control_mutex := Mutex.new()
var _snapshot_mutex := Mutex.new()
var _config = null
var _started: bool = false

var _stop_requested: bool = false
var _paused_requested: bool = false
var _requested_multiplier: float = 1.0
var _pending_exact_ticks: int = 0
var _compute_limit_reached: bool = false
var _command_revision: int = 0

var _latest_snapshot: Dictionary = {}
var _snapshot_sequence: int = 0

func start(config, start_paused: bool = false, initial_multiplier: float = 1.0) -> Error:
	assert(not _started, "Interactive runtime can only be started once")
	assert(config != null)
	assert(initial_multiplier > 0.0)
	config.validate()
	_config = config
	_paused_requested = start_paused
	_requested_multiplier = initial_multiplier
	_thread = Thread.new()
	var error: Error = _thread.start(Callable(self, "_thread_main"))
	if error == OK:
		_started = true
	return error

func stop() -> void:
	if not _started:
		return
	_control_mutex.lock()
	_stop_requested = true
	_command_revision += 1
	_control_mutex.unlock()
	_thread.wait_to_finish()
	_started = false

func set_paused(value: bool) -> bool:
	_control_mutex.lock()
	if _compute_limit_reached and not value:
		_control_mutex.unlock()
		return false
	if _paused_requested != value:
		_paused_requested = value
		if not value:
			_pending_exact_ticks = 0
		_command_revision += 1
	_control_mutex.unlock()
	return true

func toggle_paused() -> bool:
	_control_mutex.lock()
	if _compute_limit_reached:
		_control_mutex.unlock()
		return false
	_paused_requested = not _paused_requested
	if not _paused_requested:
		_pending_exact_ticks = 0
	_command_revision += 1
	_control_mutex.unlock()
	return true

func set_requested_multiplier(value: float) -> void:
	assert(value > 0.0)
	_control_mutex.lock()
	if not is_equal_approx(_requested_multiplier, value):
		_requested_multiplier = value
		_command_revision += 1
	_control_mutex.unlock()

# Exact stepping is useful for deterministic debug and verifies that threaded
# scheduling does not alter the scientific trajectory. It is accepted only
# while paused and never bypasses the computational guard.
func request_exact_ticks(count: int = 1) -> bool:
	assert(count >= 0)
	if count == 0:
		return true
	_control_mutex.lock()
	if not _paused_requested or _compute_limit_reached:
		_control_mutex.unlock()
		return false
	_pending_exact_ticks += count
	_command_revision += 1
	_control_mutex.unlock()
	return true

func latest_snapshot() -> Dictionary:
	_snapshot_mutex.lock()
	var result: Dictionary = _latest_snapshot
	_snapshot_mutex.unlock()
	return result

func _thread_main() -> void:
	# The worker owns the backend exactly as P1 owned SimulationEngine. P2 moves
	# the dependency behind ISimulationBackend without changing thread ownership.
	var backend = LegacySimulationBackendScript.new(_config)
	backend.execute(SimulationCommandScript.seed_ancestor())
	var clock = InteractiveClockScript.new(float(_config.tick_dt_min))
	var control: Dictionary = _control_state()
	clock.set_requested_multiplier(float(control["requested_multiplier"]))

	var now_usec: int = Time.get_ticks_usec()
	var last_loop_usec: int = now_usec
	var last_snapshot_usec: int = now_usec
	var speed_window_started_usec: int = now_usec
	var speed_window_ticks: int = 0
	var actual_ticks_per_second: float = 0.0
	var last_overload_usec: int = -OVERLOAD_VISIBLE_USEC
	var seen_overload_events: int = 0
	var seen_command_revision: int = int(control["revision"])
	var was_paused: bool = bool(control["paused"])
	_publish_snapshot(backend, clock, actual_ticks_per_second, now_usec, last_overload_usec, true)

	while true:
		now_usec = Time.get_ticks_usec()
		control = _control_state()
		if bool(control["stop"]):
			break

		var delta_sec: float = maxf(0.0, float(now_usec - last_loop_usec) / 1000000.0)
		last_loop_usec = now_usec
		# Attribute elapsed wall time to the control state that was active while
		# that interval elapsed. New commands take effect only at this tick
		# boundary, never retroactively across a long biological step.
		if not was_paused:
			clock.add_wall_time(delta_sec)
		clock.set_requested_multiplier(float(control["requested_multiplier"]))
		var paused: bool = bool(control["paused"])
		if paused != was_paused:
			# Pausing discards queued wall-clock demand; resuming begins at the
			# current wall time. Neither operation changes a biological tick.
			clock.clear_backlog()
		was_paused = paused

		if clock.overload_events != seen_overload_events:
			seen_overload_events = clock.overload_events
			last_overload_usec = now_usec

		var exact_tick: bool = _claim_exact_tick() if paused else false
		var automatic_tick: bool = not paused and clock.can_consume_tick()
		var stepped: bool = false
		var reached_compute_limit: bool = false
		if exact_tick or automatic_tick:
			var delta: Dictionary = backend.execute(SimulationCommandScript.advance_ticks(1))
			if automatic_tick:
				clock.consume_tick()
			speed_window_ticks += 1
			stepped = true
			if int(delta["population_after"]) >= int(_config.max_cells):
				_mark_compute_limit_reached()
				clock.clear_backlog()
				last_overload_usec = Time.get_ticks_usec()
				reached_compute_limit = true

		var completed_usec: int = Time.get_ticks_usec()
		var speed_elapsed_usec: int = completed_usec - speed_window_started_usec
		if speed_elapsed_usec >= SPEED_MEASURE_USEC:
			actual_ticks_per_second = (
				float(speed_window_ticks) * 1000000.0 / maxf(1.0, float(speed_elapsed_usec))
			)
			speed_window_started_usec = completed_usec
			speed_window_ticks = 0

		var command_changed: bool = int(control["revision"]) != seen_command_revision
		if command_changed:
			seen_command_revision = int(control["revision"])
		var snapshot_due: bool = completed_usec - last_snapshot_usec >= _snapshot_interval_usec(clock.requested_multiplier)
		if snapshot_due or command_changed or exact_tick or reached_compute_limit:
			_publish_snapshot(
				backend,
				clock,
				actual_ticks_per_second,
				completed_usec,
				last_overload_usec,
				true
			)
			last_snapshot_usec = completed_usec

		if not stepped:
			OS.delay_usec(IDLE_DELAY_USEC)

	_publish_snapshot(
		backend,
		clock,
		0.0,
		Time.get_ticks_usec(),
		last_overload_usec,
		false
	)

func _control_state() -> Dictionary:
	_control_mutex.lock()
	var result: Dictionary = {
		"stop": _stop_requested,
		"paused": _paused_requested,
		"requested_multiplier": _requested_multiplier,
		"compute_limit": _compute_limit_reached,
		"revision": _command_revision
	}
	_control_mutex.unlock()
	return result

func _claim_exact_tick() -> bool:
	_control_mutex.lock()
	var claimed: bool = _paused_requested and not _compute_limit_reached and _pending_exact_ticks > 0
	if claimed:
		_pending_exact_ticks -= 1
	_control_mutex.unlock()
	return claimed

func _mark_compute_limit_reached() -> void:
	_control_mutex.lock()
	_compute_limit_reached = true
	_paused_requested = true
	_pending_exact_ticks = 0
	_command_revision += 1
	_control_mutex.unlock()

func _publish_snapshot(
	backend,
	clock,
	actual_ticks_per_second: float,
	now_usec: int,
	last_overload_usec: int,
	worker_active: bool
) -> void:
	var snapshot: Dictionary = backend.capture_visual_snapshot()
	var control: Dictionary = _control_state()
	_snapshot_sequence += 1
	snapshot["runtime"] = {
		"sequence": _snapshot_sequence,
		"worker_active": worker_active,
		"paused": bool(control["paused"]),
		"compute_limit_reached": bool(control["compute_limit"]),
		"requested_multiplier": float(control["requested_multiplier"]),
		"requested_ticks_per_second": float(clock.requested_ticks_per_second()),
		"actual_ticks_per_second": actual_ticks_per_second,
		"achieved_multiplier": float(clock.achieved_multiplier(actual_ticks_per_second)),
		"backlog_min": float(clock.backlog_min),
		"backlog_ticks": float(clock.backlog_ticks()),
		"maximum_backlog_min": float(clock.maximum_backlog_min()),
		"unserved_requested_min": float(clock.unserved_requested_min),
		"overload_events": int(clock.overload_events),
		"overloaded": (
			now_usec - last_overload_usec <= OVERLOAD_VISIBLE_USEC
			or clock.backlog_ticks() >= 1.0
		)
	}
	_snapshot_mutex.lock()
	_latest_snapshot = snapshot
	_snapshot_mutex.unlock()

func _snapshot_interval_usec(multiplier: float) -> int:
	if multiplier >= 1000.0:
		return SNAPSHOT_EXTREME_USEC
	if multiplier >= 100.0:
		return SNAPSHOT_ACCELERATED_USEC
	return SNAPSHOT_NORMAL_USEC
