# P1 — Responsive interactive runtime

## Purpose

P1 prevents an expensive biological tick from blocking Godot input and
rendering. It changes interactive scheduling and presentation only. The M10
kernel, tick duration, phase order, RNG consumption, scientific state and
checksums remain unchanged.

P0 showed why this boundary is required. On the available target laptop, the
interactive run plateaued near 7.5 ticks/s around 64 cells whether 1× or 10×
was requested, while requesting 100× made the UI non-responsive. The former UI
called `simulation.step()` directly from `_process()`, so its nominal 12 ms
budget could be checked only after a non-preemptible tick had already blocked
the frame.

## Ownership and thread boundary

`InteractiveSimulationRuntime` owns one dedicated `Thread`. Only that worker
creates and accesses the authoritative `SimulationEngine` instance. The scene
tree never reads cells, fields, genomes, RNG state or event storage directly.

The worker periodically constructs a `VisualSnapshot` containing only fresh
scalar values, dictionaries and packed arrays. A snapshot is immutable by
contract after publication. The main thread briefly acquires a mutex to obtain
the latest snapshot reference, then renders without holding a lock.

The UI sends pause and requested-speed changes through a separate control
mutex. Commands become effective between complete biological ticks. No command
can interrupt a tick or reorder a scientific phase.

## Clock and backpressure

The biological tick remains exactly `0.10` simulated minutes. With the default
clock definition:

| Requested clock | Requested ticks/s |
|---:|---:|
| 1× | 10 |
| 10× | 100 |
| 100× | 1,000 |
| 1,000× | 10,000 |
| 5,000× | 50,000 |

`InteractiveClock` converts elapsed wall time into pending simulated minutes.
Pending demand is bounded to a half-second horizon at the requested clock. If
hardware cannot keep up, P1 does not skip biological ticks and does not let an
unbounded queue grow. It exposes:

- requested and achieved ticks/s;
- requested and achieved clock multiplier;
- pending backlog in simulated minutes and ticks;
- cumulative unserved wall-clock demand;
- overload events and a visible overload state.

Unserved demand is computational telemetry, not elapsed biological time. Only
completed exact ticks advance `simulation_time_min`.

Pausing clears pending wall-clock demand. Resuming begins at the current wall
time, so time spent paused is never replayed as simulation debt.

## Visual sampling

Snapshot production is wall-clock sampled independently of biological ticks:

- below 100×: at most every 100 ms;
- from 100× through 999×: at most every 250 ms;
- at 1,000× and above: at most every 500 ms.

Every biological tick is still executed. Lower visual cadence under load
changes only observation frequency. The current snapshot copies the glucose
field and compact cell position/volume arrays; full presentation contracts and
delta streams remain P2 work.

## Safety behavior

The interactive 256-cell guard remains explicitly computational. Reaching it
pauses the worker at a completed tick boundary and cannot be resumed silently.
It is never reported as ecological carrying capacity or stationary phase.

The project now starts at 1×. Requesting a faster clock is allowed even when it
is unattainable; the UI remains responsive and reports the achieved rate and
backlog honestly.

Stopping the scene requests worker shutdown and joins the thread. Because M10
ticks are not yet internally preemptible, shutdown may wait for the current
tick to finish, but render/input are not blocked during ordinary operation.

## Scientific equivalence gate

The P1 test starts the worker paused, requests an exact number of tick-boundary
steps and compares its result with the ordinary single-threaded M10 kernel. The
checksum, event count and population must be exact. It also verifies bounded
backpressure, visible unserved demand, immutable prior snapshots and
non-blocking speed commands.

Run the focused gate with:

```bash
godot --headless --path . --script tests/p1_responsive_runtime_tests.gd
```

P1 does not claim increased biological throughput. P3–P8 address the measured
hot paths and multiscale execution. P1 guarantees that insufficient throughput
is handled as an observable runtime condition instead of a frozen UI.

## Exit criteria

- the scene-tree thread never calls `SimulationEngine.step()`;
- the worker is the sole owner of authoritative simulation state;
- prior visual snapshots remain unchanged after later ticks;
- exact worker stepping matches the M10 checksum and event/population result;
- pause and speed commands are accepted without waiting for a biological tick;
- overload has bounded backlog and visible unserved demand;
- the 256-cell safety guard pauses at a complete tick boundary;
- the full M0–M10 regression suite and P0 baseline smoke gate continue to pass;
- a target-workstation run at 100× remains interactive even when achieved speed
  is far below the request.
