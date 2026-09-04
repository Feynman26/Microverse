# P2 — Runtime contracts and structural telemetry

## Purpose

P2 creates the migration seam that lets Architecture v2 replace M10 internals
incrementally. It is not a throughput optimization. Its acceptance criterion is
that callers stop owning or mutating `SimulationEngine` directly while the
legacy adapter produces the exact same scientific trajectory.

P1 already established responsive worker ownership. P2 preserves that thread
boundary and places `ISimulationBackend` inside it. Future P3 stores can
therefore replace the adapter without changing the interactive clock, UI or
experiment orchestration.

## Versioned contracts

All contracts begin at schema version 1:

| Contract | Responsibility |
|---|---|
| `ISimulationBackend` | authoritative command execution, visual snapshots, state summary and telemetry |
| `SimulationCommand` | immutable-by-contract scientific input |
| `SimulationDelta` | compact receipt for one committed command |
| `VisualSnapshot` | presentation-only state with no authoritative references |
| structured telemetry | phase timers plus deterministic work counters |

Supported P2 scientific commands are:

- seed one ancestor at an explicit or default position;
- advance an exact number of complete biological ticks;
- apply an external environment schedule at an explicit tick;
- execute a deterministic serial-transfer intervention.

Interactive pause, requested clock and shutdown remain runtime-control
commands. They do not enter the scientific backend because they cannot change
the result of a fixed number of completed ticks.

## Delta semantics

A delta records command sequence, command kind, before/after tick and simulated
time, before/after population, before/after event count and small
command-specific details. It deliberately does not deep-copy cells, fields or
the event log. This keeps P2 from adding a large allocation tax to every tick.

Command and delta sequences are operational metadata. They are excluded from
scientific checksums, RNG state and snapshots used for exact replay.

## Legacy adapter

`LegacySimulationBackend` is the only P2 implementation. It owns one unchanged
M10 `SimulationEngine` and translates the versioned contracts into existing
operations. The interactive worker and `ExperimentRunner` mutate scientific
state only through this backend.

M8 analytics temporarily use `legacy_inspection_state()`. This is an explicitly
read-only migration seam, not permission to mutate the engine. P3 replaces it
with indexed/dense observation views as stores move behind the interface.

## Structural telemetry

The P0 wall-clock profiler retains its backward-compatible phase report and now
also emits a versioned envelope containing:

- phase timer distributions;
- completed ticks;
- cell visits by intracellular, lifecycle, secretion and mechanics phases;
- basal cell/resource pairs;
- secondary cell/metabolite pairs;
- eligible field/lattice slots;
- extracellular reaction/lattice slots.

Counters are collected only when performance profiling is enabled. They remain
outside biology and allow P3 to report cost per unit of scientific work rather
than only total milliseconds.

## Scientific and performance boundary

P2 does not change phase order, equations, RNG consumption, time step,
diffusion, transport, expression, metabolism, mutation, division or mechanics.
It also does not claim faster ticks. The P0 benchmark remains the speed oracle;
P2 merely adds structural counters to the same report.

The backend adds compact command/delta dictionaries. P3 must remove or batch
those allocations if profiling shows they are material. No native extension or
approximation is authorized in P2.

## Exit gate

P2 is complete only when:

- interactive execution owns a backend rather than `SimulationEngine`;
- `ExperimentRunner` routes seeds, environments, ticks and interventions
  through versioned commands;
- command, delta, visual snapshot and telemetry schemas are validated;
- direct M10 execution and adapter execution have exact checksum and event-log
  equivalence;
- serial transfer remains deterministic;
- telemetry on/off has exact scientific equivalence;
- P1 responsiveness and the complete M0–M10 regression suite still pass.
