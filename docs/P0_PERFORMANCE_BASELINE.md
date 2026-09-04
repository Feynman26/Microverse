# P0 — Performance baseline and instrumentation

## Purpose

P0 measures the M10 engine before Architecture v2 changes execution semantics. It is observational: it does not optimize, reorder phases, change stochastic streams or alter biological state.

## What is measured

The engine records wall-clock microseconds for:

- diffusion;
- basal transport;
- secondary transport;
- protein secretion;
- extracellular catalysis;
- total intracellular biology, plus per-cell samples for expression, proteome
  constraint, metabolism, DNA replication, maintenance, damage/repair and
  validation;
- death/lysis;
- division/mutation;
- mechanics;
- invariant checks;
- complete tick.

The benchmark also records throughput, Godot/OS/CPU manifest, static memory, peak static memory, object counters, final population, final checksum, events and mechanics diagnostics.

Godot does not expose a reliable cross-platform per-phase allocation counter in GDScript. P0 therefore reports memory and object deltas and uses external Godot profiler captures for allocation call stacks when an interactive diagnosis is required. The report must not label memory deltas as exact allocation counts.

## Frozen scenarios

All standard scenarios use the M10 model, a 64×64 world, deterministic non-overlapping grid placement, disabled mutation and populations:

`1, 16, 64, 256, 1000`

The population is instantiated at the beginning rather than grown. This isolates steady per-tick scaling from the unrelated time required to evolve/grow each population. P0 does not interpret these fixtures as biological experiments.

## Run the complete baseline

From the repository root with Godot 4.7.2:

```bash
GODOT_BIN=godot bash benchmarks/run_p0.sh p0-baseline.json
```

Or run a custom diagnostic:

```bash
godot --headless --path . --script benchmarks/p0_baseline.gd -- \
  --warmup=5 --ticks=20 --populations=16,256 --output=p0-custom.json
```

Warmup samples are discarded. The output contains p50/p95/max phase time over measured ticks and is printed between `P0_RESULT_JSON_BEGIN` and `P0_RESULT_JSON_END` even when no file is requested.

Whole-engine phase samples occur once per tick. Names prefixed with `cell_`
occur once per living cell per tick; their totals reveal aggregate subsystem
cost while their percentiles characterize individual-cell variability.

## CI smoke gate

CI runs a deliberately short `1,16,64` smoke benchmark. It validates that the runner remains executable and emits every required scenario. It does **not** compare absolute CI timing because shared runners are noisy and differ from the target laptop.

## Target workstation procedure

1. Plug in AC power and select the same performance profile for every run.
2. Close unrelated heavy applications.
3. Record the exact commit and Godot 4.7.2 version.
4. Run the complete suite three times.
5. Keep all raw JSON reports; do not average away an anomalous run silently.
6. Use medians across runs to rank hotspots and p95/max to detect stalls.
7. Capture interactive Godot profiler data separately at requested clocks 1×, 10× and 100×.
8. Record render frame p50/p95, achieved ticks/s and backlog from the UI session.

The headless benchmark measures simulation throughput. Render frame time must be captured in an interactive run and must not be inferred from headless results.

## Interpretation rule

P0 may identify bottlenecks but cannot approve a solution. P1 scope is selected only after the target-workstation evidence exists. Any later optimization must compare against the same frozen scenarios and pass the complete M0–M10 regression suite.
