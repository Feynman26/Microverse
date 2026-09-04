# P3-A — Dense metabolic execution

## Purpose

P3-A is the first throughput-oriented Architecture v2 change. P0 measured
metabolism at approximately 43% of total tick time on the AMD Ryzen 5 4500U.
The M10 solver copied the complete metabolite dictionary and allocated
potential-flux, substrate-demand, scale, effective-flux and delta dictionaries
in every numerical substep: six times per cell and tick by default.

P3-A replaces only those temporary calculations. It does not yet replace the
externally compatible cell-pool dictionary or implement batching across cells.

## Dense representation

`MetaboliteCatalog` now declares one immutable canonical order and stable
integer index for all 19 metabolites. That order is exactly the alphabetical
order historically constructed at runtime.

`CompiledReactionNetwork` converts each ordered reaction list into:

- packed substrate indices and coefficients;
- packed product indices and coefficients;
- stable reaction IDs in the supplied order.

The compiled object is created once per `SimulationEngine`. Reversed or partial
reaction lists receive their own exact order when used by controlled assays.

## Reusable workspace

Each cell owns eight derived `PackedFloat64Array` buffers for pools, snapshots,
substrate demand/scales, reaction potential/cumulative fluxes, deltas and
catalytic activities. Buffers are reused rather than allocated per metabolic
substep, and only the cumulative-flux buffer requires a full clear per tick.

The workspace is computational cache only. It is excluded from checksums,
snapshots, inheritance, division and all biological decisions. Daughters and
restored cells begin with an empty workspace that is reconstructed on demand.

## Exact reference mode

`MetabolicSolver.step_legacy_reference()` preserves the closed M10 dictionary
algorithm as the shadow oracle. `SimConfig.metabolic_use_dense_solver` selects
between the two exact execution paths for paired tests and benchmarks. The
switch is deliberately absent from authoritative snapshot data because it is
an execution strategy, not scientific state.

The dense implementation preserves:

- reaction order;
- substrate/product key order;
- fair common-snapshot allocation;
- Float64 accumulation order per metabolite and reaction;
- substep count and duration;
- nonnegativity thresholds;
- returned reaction-flux dictionary order.

## Paired benchmark

Run the same initial state and seed through both solvers:

```bash
godot --headless --path . --script benchmarks/p3_dense_metabolism.gd -- \
  --warmup=10 \
  --ticks=30 \
  --populations=1,16,64,256,1000 \
  --output=p3-dense-metabolism.json
```

Every paired scenario alternates which solver runs first on each measured tick,
asserts exact checksum equality after every tick, and verifies final event-log
equality. This controls CPU boost and thermal-order bias. The report contains
elapsed time, ticks/s and legacy/dense speedup for each population.

## Acceptance gate

P3-A is eligible to merge only when:

- dense and frozen M10 solvers match exactly across repeated substeps, partial
  networks and reversed reaction order;
- complete simulations match checksum and event history at every tested tick;
- all M0–M10, P0–P2 and snapshot/fork regressions pass;
- the paired benchmark completes without a scientific divergence;
- three sequential target-laptop runs show at least 10% lower median tick time
  at 64, 256 and 1,000 cells;
- no 1- or 16-cell median regression exceeds 5%.

If the target is not met, the dense path remains unmerged while profiling is
used to identify buffer or boundary overhead. No biological tolerance will be
relaxed to obtain a speedup.
