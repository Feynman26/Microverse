# Microverse Architecture v2

## Decision

Microverse remains a bottom-up artificial-life laboratory. Architecture v2
scales the causal M0–M10 model without introducing fitness scores, scripted
behaviours, predefined species or human-like agents.

The target causal arc is:

`chemistry → cells → multicellularity → sensing/control → learning → communication → cognition → societies`

Complexity must be paid for in modeled material, energy, time and reproductive
opportunity. A future generative model is optional and late; it cannot replace
the evolutionary path that produces an organism, sensors, memory and actuators.

## Invariants retained from M0–M10

- one authoritative scientific state;
- no explicit fitness score;
- deterministic seeded stochasticity;
- simultaneous/fair resource competition;
- nonnegativity and declared material conservation;
- persistent molecular identity across mutation;
- physical expression, maintenance, repair and DNA-copying costs;
- reconstructable lineage, mutation and intervention history;
- exact snapshots, forks and replay;
- explicit computational limits that are never reported as biology.

M0–M10 is the compatibility oracle. An optimization is accepted only as exact
equivalence or as a named approximation with measured error.

## Runtime boundaries

### Simulation kernel

Owns time, phase order, RNG streams, authoritative stores, event scheduling and
state commits. It contains no UI state.

### World domain

Owns geometry, chunked chemical fields, boundaries, environmental sources and
the spatial index.

### Biology domain

Owns identity/lifecycle, geometry, metabolite pools, expression, genomes,
damage/repair, reproduction and later evolved capabilities.

The hot path moves from nested Dictionaries to dense Structure-of-Arrays stores
and integer-indexed catalogs. Immutable identical genomes are interned and
mutations use copy-on-write.

### Systems

Systems read a stable phase snapshot, produce proposals/deltas and commit them
at deterministic barriers:

- field transport;
- membrane exchange;
- extracellular catalysis;
- expression/regulation;
- metabolism;
- maintenance/damage;
- DNA replication;
- division/inheritance;
- mechanics;
- future sensing, control, learning and communication.

### Observation

Events, metrics, traces, causal explanations and performance measurements are
derived. Observation never consumes biological RNG or modifies state.

### Presentation

Godot consumes immutable visual snapshots. It interpolates positions, renders
fields as textures/batches and sends interventions as commands. Rendering never
advances or owns biology.

## Multiscale time

Architecture v2 replaces “every system for every cell on every tick” with:

- periodic updates for fast processes;
- adaptive integration under explicit error limits;
- dirty-driven regulation and transport;
- discrete threshold events for division, death and mutation;
- local interaction updates when neighborhoods change;
- on-demand diagnostics.

Time uses an integer monotonic representation. Equal-time events have stable
ordering by phase, entity, system and local sequence. Dormant/stable entities
use analytic catch-up, but integration must stop at any internal threshold that
could change the causal path.

## Biological level of detail

| Level | Representation | Intended use |
|---|---|---|
| L0 | exact molecular pools and cohorts | focal cells, rare mutants, validation |
| L1 | exact genetics, aggregated physiology | ordinary active cells |
| L2 | rates, reserves, damage, cycle and capabilities | stable populations |
| L3 | cohort count plus distributions/moments | abundant equivalent clones |
| L4 | density/biomass/genetic-frequency fields | future very large worlds |

Resolution follows scientific relevance, not camera distance. Mutation,
heterogeneous interaction, approaching division/death, multimodality or direct
inspection promotes detail. Stable equivalent entities may be demoted.

Every transition preserves population weight, modeled C/N/P, carrier semantics,
biomass, DNA/RNA/protein inventory, cycle state, genetic frequencies, genealogy
and RNG position. Every run reports time/population per LOD, transitions,
residuals and estimated error.

## Determinism and parallelism

Random streams derive from world seed, entity/cohort, system and event counter.
Thread scheduling, rendering, diagnostics and LOD changes cannot perturb
unrelated biological randomness.

Jobs write private buffers or exclusive partitions. Global reductions and
commits follow stable order. Supported modes are interactive, accelerated
interactive, headless throughput, exact debug and batch experiments.

The interactive engine reports requested versus achieved speed. It reduces
visual sampling or effective simulated speed instead of blocking the render
thread to catch up.

## Scalable physics and biology

- Chemical fields use contiguous double buffers and active chunks/halos.
- Uniform inactive chunks can sleep without dropping nonzero chemistry.
- Transport preserves simultaneous proposals and fair local allocation.
- Mechanics maintains a persistent spatial index, contact islands and sleeping.
- Metabolism uses compiled stoichiometry, cached catalytic activities and
  adaptive/batched integration.
- Regulation uses a compiled dependency graph and dirty propagation.
- Molecular decay and other stable dynamics use analytical catch-up.
- Replication stops at resource/threshold events and creates mutations only
  through the copying mechanism.
- Death/lysis applies conservative release packages.

Native C++ GDExtension or compute shaders are considered only after profiling
identifies a stable hot kernel and an exact reference implementation exists.

## Migration gates

1. **P0 — baseline:** observational phase timings, frozen population scenarios,
   memory/object counters, headless throughput and target-PC render capture.
2. **P1 — responsiveness:** simulation/render separation, honest backlog and
   visual sampling without scientific changes.
3. **P2 — contracts:** backend interface, commands/deltas, phase counters and a
   legacy adapter.
4. **P3 — dense data:** indexed catalogs, stores, interned genomes, compiled
   caches and reused buffers.
5. **P4 — incremental mechanics:** persistent index, contact islands and sleep.
6. **P5 — chunked fields:** contiguous buffers, active tiles and stable
   adaptive diffusion.
7. **P6 — exact multiscale scheduler:** integer clock, event buckets, dirty sets,
   threshold events and deterministic streams.
8. **P7 — adaptive/batched intracellular systems:** compiled metabolism and
   regulation with exact and characterized approximate modes.
9. **P8 — Biological LOD:** L0–L3 cohorts, conservative transitions,
   genealogy compression and error ledger.
10. **P9 — v2 persistence/laboratory:** incremental snapshots, M10 migration,
    exact forks and causal comparison.
11. **P10 — platform closure:** complete M0–M10 equivalence, performance
    evidence, documented limits and release readiness.

New biology begins only after P10. Every gate is a separate, reviewable change
with evidence and a rollback path.

## Global completion criteria

Architecture v2 is ready when:

- interactive use stays responsive;
- interactive, headless and batch modes share one kernel;
- M0–M10 passes under v2;
- snapshots/forks/replay remain valid;
- performance regression checks exist;
- cost scales with active regions/entities;
- L0–L3 transitions are conservative and auditable;
- computational limits and approximation error remain visible;
- a focal organism can be inspected molecularly without charging that cost to
  every organism.

P0 is therefore the only authorized first implementation step. Its measurements
select P1; P0 itself must not optimize or change causal semantics.
