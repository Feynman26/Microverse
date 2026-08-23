# Microverse Architecture

## 1. Purpose

Microverse is an individual-based artificial-life simulator for a deliberately small digital microbial world. Its purpose is not to reproduce one named bacterium quantitatively. It is to create a mechanistic substrate in which metabolism, regulation, heredity, mutation, ecological interaction and stochasticity can generate autonomous evolutionary outcomes.

The core design rule is:

> The engine defines possible physics and chemistry. Genomes determine how cells exploit those possibilities. Ecological and evolutionary outcomes are measured after they happen; they are not scripted as goals.

## 2. System boundaries

### Engine-defined and non-evolvable

- dimensionality and geometry of the chamber;
- numerical integration rules;
- diffusion law and boundary conditions;
- primitive molecular identities/properties;
- allowed low-level interaction grammar;
- conservation/invariant checks;
- heredity mechanics and mutation operators;
- random-number generator semantics;
- event and snapshot formats.

### Evolvable

Introduced progressively from M3 onward:

- promoter strengths and regulatory binding;
- protein abundance and stability;
- catalytic affinity/specificity;
- transporter affinity/capacity;
- repair investment and mutation rate;
- gene duplication/deletion;
- secretion and uptake capabilities;
- adhesion/motility traits;
- eventually mobile genetic elements.

### Explicitly not encoded as cell goals

The simulation must not expose `fitness`, `cooperate`, `cheat`, `predate`, `hibernate`, `seek_food`, `become_resistant` or equivalent goal primitives to the cell model. Such labels may exist only in analytics that classify an observed trajectory after the fact.

## 3. Runtime layers

### Layer A — simulation clock and deterministic stochasticity

`SimulationEngine` owns the tick counter, virtual time, deterministic RNG, ordering of simulation phases and event generation. Rendering never advances biological state directly.

For a fixed code version, configuration, initial snapshot, RNG state and intervention stream, the trajectory must be reproducible.

### Layer B — spatial world

The world is a 2D lattice of chemical concentration fields. M0-M2 use reflecting/no-flux boundaries. Later experiment modes may add chemostat-style inlets/outlets without changing cell semantics.

Each chemical is represented by an independent `ChemicalField`, making it possible to add metabolites without changing the diffusion solver.

### Layer C — cells

Cells are data objects, not Godot scene nodes. This is essential for performance and for headless execution. A cell stores position, lineage metadata and physiological state. In later milestones it also owns a genome, transcript state and proteome.

### Layer D — intracellular biology

M0-M2 use compressed transport, respiration, ATP, precursor, ROS/damage, repair, growth and division. M3-M7 progressively replace fixed ancestral scales with genome/proteome-derived kinetics.

### Layer E — ecology/evolution

There is no separate ecological equation such as Lotka-Volterra. Population dynamics arise from individual uptake, growth, death, division and spatial interaction. Lineages are analytical objects reconstructed from birth/division/mutation history.

### Layer F — observer/experimenter

The UI is an observer and perturbation console. It can change environmental boundary conditions, create snapshots, fork worlds, sample/freeze lineages and later perform controlled genetic manipulations. These interventions are recorded as events.

## 4. Tick semantics

The simulation uses explicit phases to avoid accidental code-order selection.

1. Apply scheduled experimental interventions/boundary conditions.
2. Advance extracellular diffusion.
3. Let all living cells sample the same post-diffusion environment and submit uptake requests.
4. Allocate scarce local resources proportionally between simultaneous requests.
5. Apply membrane uptake.
6. Update regulatory state (M5+).
7. Update transcription/translation (M5+).
8. Execute intracellular reactions.
9. Pay maintenance costs.
10. Update ROS, damage and repair.
11. Synthesize biomass/grow.
12. Update DNA replication/cell-cycle state (expanded in M3+).
13. Resolve physiological deaths and lysis.
14. Resolve physical overlap/mechanics (M6+).
15. Resolve division and inherited partitioning.
16. Apply mutations to daughter genomes (M3+).
17. Record events and analytics.
18. Validate invariants.

This phased model is intentionally different from naïvely calling `cell.step()` sequentially. Sequential direct consumption would make array position an unintended heritable advantage under scarcity.

## 5. Numerical strategy

### Diffusion

M0 uses an explicit 5-point finite-difference stencil with no-flux boundaries:

`C_next = C + alpha * (left + right + up + down - 4C)`

where `alpha = D * dt / dx^2`.

For the 2D stencil the configuration requires `alpha <= 0.25`. Violating this is an error, not a warning. This solver is intentionally simple and transparent for the initial model.

If later biological timescales require much larger `D/dx^2`, do not silently enlarge `dt`. Prefer one of:

- multiple diffusion substeps per biological tick;
- semi-implicit/implicit diffusion;
- quasi-steady diffusion solve.

The chosen transition must be covered by mass-conservation and convergence tests.

### Resource competition

Cells first request uptake from an unchanged environmental snapshot. For each lattice site and chemical, if total requested uptake exceeds available amount, each request receives the same proportional scaling factor.

This preserves symmetry between physiologically identical competitors sharing a site. Spatial position and physiology may cause different requests; iteration order may not.

### Intracellular model

Continuous pools are deterministic conditional on the stochastic events already sampled. M0-M2 do not attempt molecule-by-molecule stochastic simulation. Later expression noise will use explicitly seeded stochastic terms rather than global random calls.

## 6. State ownership

### Simulation state

Authoritative state must live only in simulation-domain objects:

- world chemical fields;
- cells;
- genome/proteome state;
- RNG state;
- clock;
- event sequence;
- scheduled interventions.

UI nodes must not own authoritative biological values.

### Derived state

Charts, colours, lineage labels, inferred ecotypes, novelty scores and explanations are derived analytics and can be recalculated.

## 7. Persistence

Persistence will have two complementary forms.

### Event log

Compact append-only semantic events:

- world created;
- ancestor seeded;
- birth/division/death;
- mutation;
- intervention;
- lineage freeze/reintroduction;
- detected ecological transition.

### Snapshots

Periodic complete state needed to resume efficiently:

- code/model schema version;
- configuration;
- tick and virtual time;
- RNG state;
- all chemical fields;
- all living cells and genomes;
- next IDs;
- pending interventions.

Replay starts from a snapshot and replays the subsequent event/intervention stream. Snapshot hashes will later support exact-fork verification.

## 8. Genome/protein architecture target

The first ancestral traits in `CellState` are scaffolding. They must not become permanent hard-coded phenotype sliders.

Target gene representation:

- stable gene ID for lineage analysis;
- promoter/basal expression parameters;
- zero or more regulatory binding sites;
- coding sequence/signature;
- protein localization;
- degradation/stability parameter;
- optional secretion/membrane localization tags produced by evolved sequence features.

Target protein representation:

- digital molecular signature;
- concentration/abundance;
- catalytic affinity to reaction motifs;
- binding affinity to regulator/promoter motifs;
- membrane/secreted/cytosolic localization;
- energetic/material synthesis cost;
- degradation rate.

Protein/reaction matching should support promiscuity: a protein can weakly catalyse several compatible reactions. This creates a route from mutation to side activity to duplication/divergence without a bespoke `new_function` mutation.

## 9. Chemistry grammar target

The world will use a finite primitive chemistry rather than unrestricted symbolic reactions. The initial target is roughly 15-30 metabolites spanning external resources, energetic currencies, carbon intermediates, biomass precursors, waste/stress compounds and potential signalling molecules.

Reaction definitions include:

- substrate stoichiometry;
- product stoichiometry;
- motif/signature used for catalytic matching;
- thermodynamic/directionality class;
- catalytic ceiling;
- compartment eligibility.

A cell does not automatically possess a reaction because it exists in the catalog. It requires sufficient catalytic activity from its current proteome.

## 10. Evolution target

Mutation operators are introduced in increasing structural risk:

1. scalar regulatory/kinetic mutation;
2. protein-signature point mutation;
3. regulatory-site mutation;
4. gene duplication;
5. gene deletion;
6. rearrangement/inversion;
7. mobile elements/horizontal transfer.

Every operator must preserve parseable genomes. Mutation probabilities are parameters of the replication machinery and later become partly evolvable through repair investment.

There is no mutation operator named after an adaptive result.

## 11. Spatial mechanics target

M6 replaces the current non-interacting cell positions with disk-like particles that grow, overlap transiently and undergo deterministic mechanical relaxation/pushing. The solver must prevent cells from occupying identical physical space while avoiding preferential directions or ID-order biases.

Adhesion and active motility are later protein-mediated properties. A biofilm, cluster or chemotactic trajectory is therefore an outcome, not a primitive state.

## 12. Analytics target

Analytics are observational and never feed survival unless an explicit experiment intervention uses them.

Planned detectors:

- persistent lineage expansion;
- stable polymorphism/coexistence;
- resource-use innovation;
- new secretion/uptake coupling;
- cross-feeding candidate;
- genome expansion/reduction;
- mutation-rate shift;
- population bottleneck/extinction;
- repeated environment-conditioned regulatory phenotype;
- candidate signalling relationship.

Detectors produce hypotheses, not truths. The experimenter can test causality by forking/reverting/competition assays.

## 13. Explainability strategy

`Explain this` will be a causal comparison pipeline, not free-form guessing. For a focal cell/lineage it should:

1. identify a comparison ancestor or sister lineage;
2. enumerate genomic differences;
3. replay or clone with selected mutations reverted when feasible;
4. quantify changes in protein abundance, fluxes, ATP, damage, division time and descendants;
5. report uncertainty when changes are epistatic or history-dependent.

## 14. Performance strategy

Correctness precedes optimization. Likely optimization order:

1. profile headless engine;
2. reduce allocation in hot loops;
3. use packed arrays/struct-of-arrays for large cell populations;
4. update diffusion at its own substep cadence;
5. cache reaction-affinity neighborhoods;
6. multithread only deterministic, commutative phases;
7. consider native GDExtension only after profiling proves GDScript inadequate.

The initial world remains intentionally small so architecture can be validated before performance complexity is introduced.

## 15. Non-negotiable invariants

- No materially negative chemical concentration.
- Closed diffusion conserves each extracellular chemical within numerical tolerance.
- A division conserves partitioned intracellular pools except explicit division cost/dissipation.
- Same snapshot + RNG state + intervention stream reproduces state checksum for the same model version.
- Cell IDs are never reused.
- Parentage is immutable.
- UI frame rate does not alter biological trajectory.
- No hidden fitness scalar participates in survival or reproduction.
