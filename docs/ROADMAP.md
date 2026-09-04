# Microverse Roadmap

This roadmap is gate-driven. A milestone is complete only when its exit criteria pass; later biological complexity must not be added on top of unresolved numerical or causal defects.

## Architecture v2 prerequisite — P0 performance baseline

Before M11 adds new biology, P0 instruments the closed M10 engine and freezes
standard scaling scenarios. P0 is observational only: it may measure phase
costs, throughput, memory/object counters and frame behavior, but it may not
change scientific execution semantics. `docs/P0_PERFORMANCE_BASELINE.md`
defines the runner and evidence required to select P1 work.

## M0 — Numerical world foundation

### Deliverables

- deterministic RNG wrapper;
- explicit simulation clock independent of rendering;
- 2D scalar chemical field;
- no-flux boundaries;
- configurable diffusion coefficients;
- mass/nonnegativity invariants;
- headless test runner;
- minimal observer UI.

### Exit criteria

- closed diffusion conserves mass within tolerance;
- concentrations never become materially negative;
- same seed/configuration produces identical checksum;
- UI speed and frame rate do not change simulation state for a fixed number of ticks;
- unstable diffusion configuration fails validation.

## M1 — Basal single-cell physiology

### Deliverables

- extracellular glucose and oxygen;
- saturable membrane uptake;
- fair simultaneous resource allocation;
- compressed respiration;
- ATP and precursor pools;
- maintenance cost and energy debt;
- ROS/damage/repair;
- biomass growth;
- physiological death and partial lysis.

### Exit criteria

- cell survives/grows only when a usable energetic route exists;
- starving cells eventually die;
- identical cells in an identical site receive equal resource under scarcity;
- increasing maintenance cost monotonically increases energetic burden in controlled tests;
- no intracellular pool becomes negative.

## M2 — Division, lineage and observation

### Deliverables

- division threshold and explicit cost;
- two new daughter IDs per division;
- seeded asymmetric partitioning;
- immutable parentage/generation;
- birth/division/death event log;
- minimal time controls and field/cell renderer;
- state checksum.

### Exit criteria

- resource-fed ancestor can generate descendants;
- daughter pools sum to parent post-division pools within tolerance;
- genealogy is reconstructable;
- repeated same-seed run gives identical division history.

## M3 — Heritable genome and mutation

### Goal

Replace fixed cell trait scales with inherited molecular information.

### Deliverables

- `Genome`, `Gene`, mutation-event data models;
- stable gene IDs distinct from sequence/signature;
- genome replication at division;
- seeded point mutation of promoters/protein signatures;
- mutation event log;
- genotype hash;
- clone/freeze/reintroduce primitives;
- neutral mutation mode for validation.

### Exit criteria

- mutation-free daughters are genotype-identical;
- mutation frequency matches configured distribution statistically;
- same seed reproduces exact mutations;
- neutral mutations do not alter physiology;
- no direct fitness variable exists.

## M4 — Evolvable proteins and metabolism

### Goal

Move phenotype from scalar traits to molecular mechanism.

### Deliverables

- finite metabolite/reaction catalog;
- protein digital signatures;
- reaction motifs/signatures;
- affinity function with tunable landscape ruggedness;
- catalytic flux from protein abundance/affinity;
- alternative metabolic route including oxygen-independent low-yield path;
- biomass precursor chemistry;
- explicit synthesis costs.

### Research gate before acceptance

Characterize the digital catalytic landscape by random-walk experiments:

- distribution of affinities for random proteins;
- probability a point mutation materially changes primary activity;
- frequency of useful side activities;
- correlation between neighboring genotypes;
- probability duplicated copies can diverge without immediate lethality.

### Exit criteria

- metabolism functions without hard-coded `respiration_scale`/`growth_scale` phenotype controls;
- at least two viable metabolic strategies exist under different environments;
- catalytic promiscuity is neither nearly zero everywhere nor universal;
- energy/material accounting remains bounded and nonnegative.

## M5 — Gene expression and regulatory evolution

### Deliverables

- mRNA and protein abundance;
- transcription/translation degradation;
- energetic/material expression costs;
- regulatory proteins and binding motifs;
- activator/repressor occupancy;
- promoter mutation;
- expression noise using deterministic RNG;
- environment-sensitive regulation through molecular sensors rather than named responses.

### Exit criteria

- an inducible response can be built from generic mechanisms without special-case code;
- regulatory mutations alter timing/abundance through the network;
- clones exhibit seeded phenotypic heterogeneity when noise is enabled;
- expression costs create measurable trade-offs;
- regulation can evolve differently in stable versus fluctuating test environments.

## M6 — Physical cells and spatial ecology

### Deliverables

- disk-like cell geometry;
- radius/volume relation;
- overlap detection;
- deterministic mechanical relaxation/pushing;
- division placement without persistent overlap;
- surface boundary interaction;
- optional adhesion primitive controlled by molecular phenotype;
- spatial indexing for performance.

### Exit criteria

- cells cannot occupy the same physical volume after relaxation;
- mechanics do not systematically favor low/high IDs or a coordinate direction;
- growing colonies generate nutrient gradients;
- clonal sectors can emerge from local reproduction;
- mechanics remains deterministic for a fixed seed and state.

## M7 — Ecological innovation substrate

### Goal

Create enough low-level possibility for diversification without coding interaction labels.

### Deliverables

- expanded 15-30 metabolite chemistry;
- diffusible waste/intermediates;
- secretion/export transport;
- uptake of secondary metabolites;
- extracellular catalytic reactions;
- lysis release of reusable compounds;
- toxin-like chemistry as ordinary damaging compounds, not a combat API;
- potential signal compound that has no built-in semantic role;
- gene duplication/deletion.

### Exit criteria

Controlled precursor experiments demonstrate the engine can support, without dedicated behavior functions:

- competition by depletion;
- cross-feeding;
- extracellular public-good production;
- exploitation by a nonproducer;
- detoxification benefit;
- antagonism through a secreted damaging compound.

These are capability tests, not evidence they will necessarily evolve spontaneously.

## M8 — Open evolutionary experiments

### Deliverables

- long-run headless experiment runner;
- batch seeds/replicates;
- environment schedules;
- bottlenecks/serial transfer;
- resource gradients;
- fluctuating oxygen/nutrient conditions;
- lineage frequency analytics;
- genotype/phenotype diversity metrics;
- resource-use profile metrics;
- automatic extinction and persistent-coexistence detection.

### Core experiment families

1. stable abundant resource;
2. stable scarce resource;
3. oxygen fluctuation;
4. feast/famine;
5. spatial gradient;
6. periodic bottleneck;
7. waste-rich secondary niche;
8. public-good opportunity;
9. environmental shock/rescue;
10. replicate histories from adjacent seeds.

### Exit criteria

- replicate runs can diverge stochastically while remaining individually reproducible;
- lineage selection changes with environment;
- at least one environment supports persistent diversification in controlled evolutionary trials;
- extinctions are allowed and recoverable only through explicit snapshot/reseeding.

## M9 — Experimental laboratory tooling

### Deliverables

- versioned snapshot format;
- exact fork with RNG state;
- freezer/strain library;
- clone selected cell;
- mutation revert/knockout experiment;
- lineage competition assay;
- ancestor-versus-evolved assay;
- timeline/replay;
- side-by-side fork comparison;
- `Investigate` workflow for detected transitions.

### Exit criteria

- a pre-event world can be forked into two bit-identical states;
- one controlled intervention can be applied to only one branch;
- branch differences are attributable to intervention plus subsequent stochastic path;
- saved experiments record model version, configuration, seed, state hash and interventions.

## M10 — Evolved mutation rate and genome architecture

### Deliverables

- explicit DNA-repair investment;
- mutation rate coupled to repair/fidelity machinery;
- cost of fidelity;
- gene rearrangement/inversion;
- genome-size replication cost;
- deletion bias configurable through mechanism rather than outcome;
- lineage-level mutation-rate analytics.

### Exit criteria

- high fidelity is costly;
- low fidelity increases both beneficial opportunity and deleterious load;
- mutator lineages can transiently succeed in a changing environment but are not universally optimal;
- genome reduction/expansion has mechanistic costs.

## M11 — Motility, sensing and information

### Deliverables

- generic membrane receptors;
- intracellular signal proteins;
- active motor protein consuming ATP;
- run/tumble-like motion constructed from generic motor switching, not `chemotaxis()`;
- regulatory coupling from sensed molecules to motor state;
- evolving receptor/signature affinity.

### Exit criteria

- directed migration can emerge from a generic sensor-regulator-motor circuit;
- motor activity has energetic cost;
- movement can be maladaptive in homogeneous environments;
- a formerly metabolic/waste compound can in principle acquire informational value if sensing it predicts future local conditions.

## M12 — Horizontal transfer and selfish genetic elements

This is deliberately late because it adds another replication level.

### Deliverables

- mobile genome segment representation;
- transfer/contact mechanism;
- replication/transmission cost;
- compatibility/insertion rules;
- mobile-element lineage tracking;
- optional parasitic element that can persist without host-beneficial function.

### Exit criteria

- host and mobile-element evolutionary interests can diverge;
- transfer is spatially/local-mechanism dependent;
- identical initial state/seed reproduces transfer history;
- no hard-coded `beneficial_plasmid` category exists.

## Cross-cutting workstreams

### A. Reproducibility

Every milestone extends snapshot/checksum coverage. RNG calls must occur only through the simulation RNG service.

### B. Scientific auditability

Every abstraction, equation and parameter is documented with status:

- arbitrary/compressed;
- biologically motivated;
- calibrated;
- experimentally validated within the simulation.

### C. Performance

Benchmark headless ticks/s at standard world sizes after each milestone. Optimization cannot alter deterministic semantics without an explicit architecture decision.

### D. UX

UI should progressively expose:

- field heatmaps;
- cell physiology;
- genome/proteome;
- regulatory network;
- metabolic fluxes;
- lineage tree;
- experiment timeline;
- fork comparison.

The UI must never become the authoritative state store.

## Definition of an interesting emergent result

A result is promoted from visual curiosity to a Microverse finding only if:

1. it was not directly encoded as a named goal/behavior;
2. it persists long enough to distinguish it from transient noise;
3. the originating run is exactly reproducible from saved state/seed;
4. lower-level causal variables can explain a plausible mechanism;
5. a fork, revert, competition assay or environmental perturbation can test that mechanism;
6. alternative seeds/replicates establish whether it is common, contingent or rare.
