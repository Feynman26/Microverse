# Microverse

Microverse is a small, reproducible artificial-life laboratory focused on autonomous evolution in a spatially explicit digital microbial world.

The project deliberately separates **digital physics/chemistry** from **evolvable biology**. The engine defines space, diffusion, conservation rules, possible molecular interactions, energetic costs, heredity and mutation mechanisms. It does **not** encode goals such as fitness, cooperation, competition, dormancy, predation or specialization. Those are phenomena to detect after they emerge from cell physiology, inheritance, ecology and stochasticity.

## Design principles

1. **No explicit fitness score.** Reproductive success is the evolutionary consequence.
2. **Deterministic stochasticity.** A world seed plus the same initial state and interventions must reproduce the same history.
3. **Simulation/render separation.** The simulation advances in ticks independently of the UI so it can run headless and much faster than real time.
4. **Conservation first.** Chemistry and transport must not create negative concentrations or accidental mass/energy.
5. **Mechanisms, not behaviours.** Genes, proteins, reactions, membrane transport, damage, repair and division are modelled; named strategies are not hard-coded.
6. **Observability.** Important births, deaths, mutations, lineage changes and interventions must be reconstructable.
7. **Small world, deep biology.** Initial scope is a microchamber with tens of cells and a deliberately compressed biochemical model.

## Current vertical slice: M0-M5A

The implemented causal chain is now:

`environment -> uptake -> DNA -> stochastic transcription -> mRNA -> stochastic translation -> protein cohorts -> catalytic affinity -> reaction flux -> metabolite pools -> BIO -> division -> inheritance -> mutation`

Current components include:

- closed 2D microchamber;
- diffusing carbon, oxygen, nitrogen and phosphorus resources;
- explicit numerical-stability validation;
- simultaneous/fair local resource allocation among competing cells;
- explicit 12-locus heritable genome with stable locus IDs;
- promoter-code, 16-bit protein-signature and neutral-marker mutation;
- deterministic mutation genealogy and genotype fingerprints;
- 19 intracellular digital metabolites;
- 12 structurally balanced catalytic reactions;
- Hamming-distance protein/reaction affinity landscape;
- catalytic promiscuity and side activity;
- order-independent simultaneous metabolic solver;
- explicit ATP/ADP and NAD/NADH currency pools;
- C/N/P-conserving biomass precursor synthesis;
- catalytic structural biomass (`BIO`) assembly;
- physical cell volume derived from `BIO` rather than a direct growth scalar;
- oxidative and oxygen-independent metabolic routes;
- explicit ROS production, detoxification and downstream damage;
- per-cell mRNA and protein molecular state;
- stochastic transcription, translation, mRNA decay and protein decay;
- deterministic Poisson molecular noise from the universe RNG;
- explicit ATP, nucleotide and amino-acid expression costs;
- order-independent synthesis allocation under ATP/NUC/AA scarcity;
- protein sequence cohorts that preserve inherited molecular identity across coding mutations;
- stochastic but conserving mRNA/protein partition at division;
- catalysis driven by realized protein abundance rather than promoter code;
- dormant waste-recovery reactions that can become active through ordinary coding mutation;
- deterministic whole-state checksum including expression state;
- observer UI with population/genetics/resource, focal-cell chemistry, expression and flux state;
- headless invariant, landscape, metabolic and expression suites.

### Why M4 mattered

M4 made coding mutations functionally meaningful. Protein-signature mutations can change catalytic affinity and therefore reaction flux without a fitness function.

For example, ancestral W1 recovery is completely inactive across the ancestral proteome. A single appropriate bit mutation moves one protein's reaction distance from 5 to 4, producing measurable W1-recovery flux and downstream C2. Its evolutionary consequence depends on the actual environment and descendants.

### Why M5-A matters

M5-A removes the shortcut in which promoter code directly represented protein abundance. A promoter now affects transcription; transcripts generate proteins; physical protein molecules catalyse reactions.

Importantly, molecular identity persists through mutation. If daughter DNA changes from coding signature A to B, inherited A transcripts/proteins remain A until they decay. Newly expressed molecules are B. The cell may therefore temporarily contain both ancestral and mutant protein cohorts.

That creates genotype-to-phenotype lag, stochastic clone differences and sister asymmetry as ordinary molecular consequences rather than named behavior variables.

### Current abstraction boundary before M5-B

M5-A has explicit expression but not an evolvable regulatory network. Promoter strength remains a basal heritable scalar. M5-B adds promoter binding motifs, generic regulator occupancy, activator/repressor effects and a minimal chemical/receptor bridge before M6 spatial cell mechanics.

Basal membrane transport and repair also remain transitional mechanisms rather than fully genome-derived systems.

## Requirements

- Godot 4.7.2 (standard build; no .NET required)

CI installs the same stable engine version.

## Run

Open `project.godot` in Godot and run the project.

The interactive M10 simulation runs on a dedicated worker. Rendering and input
consume immutable visual snapshots, so requesting a clock faster than the
computer can achieve does not require the scene-tree thread to execute or wait
for a biological tick. The status panel reports requested/achieved speed,
bounded backlog and unserved clock demand. The project starts at 1×.

Controls:

- `Space`: pause/resume;
- `1`: 1x experimental clock;
- `2`: 10x;
- `3`: 100x;
- `4`: 1000x;
- `5`: 5000x.

The renderer shows the extracellular glucose field, cells, virtual time, population, generation, genotype/mutation counts, environmental resources and a focal cell's BIO, ATP/ADP, carbon/waste pools, ROS/damage, total mRNA/protein, protein-cohort count, recent transcription/translation cost and dominant reaction flux.

## Headless validation

```bash
godot --headless --editor --path . --quit
godot --headless --path . --script tests/run_tests.gd
godot --headless --path . --script tests/m4_landscape_tests.gd
godot --headless --path . --script tests/m4_metabolism_tests.gd
godot --headless --path . --script tests/m5_expression_tests.gd
```

## P0 performance baseline

Architecture v2 begins with observational P0 instrumentation. Run the frozen
1/16/64/256/1000-cell headless suite with:

```bash
GODOT_BIN=godot bash benchmarks/run_p0.sh p0-baseline.json
```

See `docs/P0_PERFORMANCE_BASELINE.md` for the measurement contract, target
workstation procedure and interpretation limits. Profiling is disabled in
ordinary simulations and does not enter biological state or checksums.

P1 uses the P0 evidence to separate interactive scheduling from presentation
without changing M10 biology. See `docs/P1_RESPONSIVE_RUNTIME.md` for worker
ownership, immutable snapshot, clock/backpressure and exact-equivalence rules.

The current gates include diffusion/nonnegativity; fair resource allocation; exact genome and molecular inheritance; reaction-order independence; digital C/N/P conservation; ATP+ADP and NAD+NADH conservation; sparse-but-connected catalytic-landscape statistics; aerobic versus hypoxic metabolism; one-mutation activation of a dormant metabolic route; catalytic BIO assembly; exact stochastic expression replay; clone phenotypic divergence across seeds; expression structural/energy accounting; gene-order independence; mutation-to-protein lag; stochastic sister partition; Poisson statistics; resource-scarcity scaling; resource-supported division; extinction without replenishable energy; and complete same-seed state/history reproducibility.

## Documentation

- `docs/ARCHITECTURE.md` — system boundaries, tick semantics, target genome/protein/chemistry architecture, persistence, analytics and performance rules.
- `docs/ARCHITECTURE_V2.md` — scalable multiscale target, Biological LOD,
  deterministic parallelism and P0–P10 migration gates.
- `docs/P1_RESPONSIVE_RUNTIME.md` — interactive worker ownership, immutable
  visual snapshots, honest clock/backpressure and equivalence gate.
- `docs/SCIENTIFIC_MODEL.md` — biological abstractions, equations, assumptions, limitations and validation hierarchy.
- `docs/M3_GENETICS.md` — gene/genome representation, inheritance and mutation semantics.
- `docs/M4_CATALYTIC_LANDSCAPE.md` — catalytic-affinity landscape design and statistical gate.
- `docs/M4_METABOLIC_INTEGRATION.md` — intracellular pools, flux solver, conservation semantics and BIO growth.
- `docs/M5_EXPRESSION_CORE.md` — stochastic mRNA/protein dynamics, molecular sequence cohorts, expression material accounting and genotype-to-phenotype lag.
- `docs/ROADMAP.md` — gate-driven implementation through regulation, spatial ecology, open evolution, laboratory tooling, evolved mutation rate, sensing/motility and horizontal transfer.
