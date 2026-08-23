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

## Current vertical slice: M0-M4

The implemented causal chain is now:

`environment -> uptake -> genome -> catalytic affinity -> reaction flux -> metabolite pools -> ATP/precursors/BIO -> division -> inheritance -> mutation`

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
- proteome-abundance opportunity cost;
- dormant waste-recovery reactions that can become active through ordinary coding mutation;
- division into two new daughter identities with stochastic molecular partitioning;
- deterministic seeded RNG and whole-state checksum;
- observer UI with population/genetics/resource and focal-cell chemistry/flux state;
- headless invariant, landscape and metabolic-integration suites.

### Why M4 matters

The genome is no longer merely inherited information. Protein-signature mutations can change catalytic affinity and therefore reaction flux.

For example, ancestral W1 recovery is completely inactive across the ancestral proteome. A single appropriate bit mutation in one protein signature moves its reaction distance from 5 to 4, producing measurable W1-recovery flux and downstream C2. The engine does not label that mutation beneficial; its evolutionary consequence depends on resource availability, competing pathways, expression cost and descendants.

This is the first complete Microverse path from molecular mutation to potentially selectable physiology without a fitness function.

### Current abstraction boundary before M5

M4 still uses promoter code as a constitutive abundance proxy. M5 replaces that bridge with explicit mRNA/protein abundance, transcription, translation, degradation, regulation and seeded expression noise. Basal membrane transport and repair are also still transitional mechanisms rather than fully genome-derived systems.

## Requirements

- Godot 4.7.2 (standard build; no .NET required)

CI installs the same stable engine version.

## Run

Open `project.godot` in Godot and run the project.

Controls:

- `Space`: pause/resume;
- `1`: 1x experimental clock;
- `2`: 10x;
- `3`: 100x;
- `4`: 1000x;
- `5`: 5000x.

The renderer shows the extracellular glucose field, cells, virtual time, population, generation, genotype count, mutation count, environmental resources and a focal cell's BIO, ATP/ADP, central carbon pools, waste pools, ROS/damage and dominant reaction flux.

## Headless validation

```bash
godot --headless --editor --path . --quit
godot --headless --path . --script tests/run_tests.gd
godot --headless --path . --script tests/m4_landscape_tests.gd
godot --headless --path . --script tests/m4_metabolism_tests.gd
```

The M4 gates include reaction-order independence; nonnegative chemistry; digital C/N/P conservation; ATP+ADP and NAD+NADH conservation; sparse-but-connected catalytic-landscape statistics; zero-expression/zero-flux controls; aerobic versus hypoxic pathway behavior; one-mutation activation of a dormant metabolic route; catalytic BIO assembly; explicit expression-energy cost; resource-supported division; starvation extinction; and exact same-seed history/checksum reproducibility.

## Documentation

- `docs/ARCHITECTURE.md` — system boundaries, tick semantics, target genome/protein/chemistry architecture, persistence, analytics and performance rules.
- `docs/SCIENTIFIC_MODEL.md` — biological abstractions, equations, assumptions, limitations and validation hierarchy.
- `docs/M3_GENETICS.md` — gene/genome representation, inheritance and mutation semantics.
- `docs/M4_CATALYTIC_LANDSCAPE.md` — catalytic-affinity landscape design and statistical gate.
- `docs/M4_METABOLIC_INTEGRATION.md` — intracellular pools, flux solver, conservation semantics, BIO growth and current M4/M5 boundary.
- `docs/ROADMAP.md` — gate-driven implementation through regulation, spatial ecology, open evolution, laboratory tooling, evolved mutation rate, sensing/motility and horizontal transfer.
