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

## Current vertical slice: M0-M3

The implemented causal foundation is:

`environmental nutrient -> membrane uptake -> metabolism -> ATP -> maintenance/repair/biomass -> growth -> division -> inherited genome -> mutation`

Current components:

- closed 2D microchamber;
- glucose and oxygen diffusion;
- explicit numerical stability validation;
- simultaneous/fair local resource allocation;
- saturable nutrient transport;
- compressed respiration and ATP generation;
- maintenance, energetic debt, ROS, damage and repair;
- precursor-dependent biomass growth;
- division into two new daughter identities with stochastic partitioning;
- explicit 12-locus heritable genome;
- deep-copy inheritance with no sibling/ancestor aliasing;
- promoter-code, protein-signature and neutral-marker mutations;
- deterministic mutation event history with genotype fingerprints;
- deterministic seeded RNG and whole-state checksum;
- genotype-diversity/mutation counters in the observer UI;
- headless invariant/statistical tests.

### Important M3 boundary

The M3 genome is **heritable but intentionally not yet adaptive**. Promoter and protein-signature mutations are molecular information only. Basal physiology still uses transitional fixed scales from M0-M2.

M4 gives protein signatures functional meaning through a validated catalytic-affinity/reaction landscape. M5 gives promoters and regulatory motifs functional meaning through costly gene expression and regulation. This staged design lets us prove inheritance and neutral evolution before introducing selection on molecular function.

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

The renderer currently shows the extracellular glucose field, cells, virtual time, population, generation, live genotype count, mutation-event count, resource totals, seed and deterministic state checksum.

## Headless validation

```bash
godot --headless --editor --path . --quit
godot --headless --path . --script tests/run_tests.gd
```

Tests cover numerical conservation, fair simultaneous competition, division accounting, starvation death, exact same-seed history, genome deep-copy inheritance, parental/sibling immutability, exact mutation replay, seed-contingent mutation history, neutral-marker physiology, configured mutation frequency and mutation ancestry records.

## Documentation

- `docs/ARCHITECTURE.md` — system boundaries, tick semantics, target genome/protein/chemistry architecture, persistence, analytics and performance rules.
- `docs/SCIENTIFIC_MODEL.md` — biological abstractions, equations, assumptions, limitations and validation hierarchy.
- `docs/M3_GENETICS.md` — exact M3 gene/genome representation, inheritance, mutation semantics, neutrality gate and M4/M5 bridge.
- `docs/ROADMAP.md` — gate-driven implementation from M0 through open evolution, laboratory tooling, evolved mutation rate, motility/information and horizontal transfer.
