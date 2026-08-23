# Microverse

Microverse is a small, reproducible artificial-life laboratory focused on autonomous evolution in a spatially explicit digital microbial world.

The engine defines **physical and chemical possibility**; it does not assign fitness, cooperation, competition, dormancy, predation, signalling or specialization. Those labels belong to later analysis of outcomes, not to cell-state flags.

## Design principles

1. **No explicit fitness score.** Reproductive success is the evolutionary consequence.
2. **Deterministic stochasticity.** Same complete state + RNG state + interventions reproduces the same trajectory.
3. **Simulation/render separation.** Biology advances independently of UI frame rate.
4. **Conservation first.** Numerical or material inconsistencies fail loudly.
5. **Mechanisms, not behaviours.** Molecular rules exist; named strategies do not.
6. **Observability.** Births, deaths, mutations, lineages and molecular transitions must be reconstructable.
7. **Small world, deep biology.** Scope is deliberately narrow enough to inspect causality.

## Current development vertical slice: M0-M5

The branch under validation implements:

`environment -> uptake -> intracellular ligand state -> regulation -> mRNA -> protein cohorts -> catalytic affinity -> reaction flux -> ATP/precursors/BIO -> division -> inheritance -> mutation -> delayed phenotype replacement`

Core capabilities now include:

- closed 2D microchamber and diffusing C/O/N/P resources;
- numerically stable conservative diffusion;
- simultaneous fair allocation of scarce local resources;
- 12-locus heritable genome with stable locus IDs;
- mutable promoter strength, protein signature, regulatory motif and neutral marker;
- deterministic genotype fingerprints and mutation genealogy;
- 19 intracellular digital metabolites and 12 balanced reactions;
- Hamming-distance catalytic promiscuity;
- simultaneous order-independent metabolic flux allocation;
- explicit ATP/ADP and NAD/NADH currencies;
- C/N/P-conserving precursor and BIO assembly;
- cell volume derived from structural BIO;
- aerobic, oxygen-independent and dormant waste-recovery metabolic possibilities;
- explicit ROS production/detoxification and damage;
- **mRNA abundance per locus**;
- **persistent protein cohorts identified by the signature present when they were translated**;
- transcription, translation and first-order molecular turnover;
- ATP, AA and NUC opportunity cost of expression;
- material recycling on macromolecular turnover;
- C/N/P accounting extended across free metabolites + mRNA/protein material;
- generic protein-to-promoter activation/repression by signature affinity;
- generic intracellular ligand allostery with no named sensor classes;
- ligand signatures deliberately outside the ancestral allosteric radius;
- seeded expression noise and stochastic but conservative molecular partition at division;
- genotype-to-phenotype lag after coding mutation because inherited old protein cohorts persist until turnover;
- observer UI exposing genotype, metabolite pools, expression burden, protein-cohort count and dominant flux;
- Godot 4.7.2 headless CI gates.

### Why protein cohorts matter

A DNA mutation must not retroactively mutate proteins that already exist.

M5 therefore tracks protein abundance by both locus and molecular signature. If a daughter inherits protein encoded by signature `A` and its DNA mutates to signature `B`, its immediate state is:

`DNA = B`

`proteome = inherited A + newly translated B`

Over time, A decays and B accumulates. Consequently mutations can have delayed, buffered or transient effects without a dedicated phenotype-delay rule.

### Why chemical sensing is not a behaviour API

Each digital metabolite has a 16-bit ligand signature. A protein can become allosterically sensitive to a molecule only when mutation brings its protein signature into the configured binding neighborhood. The engine does not know that glucose is food, ROS is stress, or `X` is a signal.

A molecule becomes **information** only when evolved binding and regulatory wiring make its concentration predictive/useful downstream.

Controlled M5 assays include a semantically neutral `X` compound that changes biomass-reaction flux only after a constructed molecular circuit makes X-compatible protein binding regulate the relevant enzyme. This validates possibility without hard-coding the ecological interpretation.

## Requirements

- Godot 4.7.2, standard build; no .NET required.

## Run

Open `project.godot` in Godot and run the project.

Controls:

- `Space`: pause/resume;
- `1`: 1x experimental clock;
- `2`: 10x;
- `3`: 100x;
- `4`: 1000x;
- `5`: 5000x.

## Headless validation

```bash
godot --headless --editor --path . --quit
godot --headless --path . --script tests/run_tests.gd
godot --headless --path . --script tests/m4_landscape_tests.gd
godot --headless --path . --script tests/m4_metabolism_tests.gd
godot --headless --path . --script tests/m5_expression_tests.gd
godot --headless --path . --script tests/m5_sensing_tests.gd
```

M5 is not considered complete merely because these new suites pass. All earlier M0-M4 gates must remain green as continuity constraints.

## Documentation

- `docs/ARCHITECTURE.md` — system boundaries and simulation semantics.
- `docs/SCIENTIFIC_MODEL.md` — abstractions, equations, assumptions and validation hierarchy.
- `docs/M3_GENETICS.md` — gene/genome representation and mutation semantics.
- `docs/M4_CATALYTIC_LANDSCAPE.md` — catalytic landscape design/statistical gate.
- `docs/M4_METABOLIC_INTEGRATION.md` — explicit chemistry, flux and BIO growth.
- `docs/M5_EXPRESSION_REGULATION.md` — expression, protein cohorts, regulation, allostery, sensing and phenotype lag.
- `docs/ROADMAP.md` — gate-driven path through spatial ecology, open evolution, experimental tooling, evolved heredity, motility and horizontal transfer.
