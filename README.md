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

## Initial vertical slice

The first implementation milestone establishes the minimum causal chain:

`environmental nutrient -> membrane uptake -> metabolism -> ATP -> maintenance/biomass -> growth -> division`

It also includes a spatial diffusion field, seeded deterministic RNG, headless simulation, invariant checks and a minimal renderer. Evolution, regulatory networks and evolvable reaction affinity are introduced only after the basal physiology is numerically stable.

Detailed architecture and scientific assumptions live under `docs/` as the implementation is built.
