# ADR 0001 — Deterministic phased simulation instead of per-cell sequential stepping

- Status: Accepted
- Date: 2026-08-22

## Context

Microverse is intended to study evolutionary outcomes. Small implementation asymmetries can therefore become amplified by selection and be mistaken for biology. A naïve object-oriented loop such as:

`for cell in cells: cell.step(world)`

lets earlier cells mutate the shared environment before later cells observe it. Under resource scarcity, array order then becomes an unmodelled competitive trait.

Rendering order, insertion order, ID order, hash/dictionary iteration, and thread scheduling must not accidentally become sources of selective advantage.

## Decision

The simulation uses a deterministic phased tick. Shared-state operations are split into observation/request and resolution/application phases.

For membrane uptake:

1. diffuse environment;
2. all cells sample the same environmental state;
3. all cells submit uptake requests;
4. requests are aggregated by chemical and spatial site;
5. if demand exceeds supply, each request receives the same proportional scaling factor;
6. allocated amounts are removed and applied to cells;
7. intracellular physiology proceeds.

Future shared interactions must follow the same principle where simultaneous semantics are biologically intended. Examples include extracellular enzyme access, toxin exposure, public-good uptake, and physical collision resolution.

Stochastic processes use the simulation-owned seeded RNG. If an algorithm deliberately uses asynchronous random updates, that choice must be explicit, seeded, documented, and tested for sensitivity because asynchronous dynamics can itself change ecological/evolutionary outcomes.

## Consequences

### Positive

- Removes a major hidden selection artifact.
- Makes identical competitors symmetric under identical local conditions.
- Improves reproducibility and causal interpretation.
- Makes future parallelization conceptually safer for commutative phases.

### Costs

- Requires request buffers and aggregation.
- Some processes are less convenient than a single `cell.step()` method.
- The engine must distinguish biologically simultaneous and deliberately asynchronous mechanisms.

## Rejected alternatives

### Fixed sequential cell update

Rejected because stable cell/insertion order can become a persistent resource-access advantage.

### Shuffle cells each tick and update sequentially

Better than fixed order, but still makes resource allocation depend on sampled update sequence. This can be a valid biological modeling choice for some asynchronous processes, but it is not appropriate as an accidental default for simultaneous uptake.

### Give every cell its requested resource and clamp environment afterward

Rejected because it creates matter when aggregate demand exceeds availability and changes competition strength unpredictably.

## Validation

The M0-M2 test suite includes a deliberately supply-limited case with two identical cells in the same lattice site. Both must receive equal substrate and aggregate uptake must not exceed the local available pool.
