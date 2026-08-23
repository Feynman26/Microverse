# M6-A physical-cell core

## Purpose

M6 makes space causal. A cell is no longer a dimensionless sampler located at a point: it occupies a finite disk whose size follows structural biomass, collides with neighboring disks, is confined by chamber walls, and is mechanically displaced by crowding.

This first M6 slice deliberately uses an O(N^2) reference contact search for populations up to the current 64-cell cap. Correct, order-independent mechanics are validated before adding a spatial acceleration structure.

## Geometry

The 2D simulation treats structural cell volume as occupied area for mechanics. Radius therefore scales as

`r(V) = r0 * sqrt(V / V0)`

where `V0` is ancestor volume and `r0` is the configured ancestor radius in grid coordinates. The default `r0 = 0.45` gives an ancestral diameter slightly below one chemical lattice spacing while allowing continuous sub-grid positions.

Walls constrain the complete disk:

`r <= x <= width - 1 - r`

`r <= y <= height - 1 - r`

rather than merely clamping the cell center to the field indices.

## Contact solver

Each relaxation iteration is Jacobi/snapshot based:

1. snapshot all living-cell positions;
2. detect all overlapping disk pairs from that same snapshot;
3. compute an equal-and-opposite correction for each pair;
4. accumulate all corrections without moving any cell;
5. apply accumulated displacement simultaneously;
6. project disks back inside chamber walls;
7. repeat until no material overlap remains or the fixed iteration budget is exhausted.

For pair `i,j`, with overlap `delta = ri + rj - distance`, each receives half of a relaxed correction along the pair normal. No cell ID affects correction magnitude, direction, viability, or resource access. IDs are used only to canonicalize floating-point accumulation order for exact replay.

Exact coincident centers are geometrically singular. Production division already supplies a seeded nonzero daughter axis. For externally constructed degenerate states, the solver derives a deterministic fallback direction from pair geometry/location only, not from identity or insertion order.

## Tick integration

After intracellular growth, deaths, and divisions, the simulation relaxes the living population mechanically. Thus a biomass increase changes physical radius and can push neighbors before the next tick samples extracellular resources.

Spatial ecology remains downstream of ordinary chemistry: mechanics does not contain fitness, kin preference, biofilm state, cooperation, competition, or motility rules.

## M6-A validation

The dedicated headless gate checks:

- square-root radius scaling;
- complete-disk wall confinement;
- material-overlap removal;
- equal/opposite symmetry for an isolated pair;
- insertion-order invariance of a controlled multi-cell equilibrium;
- local daughter placement followed by non-overlapping relaxation;
- approximate center-of-volume conservation through division mechanics.

## Remaining M6 work

M6-A is the correctness reference. Remaining M6 scope includes a neighbor/spatial index, explicit tests of colony-generated resource gradients and spatial ecological consequences, lineage-sector visualization, and optional molecularly derived adhesion without a named biofilm state.
