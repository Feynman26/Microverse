# M6-C spatial ecology causal gate

M6-A established finite disk geometry and deterministic contact relaxation. M6-B replaced the all-pairs broad phase with a spatial hash while preserving exact contact semantics. M6-C asks the remaining scientific question: does physical position actually change ecology, rather than merely changing where circles are drawn?

## Causal chain under test

The production path is:

`cell position -> local lattice resource -> membrane uptake -> intracellular chemistry/expression -> BIO/ATP/viability -> local division -> mechanics -> descendant positions`

No fitness, colony, sector, competition or resource-preference variable is introduced. Spatial effects arise from the existing local sampling, finite diffusion, local birth and finite-area mechanics.

## Blocking experiments

`tests/m6_spatial_ecology_tests.gd` contains four deterministic capability gates.

### 1. Growing-colony depletion gradient

A mutation-free ancestor grows in an initially uniform finite-resource chamber with slow but nonzero diffusion. After descendants are present, glucose is compared between lattice sites near living cells and sites far from every living cell. Acceptance requires lower near-colony glucose and a materially nonuniform field.

This test is deliberately based on consumption by a growing colony rather than a pre-painted gradient.

### 2. Position-only ecological counterfactual

Two simulations share the same configuration, RNG seed, ancestor genotype, molecular initialization and fixed resource field. The only changed variable is ancestor position relative to a nondiffusing resource patch.

The inside-patch lineage must remain viable and reproduce while the otherwise identical outside-patch lineage reaches ordinary energy-failure extinction. This is the direct M6 causal gate that physical position can change a fixed-horizon ecological outcome.

### 3. Clonal-sector capability

Two neutral-marker founders undergo three controlled local generations using the production daughter-placement and mechanical-relaxation code. The markers have no physiological effect and exist only to identify descendants in the assay.

Acceptance requires exact marker inheritance, no material overlap and two coherent descendant neighborhoods without any kin-attraction, sector-assignment or lineage-grouping rule.

The forced divide-ready states in this test isolate the geometry/local-reproduction capability. They are not a population-growth or fitness experiment.

### 4. Same-seed spatial replay

Two identical localized-resource runs must reproduce the complete checksum, event history, surviving IDs and every final x/y coordinate exactly. Because cell position is already included in the cell/world checksum, this extends the deterministic replay contract through M6 mechanics and local ecological feedback.

## Interpretation boundary

Passing M6-C means space is causally active and reproducible. It does not claim spontaneous ecological diversification, cross-feeding, public goods, cheating, biofilm evolution or regulatory selection. Those interpretations require later M7/M8 mechanisms and replicated experiments.

Adhesion remains optional and must be added only as a generic molecularly controlled physical interaction if a later experiment requires it; M6 does not introduce a named `biofilm` state.
