# M4-B — Evolvable intracellular metabolism

## Purpose

M4-B replaces the provisional M0-M3 rule

`glucose -> fixed respiration scalar -> ATP + precursor -> direct volume growth`

with an explicit causal chain:

`genome -> protein signatures/promoter abundance -> catalytic affinities -> reaction fluxes -> metabolite pools -> BIO -> physical volume -> division`.

No reaction carries a fitness value. A mutation is useful, neutral, or harmful only through the chemistry it changes in the current environment.

## Authoritative intracellular state

Every cell owns one metabolite ledger containing 19 pools:

- external-resource forms: `G`, `O2`, `NH4`, `P`;
- carbon intermediates: `C3`, `C2`;
- reusable wastes: `W1`, `W2`, `CO2`;
- structural precursors: `AA`, `LIP`, `NUC`;
- assembled structural material: `BIO`;
- energetic/redox currencies: `ATP`, `ADP`, `NAD`, `NADH`;
- stress/information substrate: `ROS`, `X`.

`BIO` is not a display statistic. Cell volume is derived from it:

`volume = BIO / biomass_units_per_volume`.

Thus growth is no longer an independent phenotype knob.

## Reaction network

The current M4 network has twelve reactions:

1. R01 carbon activation;
2. R02 oxidative carbon processing;
3. R03 oxidative phosphorylation;
4. R04 fermentative redox relief;
5. R05 W1 recovery;
6. R06 amino precursor synthesis;
7. R07 lipid precursor synthesis;
8. R08 nucleotide precursor synthesis;
9. R09 oxidative damage control;
10. R10 W2 overflow;
11. R11 W2 oxidative recovery;
12. R12 structural biomass assembly.

Every reaction is validated for digital C/N/P structural balance. ATP/ADP and NAD/NADH are modeled as separate energetic/redox currency pairs and are conserved by catalytic reactions.

## Catalysis

For a gene with protein signature `p` and reaction signature `r`, catalytic affinity depends on Hamming distance:

`A(p,r) = exp(-0.70 d)` for `d <= 4`, otherwise `0`.

M4 uses promoter code as a temporary constitutive abundance proxy:

`activity = promoter_strength * affinity * catalytic_ceiling`.

M5 replaces this proxy with explicit transcription, translation, degradation, regulation and stochastic molecular abundance. The catalytic-affinity geometry remains.

## Flux solver

Reaction execution order must not become biochemistry. Each metabolic substep therefore uses a synchronous allocation algorithm.

### Phase 1 — immutable snapshot

All reactions see the same intracellular pool snapshot. Each computes a potential flux from:

- genome-derived catalytic activity;
- cell volume;
- global metabolic rate scale;
- limiting-substrate saturation.

### Phase 2 — shared-substrate allocation

For each metabolite, total requested consumption across all reactions is calculated. If demand exceeds supply, one scale factor is computed from the immutable snapshot.

### Phase 3 — simultaneous application

Each reaction receives the most restrictive scale among its substrates. All substrate/product deltas are accumulated and only then applied to the authoritative pools.

Reversing the reaction-array order is an explicit CI invariant and must not change the resulting state beyond numerical tolerance.

## Environmental resources

M4 exposes four diffusing environmental fields:

- glucose/carbon;
- oxygen;
- reduced nitrogen;
- phosphorus.

Their membrane uptake is intentionally still basal scaffolding. The evolutionary gate in M4 concerns intracellular catalytic metabolism. Evolvable membrane machinery and regulation are introduced in later milestones rather than mixing two major genotype-to-phenotype transitions at once.

Competition between cells for every transported field uses the same simultaneous request/allocation protocol already established in M0-M2.

## Energetic costs outside catalytic metabolism

Maintenance, basal repair, constitutive proteome burden and cell division spend ATP by converting it to ADP. They do not delete ATP from the ledger without a corresponding low-energy currency.

The M4 promoter-to-abundance bridge has an explicit ATP opportunity cost. This prevents promoter mutations that increase expression from being universally free. M5 will replace this compressed cost with explicit mRNA/protein synthesis and turnover.

## ROS and damage

R03 generates ROS as a metabolic consequence of oxidative phosphorylation. R09 can consume ROS at ATP cost. Residual ROS contributes to damage. A small spontaneous decay term remains as nonenzymatic chemistry; basal damage repair remains transitional physiology until repair capacity is made protein-derived.

## Deliberately dormant innovations

R05 and R11 are inaccessible to the complete ancestral proteome, not merely to one designated locus. Their target signatures are exactly five bits away from one ancestral coding signature while the active radius is four.

Therefore one ordinary coding mutation can create low but nonzero activity. No mutation operator knows that the new reaction is useful.

For R05, CI demonstrates the full causal chain:

`one bit flip -> distance 5 to 4 -> nonzero W1 recovery flux -> increased C2`.

Whether that increases descendants depends on environment and competition.

## Conservation boundaries

Within `MetabolicSolver.step`:

- structural C/N/P are conserved;
- ATP + ADP is conserved;
- NAD + NADH is conserved;
- no metabolite may become negative.

At the whole-cell level, ATP is converted to ADP by maintenance, repair, expression burden and division. ROS may decay spontaneously. Death currently returns only directly transported small-resource pools to the environment. General organic lysis/recycling is deferred to M7 and must be implemented explicitly rather than assumed.

## Validation gates

M4-B CI must keep all of the following green:

- baseline simulation/genetic invariants;
- M4 catalytic-landscape characterization;
- reaction-order independence;
- C/N/P conservation;
- ATP/ADP conservation;
- NAD/NADH conservation;
- zero-expression/zero-flux control;
- aerobic versus hypoxic pathway response;
- one-mutation activation of dormant W1 recovery;
- catalytic BIO assembly;
- explicit expression-energy cost;
- resource-supported growth/division;
- extinction without replenishable energy;
- exact same-seed history/checksum reproducibility.

## Known abstractions before M5

M4 is intentionally not yet a quantitative bacterial cell model. Important abstractions remain:

- promoter code is not yet mRNA/protein copy number;
- membrane transport is not yet genome-derived;
- protein synthesis does not yet consume explicit amino-acid material;
- repair capacity is still basal;
- reaction kinetics use a common compressed saturation grammar;
- no gene regulation, signalling or phenotypic noise exists yet;
- no organic secretion/cross-feeding exists yet.

These are explicit dependency boundaries, not hidden omissions. M5 should first replace promoter abundance with costly stochastic gene expression and regulation before additional ecological complexity is introduced.
