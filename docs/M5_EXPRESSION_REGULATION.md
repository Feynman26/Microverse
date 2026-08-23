# M5 — explicit expression, regulation, and phenotypic noise

## Purpose

M5 inserts a transient molecular phenotype between genotype and catalytic function. A gene no longer contributes catalytic activity merely because its promoter code is high. Living cells must transcribe mRNA, translate protein, retain or degrade those molecules, pay their resource costs, and inherit a stochastic fraction of the molecular state during division.

The authoritative causal chain becomes:

`environment -> uptake -> genome -> regulation -> mRNA -> protein -> catalytic affinity -> reaction flux -> ATP/precursors/BIO -> volume -> division`

This distinction is essential. Two cells with identical genomes can now have different protein abundances and therefore different immediate physiology. Conversely, a genotype mutation can change future molecular dynamics without instantaneously rewriting the existing proteome.

## State separation

### Genome

Inherited, discrete information:

- stable `locus_id`;
- `promoter_code` controlling basal transcription propensity;
- `protein_signature` controlling catalytic and regulatory molecular identity;
- `regulatory_signature` controlling what proteins can bind the promoter;
- `neutral_marker` for drift/lineage controls.

### ExpressionState

Transient, per-cell state:

- mRNA abundance per locus;
- protein abundance per locus.

Expression state is not part of genotype identity. It is partitioned at division and can therefore differ between clonal sisters.

## Basal expression kinetics

For locus `i`, ignoring resource limitation and regulation:

`dm_i/dt = k_tx * promoter_i - k_mdeg * m_i`

`dp_i/dt = k_tl * m_i - k_pdeg * p_i`

The default rates are calibrated so that the initial ancestral protein steady state is approximately equal to the former M4 promoter-strength proxy. This preserves continuity across the milestone while removing that proxy from living-cell physiology.

Initialization uses the deterministic basal steady-state estimate. Subsequent dynamics are explicit.

## Expression cost

Positive synthesis is not free.

Transcription consumes:

- ATP, converted to ADP;
- NUC precursor material.

Translation consumes:

- ATP, converted to ADP;
- AA precursor material.

mRNA/protein degradation returns the represented NUC/AA precursor fraction. This is a compressed turnover model, not a literal ribosome/proteasome model. Its purpose is to make expression an opportunity cost while keeping represented material auditable.

When ATP or precursor material is insufficient, all synthesis requests sharing the limiting resource are scaled from the same pre-synthesis snapshot. No locus wins because it appears earlier in an array.

## Generic regulation

M5 does not define named transcription factors, stress genes, dormancy genes, or quorum genes.

Every expressed protein is *capable* of promoter binding if its 16-bit `protein_signature` lies within the configured Hamming radius of a target gene's inherited `regulatory_signature`.

Binding weight is:

`occupancy = protein_abundance * exp(-decay * HammingDistance)`

Proteins with the high signature bit clear contribute activating occupancy; proteins with it set contribute repressing occupancy. This sign convention is a digital-chemistry rule, not a behavior label.

For a target promoter:

`net = (activation - repression) / (1 + activation + repression)`

`regulation_factor = clamp(1 + gain * net, min_factor, max_factor)`

The transcription request is then proportional to:

`promoter_strength * regulation_factor`.

Because the same `protein_signature` also determines catalytic affinities, mutation can have pleiotropic consequences: a coding change may alter catalysis, regulation, or both.

## Ancestral regulatory topology

The ancestor contains a sparse ring of exact motif matches. This is not a hard-coded adaptive circuit. It merely guarantees that the regulatory grammar is exercised from generation zero and provides a measurable baseline against which mutation can alter network topology.

The ring crosses the high-bit boundary, so some edges activate and some repress.

## Stochastic expression

Each transcription and translation request receives a bounded multiplicative noise factor from the single universe RNG. Therefore:

- same complete state + same RNG state -> exact same trajectory;
- different seeds can produce different molecular trajectories for the same genotype;
- no unseeded/random UI source can alter biology.

Direct isolated-cell compatibility tests that intentionally omit the universe RNG run noise-free. Authoritative `SimulationEngine` execution always passes the seeded RNG.

## Division

Metabolite pools continue to use the cell-volume partition ratio.

mRNA and proteins receive a separate seeded local partition perturbation per locus while conserving the exact parent total:

`daughter_A + daughter_B = parent`

for every tracked mRNA and protein abundance.

This creates a mechanistic source of clonal phenotypic divergence without a `phenotype_randomization` API.

If a daughter receives a DNA mutation after partition, inherited parental proteins can temporarily persist even though the new genome encodes a different sequence. This molecular carryover is intentional and biologically analogous to transient maternal/parental protein inheritance.

## M5 acceptance gates in this subphase

The implementation must demonstrate:

1. initial expression preserves the M4 activity scale closely enough to avoid a discontinuous model reset;
2. explicit zero protein abundance yields zero catalytic reaction flux;
3. transcription and translation have ATP/material costs;
4. generic exact-match activators and repressors move transcription in opposite directions;
5. same-seed expression noise is exactly reproducible;
6. molecular partition conserves total mRNA/protein while allowing sister differences;
7. regulatory motifs are inherited and mutably encoded;
8. full simulation chemistry + expression + genetics remains replayable from the same seed.

## Remaining M5 work after this subphase

This layer creates regulation but not yet generic environmental sensing. The next M5 gate should add molecular sensor coupling to the same regulatory machinery so an extracellular/internal chemical state can modulate regulator activity without adding named environmental-response behaviors.

A valid sensor design should preserve the following rule:

> The engine exposes molecules and binding relationships; it never exposes a semantic instruction such as `respond_to_starvation`.

M6 spatial mechanics should not start until M5 sensing/regulation behavior has passed controlled inducibility and reproducibility experiments.
