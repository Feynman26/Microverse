# M5-A — Explicit stochastic gene expression

## Purpose

M5-A removes the last direct bridge from DNA promoter value to catalytic capacity. The production causal chain is now:

`DNA promoter -> stochastic transcription -> mRNA cohort -> stochastic translation -> protein cohort -> protein/reaction affinity -> metabolic flux`.

A promoter changes how often new transcripts are proposed. It does not directly change reaction rate.

## Why molecular cohorts are required

A simple `protein_count[locus]` representation creates a serious artifact after coding mutation: if the locus protein signature changes, every previously inherited molecule would appear to change sequence instantaneously.

M5-A therefore stores expression state as:

```text
locus_id -> {
  mRNA:    {coding_signature -> amount},
  protein: {coding_signature -> amount}
}
```

The DNA has one current coding signature. Historical molecular cohorts may have older signatures.

After a coding mutation:

```text
mother DNA: A
mother protein cohort: A

          division + mutation
                  |
                  v
 daughter DNA: B
 inherited protein: A
 new transcript: B
 new protein: B

 transient proteome: A + B
```

Old molecules decay according to ordinary turnover. Nothing explicitly implements `phenotypic_lag` or `memory`; the lag follows from molecular persistence.

## Basal kinetic grammar

For gene `i`, basal transcription intensity is:

`lambda_tx,i = k_tx,max * promoter_strength_i * dt`.

For every mRNA sequence cohort `s` at that locus:

`lambda_tl,i,s = k_tl * mRNA_i,s * dt`.

Degradation intensities are:

`lambda_mdeg,i,s = k_mdeg * mRNA_i,s * dt`

and

`lambda_pdeg,i,s = k_pdeg * protein_i,s * dt`.

Discrete event proposals are sampled using the universe's deterministic Poisson RNG. M5-A uses exact Knuth sampling for small intensities and decomposes larger intensities into independent chunks.

The default compressed kinetics are chosen such that the deterministic basal steady-state protein amount approximately satisfies:

`protein_i / expression_reference_protein_count = promoter_strength_i`.

This preserves the already-validated M4 ancestral catalytic scale while replacing its direct promoter shortcut with real intermediate molecular state.

These values are digital kinetic parameters, not claims about quantitative bacterial copy numbers.

## Deterministic stochasticity

All transcription, translation, degradation, division partitioning and mutation draws use the same authoritative universe RNG.

To make replay independent of container ordering:

- genes are traversed by stable locus ID;
- mRNA/protein cohorts are traversed by sorted coding signature;
- synthesis-resource allocation is calculated from one immutable proposal snapshot.

Same seed/state therefore reproduces exact molecular noise; different seeds may separate genetically identical clones.

## Shared synthesis resources

Expression is not free.

Transcription consumes:

- ATP -> ADP;
- a small amount of `NUC` structural precursor transferred into mRNA storage.

Translation consumes:

- ATP -> ADP;
- a small amount of `AA` structural precursor transferred into protein storage.

Degradation currently returns the corresponding NUC/AA material completely to precursor pools. This is a deliberate compressed-recycling assumption for M5-A; partial recovery and waste generation can be introduced later if mechanistically useful.

### Order-independent scarcity

All genes first propose synthesis events. The engine then calculates global demand for ATP, NUC and AA from the same snapshot.

If resources are limiting, one shared scale is computed before any locus is updated:

```text
all loci propose transcription/translation
                |
                v
       total ATP/NUC/AA demand
                |
                v
       shared availability scale
                |
                v
       accepted events per locus
```

No gene obtains synthesis priority because it happens to be first in an array.

## Material accounting

Expression introduces structural matter outside the metabolite dictionary, so conservation must include both compartments:

`modeled structural matter = metabolite pools + mRNA storage + protein storage`.

mRNA storage uses the digital C/N/P structural signature of `NUC`; protein storage uses the digital structural signature of `AA`.

CI verifies that combined C/N/P is unchanged by expression/turnover under the complete-recycling assumption.

ATP is not destroyed by expression. It is converted to ADP, so ATP+ADP is conserved by expression events.

## Catalytic coupling

M5 metabolism never asks a gene for promoter strength.

For every physical protein cohort with signature `p`, reaction `r` receives a contribution proportional to:

`protein_amount / reference_count * affinity(p,r) * catalytic_ceiling_r`.

All historical protein cohorts at the locus contribute according to their own sequence. Therefore a mutant allele can coexist transiently with inherited ancestral catalytic activity.

## Division

mRNA and protein cohorts are inherited by both daughters.

Partitioning:

- exactly conserves total amount across sisters;
- uses the same overall cell partition ratio as a center;
- adds seeded stochastic deviation per molecular cohort;
- attenuates partition noise approximately as `1/sqrt(amount)` so scarce molecules fluctuate more strongly than abundant ones.

This creates immediate phenotypic asymmetry between genetically identical sisters without inventing different cell types.

## Mutation timing

Current division semantics are:

1. parent pays division ATP cost;
2. metabolite and molecular state are partitioned;
3. daughters receive deep-copy genomes;
4. daughter genomes mutate;
5. inherited molecular cohorts remain unchanged;
6. subsequent expression reads the mutated DNA.

This sequence is intentional. A mutation changes future molecular production rather than historical molecules.

## Current observable consequences

The UI reports for a focal cell:

- total mRNA;
- total protein;
- number of protein sequence cohorts;
- transcription/translation accepted in the most recent tick;
- ATP spent on expression;
- central metabolic state and dominant reaction flux.

A protein-cohort count above the number of loci can indicate transient sequence coexistence after mutation.

## Validation gates

M5-A CI proves:

1. M4 equilibrium catalytic scale is preserved at initialization;
2. identical seed reproduces exact expression trajectory;
3. different seeds can separate clone expression phenotypes;
4. combined pools + expression storage conserve digital C/N/P;
5. transcription/translation spend ATP but conserve ATP+ADP;
6. reversing gene-array order does not change stochastic expression history;
7. realized protein abundance controls catalytic flux;
8. coding mutation does not rewrite inherited protein;
9. new mutant protein appears only through new transcription/translation;
10. ancestral protein can decay after mutation;
11. stochastic sister partition conserves molecular amount while creating asymmetry;
12. Poisson sampler ensemble mean is statistically correct;
13. resource scarcity scales synthesis without negative pools;
14. whole-world same-seed replay includes exact expression state.

## Performance note

Explicit expression substantially increases per-tick work: each cell now evaluates stochastic events for every locus and molecular cohort. The current CI baseline remains tractable for the intended small-world phase, but M5-A is the first point where long M8 replicate experiments will eventually require profiling and optimization.

Optimization must preserve model semantics. Candidate future strategies include:

- event aggregation at high molecule counts;
- adaptive expression stepping when rates are slowly varying;
- sparse cohort storage;
- batched headless simulation;
- avoiding UI updates during accelerated/headless runs.

No optimization should change random-event ordering or introduce seed-dependent platform divergence without a model-version boundary.

## M5-B boundary

M5-A has expression but not regulation. Promoter strength remains a basal scalar encoded in the genome.

M5-B adds generic evolvable regulation:

- promoter binding motifs;
- regulator/promoter affinity;
- activator/repressor effects;
- protein abundance-dependent promoter occupancy;
- generic chemical/receptor modulation;
- feedback and finite history dependence;
- mutations of regulatory architecture;
- controlled fluctuating-environment comparisons.

M5-B must not introduce named APIs such as `stress_response`, `dormancy`, `quorum_sensing`, or `memory`. Those labels may only be applied later to dynamics that arise from the generic molecular network.
