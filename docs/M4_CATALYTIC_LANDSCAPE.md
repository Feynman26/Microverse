# M4 catalytic landscape and digital chemistry

## Status

This document defines the M4-A research gate. The candidate chemistry and affinity landscape are implemented, but **must not control basal physiology until the characterization tests pass**.

## 1. Why M4 needs a landscape gate

A digital evolutionary system can fail even when mutation and selection work perfectly if genotype-to-function geometry is poor.

Two pathological extremes are especially dangerous:

1. **Disconnected function space** — almost every mutation destroys function and genuinely new catalytic functions are effectively unreachable.
2. **Universal promiscuity** — almost every protein catalyses almost every reaction, so specialization and innovation have little meaning.

Microverse therefore treats the catalytic landscape itself as a scientific object that requires characterization before use.

## 2. Digital chemistry scope

The M4 candidate defines 18 primitive species.

### Structural carbon/nutrient species

- `G` — six-carbon-unit external carbon source;
- `C3` — three-carbon-unit intracellular intermediate;
- `C2` — two-carbon-unit intracellular intermediate;
- `W1`, `W2` — alternative three-carbon-unit waste/intermediate compounds;
- `CO2` — one-carbon-unit oxidized waste;
- `NH4` — one nitrogen structural unit;
- `P` — one phosphorus structural unit;
- `AA` — biomass precursor containing C2/N1;
- `LIP` — biomass precursor containing C4;
- `NUC` — biomass precursor containing C2/N1/P1.

### Energetic/redox/stress species

- `ATP`, `ADP` — high/low energy currency;
- `NAD`, `NADH` — oxidized/reduced redox currency;
- `O2` — terminal electron acceptor;
- `ROS` — oxidative stress currency.

### Semantically neutral ecological species

- `X` — a diffusible compound with no built-in role. Future sensing may give it informational relevance, but the engine never labels it a signal.

## 3. Structural conservation

The model does not claim these abstract species are real molecules. Nevertheless, structural matter cannot be created by reactions.

Each metabolite has explicit digital C/N/P units. Every `ReactionDefinition` calculates product-minus-substrate balance for those units and refuses an unbalanced reaction.

ATP/ADP and NAD/NADH are deliberately tracked as energetic/redox currencies rather than included in C/N/P structural bookkeeping. This keeps coarse-grained reactions interpretable while preventing accidental creation of biomass carbon, nitrogen or phosphorus.

## 4. Candidate reaction network

### R01 — carbon activation

`G + 2 ADP + 2 NAD -> 2 C3 + 2 ATP + 2 NADH`

Provides the first carbon/energy entry route.

### R02 — oxidative carbon processing

`C3 + NAD + ADP -> C2 + CO2 + NADH + ATP`

Converts carbon intermediate into a biomass-useful C2 pool while releasing oxidized carbon.

### R03 — oxidative phosphorylation

`NADH + O2 + 3 ADP -> NAD + 3 ATP + ROS`

High energetic yield coupled to oxidative stress.

### R04 — fermentative redox relief

`C3 + NADH + ADP -> W1 + NAD + ATP`

Lower-yield oxygen-independent redox route. The ancestor has only weak activity.

### R05 — W1 recovery

`W1 + NAD + ADP -> C2 + CO2 + NADH + ATP`

Absent in the ancestor but deliberately one coding mutation away from the active catalytic radius of ancestral locus 4.

### R06 — amino precursor synthesis

`C2 + NH4 + ATP -> AA + ADP`

### R07 — lipid precursor synthesis

`2 C2 + 2 ATP -> LIP + 2 ADP`

### R08 — nucleotide precursor synthesis

`C2 + NH4 + P + 2 ATP -> NUC + 2 ADP`

### R09 — oxidative damage control

`ROS + ATP -> ADP`

Represents costly detoxification rather than a free stress-resistance trait.

### R10 — W2 overflow

`C3 + NADH -> W2 + NAD`

Weak ancestral overflow/redox relief with no direct ATP gain.

### R11 — W2 oxidative recovery

`W2 + O2 + ADP -> C2 + CO2 + ATP + ROS`

Absent ancestrally but one coding mutation from the active radius of ancestral locus 9.

## 5. Protein/reaction signature landscape

Both coding products and reactions use 16-bit digital signatures.

Catalytic distance is Hamming distance:

`d = popcount(protein_signature XOR reaction_signature)`

Affinity is:

- `0` when `d > 4`;
- `exp(-0.70 * d)` when `d <= 4`.

Thus:

- distance 0: affinity 1.000;
- distance 1: ~0.497;
- distance 2: ~0.247;
- distance 3: ~0.122;
- distance 4: ~0.061;
- distance 5+: inactive.

This geometry is intentionally continuous locally while retaining sparse global activity.

## 6. Constitutive M4 bridge

Before M5 implements actual transcription/translation, M4 estimates a gene's constitutive catalytic contribution as:

`activity = promoter_strength * affinity * reaction_catalytic_ceiling`

This is a temporary bridge, not the final expression model.

It has two advantages:

- promoter mutations can eventually alter enzyme abundance continuously;
- coding mutations alter affinity through sequence geometry rather than direct phenotype labels.

M5 replaces promoter-as-abundance with explicit costly mRNA/protein dynamics.

## 7. Ancestral geometry

The candidate reaction signatures are positioned relative to the M3 ancestor.

### Exact ancestral catalytic matches

- locus 1 -> R01;
- locus 2 -> R02;
- locus 3 -> R03;
- locus 5 -> R06;
- locus 6 -> R07;
- locus 7 -> R08;
- locus 8 -> R09.

These are distance 0.

### Weak ancestral side routes

- locus 4 -> R04 at distance 3;
- locus 9 -> R10 at distance 3.

They are active but substantially weaker than exact matches.

### Dormant but accessible innovations

- locus 4 -> R05 at distance 5;
- locus 9 -> R11 at distance 5.

Each is initially inactive. Flipping one of the five differing bits moves distance from 5 to 4 and creates low but nonzero activity. Selection can then, in principle, amplify the new activity through subsequent coding/promoter evolution.

This gives Microverse a concrete test of innovation from side-function geometry rather than a special `unlock_reaction` mutation.

## 8. Statistical landscape gate

CI samples 4096 random protein signatures against all 11 reaction signatures.

Acceptance window:

- active random protein/reaction pairs: 2.5%–5.5%;
- mean active reactions per random protein: 0.20–0.80;
- no sampled universal catalyst;
- all one-bit neighbors of an exact match retain nonzero activity;
- both dormant waste-recovery routes become active after one appropriate bit flip.

These bounds are deliberately broad enough to tolerate sampling noise but narrow enough to reject a qualitatively different landscape.

The theoretical probability that a random 16-bit signature lies within Hamming distance 4 of one fixed reaction signature is approximately 3.84%, matching the intended sparse regime.

## 9. What passing this gate means

Passing M4-A does **not** prove autonomous metabolic innovation will occur. It proves the function landscape has three necessary properties:

1. local mutations can preserve partial function;
2. random proteins are usually inactive for any one reaction;
3. new reactions can be mutationally accessible without being pre-active.

Only after this gate passes should the current hard-coded respiration/growth scaffold be replaced by reaction flux generated from the evolving genome.

## 10. M4-B integration plan

The next vertical slice will introduce a dictionary/packed-array intracellular metabolite state containing the M4 species and replace the compressed respiration function with reaction-network execution.

Integration order:

1. add intracellular `MetabolicState` with nonnegative pools;
2. preserve compatibility tests while migrating glucose/O2/ATP/ROS;
3. initialize ADP/NAD and imported NH4/P pools;
4. calculate effective reaction capacities from genome activities;
5. execute reactions under substrate limitation;
6. create ATP and biomass precursor pools through the network;
7. replace single `precursor` growth with explicit AA/LIP/NUC requirements;
8. remove `respiration_scale` and `growth_scale` from selectable physiology;
9. add controlled oxygen-rich versus oxygen-poor comparisons;
10. only then expose coding/promoter mutations to actual selection.

## 11. Required safeguards during M4-B

- reaction scheduling must not create a systematic order advantage between reactions;
- substrate competition needs either sufficiently small substeps or a simultaneous flux-allocation scheme;
- no pool may become negative;
- C/N/P structural totals must change only through explicit environment exchange, biomass accounting or lysis;
- ATP production and consumption must be explicit;
- oxygen-independent metabolism must remain lower-yield rather than a free equivalent of respiration;
- ROS trade-off must remain coupled to oxidative ATP generation;
- M3 neutral markers remain functionless controls.
