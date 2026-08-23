# M5-B — Evolvable regulation and generic chemical sensing

## Purpose

M5-A established a physically persistent molecular phenotype:

`DNA -> stochastic mRNA cohorts -> stochastic protein cohorts -> catalysis`.

M5-B allows those physical protein cohorts to alter future transcription and to be modulated by intracellular chemistry. The engine still does not contain named biological responses.

The causal chain becomes:

`environment -> transport/metabolism -> intracellular molecule concentrations -> protein allostery -> promoter occupancy -> stochastic transcription -> mRNA/protein cohorts -> catalysis -> physiology -> reproduction`.

The aim is not to prescribe a useful response. The aim is to make useful, useless, harmful and neutral regulatory architectures all physically possible under the same generic rules.

## 1. Inherited regulatory genotype

Each `Gene` now contains:

- stable `locus_id`;
- basal `promoter_code`;
- `protein_signature` controlling the sequence identity of newly synthesized molecular cohorts;
- `regulatory_signature`, a 16-bit promoter-binding motif;
- `neutral_marker`.

The regulatory signature participates in canonical genotype identity, fingerprints and checksums.

A regulatory mutation is an ordinary one-bit mutation of this motif. It does not directly set expression, fitness, strategy or environmental response.

## 2. Regulatory binding landscape

For a physical protein cohort with signature `p` and a target promoter motif `m`, the sequence distance is Hamming distance:

`d(p,m) = popcount(p XOR m)`.

A cohort can bind only if:

`d <= regulatory_max_distance`.

Within that radius, binding affinity is:

`A_reg = exp(-regulatory_distance_decay * d)`.

The normalized abundance of the cohort is:

`P = protein_amount / expression_reference_protein_count`.

Before allosteric modulation, promoter occupancy contribution is proportional to:

`P * A_reg`.

### Activator/repressor sign

M5-B uses one digital-chemistry convention:

- signature high bit clear -> activating contribution;
- signature high bit set -> repressing contribution.

This is not a `TranscriptionFactor` class and does not encode a named biological function. Any protein cohort can regulate any promoter if sequence geometry permits.

For all compatible cohorts:

`activation = sum(positive occupancy)`

`repression = sum(negative-sign occupancy)`.

The bounded net signal is:

`net = (activation - repression) / (1 + activation + repression)`.

The target transcription multiplier is:

`R = clamp(1 + regulatory_gain * net, regulatory_min_factor, regulatory_max_factor)`.

M5-A transcription intensity is then changed from:

`lambda_tx = k_tx,max * promoter_strength * dt`

to:

`lambda_tx = k_tx,max * promoter_strength * R * dt`.

Every subsequent transcription event is still sampled by the canonical deterministic Poisson RNG. Regulation changes propensity; it does not bypass stochastic expression or expression cost.

## 3. Dormant ancestral topology

The M5-B ancestor is deliberately not born with a functional regulatory network.

Every ancestral promoter motif is at least Hamming distance 4 from every ancestral protein sequence, while the active regulatory radius is 3.

Therefore:

- ancestral M5-A basal phenotype remains the M5-B baseline;
- no promoter receives an inherited regulatory advantage merely because the feature exists;
- at least some regulatory opportunities are exactly one motif mutation outside the active radius.

This creates a mutationally accessible regulatory landscape rather than a predesigned regulatory circuit.

## 4. Generic ligand identity

Every modeled metabolite carries an independent 16-bit ligand signature.

Examples include ordinary metabolic species and the semantically neutral compound `X`. The ligand signature does not encode whether a molecule is a nutrient, toxin, waste, signal or stressor.

For protein signature `p` and molecule ligand signature `l`:

`d_ligand = d(p,l)`.

Allosteric compatibility exists only if:

`d_ligand <= allosteric_max_distance`.

Affinity is:

`A_ligand = exp(-allosteric_distance_decay * d_ligand)`.

For intracellular amount `C`:

`compatible_amount = C * A_ligand`

and bounded occupancy is:

`theta = compatible_amount / (allosteric_km + compatible_amount)`.

M5-B currently uses the strongest compatible ligand occupancy for each protein cohort. This prevents many simultaneous weak ligand matches from producing an unbounded hidden gain.

## 5. Allosteric direction

A second protein-sequence bit determines whether ligand occupancy potentiates or inhibits that protein's regulatory contribution:

- bit 14 clear -> potentiation;
- bit 14 set -> inhibition.

The multiplier is:

`A_allosteric = clamp(1 +/- allosteric_gain * theta, allosteric_min_factor, allosteric_max_factor)`.

The final promoter-occupancy contribution of a compatible protein cohort is therefore:

`occupancy = P * A_reg * A_allosteric`.

This rule acts only on the regulatory contribution in M5-B. It does not directly change catalytic activity. Later milestones may generalize allostery to catalysis only if that extension passes a separate causal gate.

## 6. Why a molecule can become information

The engine never asks whether a molecule is a signal.

A molecule becomes usable information only if a genotype creates a causal chain such as:

`molecule concentration -> compatible protein allostery -> promoter occupancy -> altered transcription -> altered proteome -> altered physiology`.

The same molecule may be:

- ignored by one genotype;
- useful to another;
- harmful to a third;
- correlated with the environment but evolutionarily irrelevant.

`X` is intentionally semantically neutral and is used in controlled tests to demonstrate this principle.

Calling an evolved `X` response "sensing" is therefore an interpretation of the causal mechanism, not an engine instruction.

## 7. Protein-cohort history and finite memory

M5-A protein cohorts retain the coding signature with which they were synthesized.

Suppose a locus changes from DNA sequence `A` to `B` at division. The daughter can transiently contain:

- inherited protein cohort `A`;
- inherited mRNA cohort `A`;
- new mRNA cohort `B`;
- later new protein cohort `B`.

M5-B regulation reads those physical cohort signatures, not the daughter's current DNA as a shortcut.

Therefore an ancestral protein can continue to regulate a promoter after the encoding DNA has mutated. Its effect disappears only through ordinary turnover or dilution.

This produces finite history dependence without a `memory` variable.

A statement such as "the cell remembers" is justified only when a controlled experiment demonstrates that two cells with the same current external condition and genome have different responses because their molecular histories differ.

## 8. No special response APIs

M5-B must not contain APIs such as:

- `respond_to_starvation`;
- `glucose_sensor`;
- `stress_response`;
- `dormancy_program`;
- `quorum_sensing`;
- `cooperation_gene`;
- `fitness_bonus`.

Permitted primitives are molecular:

- sequence signatures;
- binding distance;
- abundance;
- concentration;
- occupancy;
- transcription propensity;
- mutation;
- turnover.

Higher-level biological descriptions may be assigned only after observing dynamics generated by those primitives.

## 9. Causal validation gates

M5-B acceptance requires all prior milestone suites plus specific controls proving:

1. **Ancestral neutrality** — enabling M5-B does not create an ancestral regulatory edge or accidental ligand match.
2. **Mutational accessibility** — a one-bit promoter-motif mutation can cross the inactive/active affinity threshold.
3. **Bidirectional regulation** — generic compatible proteins can increase or decrease transcription propensity.
4. **Genotype integrity** — regulatory mutation changes the inherited motif without rewriting other gene fields or the parent genome.
5. **Chemical inducibility** — changing only a compatible molecule concentration changes promoter regulation.
6. **Semantic neutrality** — the same molecule has no effect when generic allostery is disabled.
7. **Causal propagation** — a molecule-conditioned regulatory change propagates to protein abundance and downstream reaction flux/biomass outcome.
8. **Molecular history** — inherited old-sequence proteins retain their regulatory effect after DNA sequence change.
9. **Deterministic ordering** — active regulation/allostery remains independent of gene-array ordering and exactly replayable under the same RNG state.

These gates show capability. They do not by themselves prove that adaptive regulatory networks evolve spontaneously under open evolution.

## 10. Remaining M5 research gate

The roadmap's final M5 criterion is stronger than demonstrating an inducible mechanism:

> regulation can evolve differently in stable versus fluctuating test environments.

That requires a controlled evolutionary experiment rather than another unit test.

Before entering M6, Microverse should therefore run a compact replicated M5-C experiment with at least:

- one stable chemical environment;
- one temporally fluctuating version with the same long-term mean availability;
- identical ancestral genotype distribution;
- multiple deterministic seeds per condition;
- lineage/genotype frequency tracking;
- regulatory-edge and expression-variance metrics;
- no explicit fitness score.

A difference between treatments must be evaluated from inherited regulatory architecture and reproductive outcome, not from manually labeled strategy classes.

M6 physical-cell mechanics should begin only after M5-B implementation gates are green and the scope of this M5-C evolutionary gate is explicitly resolved.
