# M5-C ribosome-limited confirmation

## Why this mechanism was added

The independent reproductive-success panel using seeds `21001..21064` showed only a +1.953 percentage-point responsive interaction under the high/anoxic environment and did not clear the predeclared McNemar or 3:1 discordance gates. Stable-high reproduction was 64/64 for both architectures. The result therefore falsified the claim that finite proteome occupancy plus expression ATP/material costs alone generate a sufficiently strong environment-dependent reproductive advantage.

The remaining biological omission was a shared translation-machinery constraint. Before this change, every mRNA cohort could propose translation independently; ATP and amino acids were shared, but ribosomal throughput itself was unlimited. Consequently, unnecessary expression could consume proteome space and material while still allowing all other mRNAs to translate in the same interval.

M5-C now adds a generic finite ribosome budget. This is not specific to R03, oxygen, or either experimental genotype. All mRNA cohorts compete simultaneously for the same translation throughput, and proposals are scaled proportionally from a common snapshot so iteration order cannot create priority.

## Frozen global parameter

`translation_capacity_fraction_of_proteome_per_min = 0.06`

With the current finite proteome capacity:

- maximum proteome = `5.0 * 160 = 800` protein units;
- shared translation capacity = `800 * 0.06 = 48` translation events/min;
- basal protein decay is `0.05/min`, corresponding to approximately 40 replacement events/min at full proteome occupancy.

The remaining throughput is therefore limited. A lineage that expresses unnecessary protein consumes capacity that cannot simultaneously be used to rebuild or maintain other proteins. Conversely, down-regulating a protein can free ribosomal capacity, but restoring that protein after an environmental transition still requires transcription, translation and time.

## Mechanistic gates

Before any selection claim is evaluated, CI must demonstrate that:

1. abundant mRNA cannot exceed shared translation throughput;
2. rejected translation proposals do not consume ATP or amino-acid material;
3. overexpression of one locus reduces translation available to another locus under a binding ribosome budget;
4. allocation remains independent of genome-array ordering;
5. all prior M0-M5-B, RNG, metabolic, expression and finite-proteome gates remain green.

## Independent confirmation panel

The final confirmation seeds are frozen as `22001..22064`. None was used in prior M5-C panels.

The biological construct and environment remain unchanged from the preceding reproductive-success panel:

- constitutive and responsive genomes are isogenic except for the inherited R03 regulatory motif;
- mutation is disabled during the assay;
- stable-high condition uses O2 = 6;
- fluctuating condition alternates O2 = 6 and O2 = 0 every 40 min;
- four phase offsets are evaluated per seed;
- reproductive success means reaching ordinary `ready_to_divide()` criteria before death;
- death before first division is a real reproductive failure, not censored data.

## Predeclared acceptance criteria

All six criteria must pass without modification after observing seeds `22001..22064`:

1. both architectures reproduce in at least 95% of stable-high trials;
2. stable-high reproductive-success difference is at most 3 percentage points in absolute value;
3. the environment-dependent responsive interaction is at least +4 percentage points;
4. the fluctuating paired assay contains at least 8 discordant outcomes;
5. continuity-corrected signed McNemar z is at least 1.96;
6. responsive-only reproductive successes are at least 3 times constitutive-only successes.

If this panel fails, the failure remains part of the M5-C evidence history. Seeds or thresholds must not be changed to rescue the same mechanism.