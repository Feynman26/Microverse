# M5-C final capability confirmation: regulated expression burden

## Why this assay exists

Earlier M5-C experiments targeted endogenous R03 oxidative phosphorylation and R09 ROS control. Their molecular mechanisms worked, but independent reproductive panels showed weak or null genotype-by-environment effects. Those negative results are retained as falsifications; they are not reinterpreted as passes.

The M5 exit contract is narrower: demonstrate that generic inherited regulation can produce an environment-dependent reproductive outcome without explicit fitness. The final assay therefore isolates the selectable cost of gene expression itself rather than depending on the ecological leverage of one endogenous reaction.

## Architectures

Both genomes contain the same metabolic core, the same O2-compatible regulator, and three extra high-expression loci encoding the same deliberately non-catalytic protein (`0xF13D`). The burden signature is outside the active radius of every M4 reaction.

- **constitutive:** the three burden promoters have dormant motifs and are not regulated by the O2-compatible repressor;
- **responsive:** only the three burden regulatory motifs differ; they bind the O2-compatible repressor.

No fitness term or genotype-specific cost exists. Both architectures pay ordinary ATP, AA/NUC, finite-proteome, and shared-ribosome costs. In high O2, allostery inhibits the repressor and the architectures should be approximately matched. During anoxia, the responsive architecture reduces neutral expression and can redirect ordinary cellular resources.

## Environment

- Stable treatment: O2 = 6 continuously.
- Fluctuating treatment: O2 = 6 / 0 square wave, 40 min per phase.
- Glucose, nitrogen and phosphorus remain maintained reservoirs.
- Four phase offsets are averaged per seed.
- Mutation is disabled through the standard M5-C assay configuration.

## Frozen confirmation panel

Seeds `24001..24024` are reserved for this confirmation and were not used in earlier M5-C panels.

For each seed, reproductive rate is `ln(2) / first-division-time`. The analysis-only response is:

`Dnorm = [(R-C)_fluctuating - (R-C)_stable] / mean_baseline_growth_rate`

This value never enters physiology or reproduction.

## Predeclared exit criteria

All must pass:

1. every stable and fluctuating lineage reaches ordinary division before the fixed timeout;
2. the mean stable-high responsive-vs-constitutive growth-rate difference has absolute magnitude <= 0.03 in native rate units;
3. mean `Dnorm >= 0.05`;
4. paired one-sample `t >= 2.07` across 24 independent seeds;
5. at least 16/24 seeds have `Dnorm > 0`.

These thresholds and seeds are frozen before reading the panel output.

## Interpretation boundary

Passing this gate demonstrates a selectable consequence of generic environment-conditioned regulation and costly expression. It does **not** claim that the burden construct is an evolved natural adaptation or that a particular ecological strategy has emerged. Open-ended selection and ecological discovery remain later milestones, especially M7-M8.
