# M5-C final capability confirmation: regulated expression burden

## Why this assay exists

Earlier M5-C experiments targeted endogenous R03 oxidative phosphorylation and R09 ROS control. Their molecular mechanisms worked, but independent reproductive panels showed weak or null genotype-by-environment effects. Those negative results are retained as falsifications; they are not reinterpreted as passes.

The M5 exit contract is narrower: demonstrate that generic inherited regulation can produce an environment-dependent reproductive outcome without explicit fitness. The final assay therefore isolates the selectable cost of gene expression itself rather than depending on the ecological leverage of one endogenous reaction.

## Architectures

Both genomes contain the same metabolic core, the same O2-compatible regulator, and three extra high-expression loci encoding the same deliberately non-catalytic protein (`0xC63F`). This burden signature was chosen only after the mechanics gate rejected an earlier signature that accidentally self-regulated. `0xC63F` is outside the active radius of every M4 reaction and outside the relevant ligand/promoter binding radii used by this construct.

- **constitutive:** the three burden promoters have dormant motifs and are not regulated by the O2-compatible repressor;
- **responsive:** only the three burden regulatory motifs differ; they bind the O2-compatible repressor.

No fitness term or genotype-specific cost exists. Both architectures pay ordinary ATP, AA/NUC, finite-proteome, and shared-ribosome costs. In high O2, allostery inhibits the repressor and the architectures should be approximately matched. During anoxia, the responsive architecture reduces neutral expression and can redirect ordinary cellular resources.

## Environment

- Stable treatment: O2 = 6 continuously.
- Fluctuating treatment: O2 = 6 / 0 square wave, 40 min per phase.
- Glucose, nitrogen and phosphorus remain maintained reservoirs.
- Four phase offsets are averaged per seed.
- Mutation is disabled through the standard M5-C assay configuration.

## First burden panel: informative but not confirmatory

Seeds `24001..24024` were frozen before their output was read. The panel produced:

- mean `Dnorm = +0.102924`;
- median `Dnorm ≈ +0.170776`;
- 17/24 positive seed effects;
- mean stable-high advantage `-0.00061652` native rate units.

However, some fluctuating lineages did not reach division, so first-division outcomes contained genuine zero/censored reproductive rates and the distribution was strongly non-Gaussian (including a large negative tail). The predeclared one-sample t criterion therefore failed (`t = 1.592926`) even though the material-effect and sign criteria passed. The original panel remains a failed predeclared gate; its t threshold is not retroactively changed.

This result motivates a separate robust confirmation rather than reinterpreting the first panel.

## Final robust confirmation — frozen before observation

Seeds `25001..25024` are reserved for the final independent confirmation and have not been used by earlier M5-C panels.

For each seed, reproductive rate remains `ln(2) / first-division-time`. Failure to divide before death/timeout remains a zero reproductive rate; it is not discarded. The analysis-only response remains:

`Dnorm = [(R-C)_fluctuating - (R-C)_stable] / mean_baseline_growth_rate`

This value never enters physiology or reproduction.

Because the first panel demonstrated that the response distribution is non-Gaussian/censored, the final confirmation uses robust location and sign criteria rather than a normal-theory t statistic.

All of the following are predeclared before reading any `25001..25024` output:

1. both architectures divide in every stable-high control seed;
2. mean stable-high responsive-vs-constitutive growth-rate difference has absolute magnitude <= 0.03 native rate units;
3. mean `Dnorm >= 0.05`;
4. median `Dnorm >= 0.05`;
5. at least 17/24 independent seeds have `Dnorm > 0`.

For 24 untied Bernoulli signs under a 50/50 null, `P(X >= 17) ≈ 0.03196` one-sided. Requiring 17 positive seeds is therefore a predeclared exact sign-test-level criterion while remaining robust to the heavy-tailed magnitude distribution exposed by the first panel.

## Interpretation boundary

Passing this gate demonstrates a selectable consequence of generic environment-conditioned regulation and costly expression. It does **not** claim that the burden construct is an evolved natural adaptation or that a particular ecological strategy has emerged. The endogenous R03/R09 null results remain scientifically meaningful. Open-ended selection and ecological discovery remain later milestones, especially M7-M8.
