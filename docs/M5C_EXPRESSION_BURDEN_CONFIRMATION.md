# M5-C expression-burden confirmation

## Purpose

Earlier M5-C experiments targeted endogenous R03 oxidative phosphorylation and R09 ROS control. Their molecular mechanisms worked, but independent reproductive panels showed weak or null genotype-by-environment effects. Those negative results are retained as falsifications.

The final M5 capability assay isolated the selectable cost of gene expression itself rather than depending on the ecological leverage of one endogenous reaction.

## Architectures

Both genomes contain the same metabolic core, the same O2-compatible regulator, and three extra high-expression loci encoding the same deliberately non-catalytic protein (`0xC63F`). This burden signature was chosen only after the mechanics gate rejected an earlier signature that accidentally self-regulated. `0xC63F` is outside the active radius of every M4 reaction and outside the relevant ligand/promoter binding radii used by this construct.

- **constitutive:** the three burden promoters have dormant motifs and are not regulated by the O2-compatible repressor;
- **responsive:** only the three burden regulatory motifs differ; they bind the O2-compatible repressor.

No fitness term or genotype-specific cost exists. Both architectures pay ordinary ATP, AA/NUC, finite-proteome, and shared-ribosome costs. In high O2, allostery inhibits the repressor. During anoxia, the responsive architecture reduces neutral expression and can redirect ordinary cellular resources.

## Environment

- Stable treatment: O2 = 6 continuously.
- Fluctuating treatment: O2 = 6 / 0 square wave, 40 min per phase.
- Glucose, nitrogen and phosphorus remain maintained reservoirs.
- Four phase offsets are averaged per seed.
- Mutation is disabled through the standard M5-C assay configuration.

## First burden panel — failed predeclared gate

Seeds `24001..24024` were frozen before their output was read. The panel produced:

- mean `Dnorm = +0.102924`;
- median `Dnorm ≈ +0.170776`;
- 17/24 positive seed effects;
- mean stable-high advantage `-0.00061652` native rate units;
- `t = 1.592926`, below the frozen `2.07` threshold.

Some fluctuating lineages did not reach division, producing genuine zero/censored reproductive rates and a strongly non-Gaussian distribution. The original panel therefore remained a failed gate; its t threshold was not changed after observation.

## Final robust confirmation — independent falsification

A second panel, seeds `25001..25024`, was frozen before observation. Failure to divide before death/timeout remained a zero reproductive rate and was not discarded.

The predeclared robust criteria were stable-control validity, stable equivalence, mean and median `Dnorm >= 0.05`, and at least 17/24 positive seeds.

Observed result:

- stable paired controls valid: `18/24`;
- mean stable-high advantage: `+0.00012381` native rate units;
- median `Dnorm = -0.096007`;
- positive seeds: `8/24`;
- fluctuating divisions: `122/192` trajectories.

The reported arithmetic mean `Dnorm = +25,665,273.684331` is invalid as an effect summary. In seeds 25009 and 25023 the baseline reference rate approached zero, so normalization generated values near `-2.0e9` and `+2.7e9`. This demonstrates that the first-division normalized ratio is itself numerically unsuitable in regimes containing extinction/non-division.

The robust confirmation therefore fails three of five criteria and, more importantly, reverses the sign prevalence observed in the first panel. The reproductive effect is **not independently confirmed**.

## Final interpretation

The M5 molecular machinery remains validated: regulation can alter transcription, protein abundance, flux, and expression cost as a consequence of generic molecular interactions. What is not established is a reproducible reproductive advantage for the tested regulatory architectures under this M5-scale assay.

No further M5-C tuning or confirmation panels are planned. Population-level selection experiments are deferred to M8, where fixed-horizon counts, survival, lineage expansion, and competition frequencies can be analyzed without dividing by near-zero first-division rates.
