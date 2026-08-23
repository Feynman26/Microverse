# M5-C final reproductive-success confirmation — frozen before output

## Why the endpoint changed

The amplified G×E first-division panel (`20001–20032`) established a material and statistically reproducible environment-dependent regulatory effect, but failed its predeclared requirement that every lineage divide. The observed interaction was `mean D_norm = +0.075556`, `t = 2.750374`, with `25/32` seed-level effects positive.

A diagnostic rerun did not reinterpret that failed gate as success. It showed that the missing divisions were not timeout censoring: cells actually died during the high/anoxic schedule. Across the 128 phase-specific fluctuating pairs, constitutive lineages failed before division more often than responsive lineages. Therefore increasing the timeout or retaining an `all_divided` requirement would discard a real reproductive phenotype.

The M5 exit gate asks whether inherited regulatory architecture changes reproductive outcome when environmental structure changes. Death before reproduction is itself such an outcome. The final assay therefore measures the probability of achieving ordinary first division before death under a fixed observation horizon.

## Biology and circuit are frozen

No production physiology, regulatory parameter, proteome capacity, promoter, protein sequence, environmental concentration, or phase duration is changed after the `20001–20032` result.

The final assay uses the same amplified isogenic competitor pair:

- both genotypes have identical high basal expression of R03 oxidative phosphorylation protein and the O2-compatible regulator;
- both therefore pay the same basal expression and finite-proteome opportunity costs;
- the only inherited difference is the R03 promoter regulatory motif;
- the responsive motif binds the generic O2-compatible repressor; the constitutive motif does not;
- no explicit fitness, survival bonus, behavior flag, or genotype-specific death rule exists.

Environment is unchanged:

- stable reference: constant `O2 = 6`;
- fluctuating treatment: `O2 = 6` / `O2 = 0`, 40 biological minutes per half-cycle;
- fluctuating phase offsets: `0, 200, 400, 600` ticks;
- glucose, nitrogen and phosphorus remain maintained reservoirs;
- mutation remains disabled.

The maximum observation horizon remains 7200 ticks. A lineage is a reproductive success only if ordinary `ready_to_divide()` criteria are reached before death. Death or surviving without division by the horizon counts as reproductive failure; no synthetic penalty is added.

## Independent confirmation panel

Final seeds are **21001–21064**. None appeared in the discovery, null, diagnostic, G×E, amplified G×E, or censored-lineage panels.

For each seed:

1. run constitutive and responsive lineages once in stable high oxygen;
2. run a paired constitutive/responsive lineage at each of the four fluctuating phase offsets;
3. record `success = 1` only when ordinary first-division criteria are reached before death/horizon.

This gives:

- 64 stable trials per genotype;
- 256 paired fluctuating trials per genotype.

Let

`pC_stable`, `pR_stable`, `pC_fluct`, `pR_fluct`

be the corresponding reproductive-success proportions.

The material genotype-by-environment interaction is

`I = (pR_fluct - pC_fluct) - (pR_stable - pC_stable)`.

For each paired fluctuating trajectory, define a responsive win when responsive succeeds and constitutive fails, and a constitutive win for the opposite discordance. Ties do not enter the McNemar statistic.

The continuity-corrected signed statistic is

`z = sign(b-c) * max(0, |b-c|-1) / sqrt(b+c)`

where `b` is responsive-only success and `c` is constitutive-only success.

## Prospectively frozen acceptance gate

All criteria are required:

1. `pC_stable >= 0.95` and `pR_stable >= 0.95`, showing both architectures remain viable in the reference environment;
2. `abs(pR_stable - pC_stable) <= 0.03`, preventing a large baseline viability difference from masquerading as G×E;
3. `I >= +0.04`, requiring at least a four-percentage-point material environment-dependent reproductive advantage;
4. at least 8 discordant fluctuating pairs, preventing inference from a trivial number of events;
5. continuity-corrected paired McNemar `z >= 1.96` in the predicted responsive direction;
6. responsive-only reproductive successes are at least three times constitutive-only successes (`b >= 3 * max(1, c)`).

No seed count, seed value, phase schedule, horizon, circuit parameter, threshold, or statistic may be changed after output from `21001–21064` is observed.

## Interpretation boundary

A pass establishes the M5 exit-gate claim narrowly: with explicit stochastic expression, generic molecular regulation, finite proteome cost and ordinary metabolism/death/division, inherited regulatory architecture changes reproductive success depending on environmental structure, without any explicit fitness function.

It does not claim long-run fixation, optimality, spontaneous discovery of this exact circuit, or universal benefit of regulation. Those remain later evolutionary experiments.
