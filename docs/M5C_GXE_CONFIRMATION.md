# M5-C final G×E confirmation — frozen before output

## Why this assay exists

M5-C must demonstrate that alternative inherited regulatory architectures have environment-dependent reproductive outcomes without an explicit fitness function.

Earlier experiments deliberately used a stronger nuisance control: stable `O2=0.5` was compared with a `6/0` square wave chosen to have the same time-averaged empty-cell saturable transport opportunity. That design generated several genuine null/falsification results. After adding a finite shared proteome and correcting an RNG stream-coupling artifact, two independent 48-seed panels still showed only a weak interaction under that mean-matched environment. The fresh `18001–18048` panel produced `mean D_norm = +0.020592`, `SD = 0.089518`, `t = 1.593683`, with 27/48 positive; all lineages divided. The material-effect threshold passed, but significance and sign-consistency gates failed.

That result is not converted into success. Instead, the final experiment asks the biologically direct genotype-by-environment question implied by the circuit.

## Production biology remains fixed

No M5 production parameter is changed for this final assay.

The already validated model contains:

- explicit stochastic transcription/translation and molecular turnover;
- catalytic flux driven only by realized protein cohorts;
- generic promoter binding and generic ligand allostery;
- finite proteome capacity (`5 × 160 = 800` protein units) with proportional, sequence-neutral allocation;
- ATP/material costs for expression;
- stream-stable one-uniform Poisson draws so a local lambda change cannot phase-shift unrelated random events;
- no explicit fitness variable.

Competitors remain isogenic except for the inherited regulatory motif controlling R03 oxidative phosphorylation.

## Environmental contrast

Reference environment:

- constant `O2 = 6`.

Fluctuating environment:

- `O2 = 6` for 40 biological minutes;
- `O2 = 0` for 40 biological minutes;
- repeated 50/50 square wave.

The four fluctuating phase offsets are `0, 200, 400, 600` ticks.

The 40-minute half-period was fixed before final seeds were observed. Production protein decay is `0.05 min^-1`, giving a characteristic lifetime of 20 minutes; 40 minutes therefore provides approximately two protein lifetime constants for the regulatory phenotype to change before the environment reverses.

This assay intentionally does **not** mean-match oxygen availability. That earlier constraint was useful to falsify trivial explanations but is not part of the M5 exit gate. Here the causal question is whether the same regulatory genotype changes its relative reproductive performance when the environment alternates between a state in which R03 is useful and one in which R03 cannot function.

Mechanistic directional hypothesis:

1. at `O2=6`, ligand binding inhibits the `0xCCCC` repressor, so the responsive R03 promoter approaches its constitutive state;
2. during sustained anoxia, the repressor becomes active and lowers responsive R03 expression;
3. because total proteome is finite, reducing an unusable oxidative-phosphorylation protein can free capacity and expression expenditure for the remainder of physiology;
4. therefore the responsive-minus-constitutive growth advantage is prospectively predicted to be larger in the high/anoxic environment than in constant high oxygen.

No code rewards that outcome; ordinary ATP, protein, metabolic flux, biomass and division determine it.

## Independent seeds and metric

Final seeds are **19001–19024**. None were used in any prior M5-C discovery, null, diagnostic or confirmation panel.

For each first-division lineage:

`r = ln(2) / T_division`.

For each seed:

`A_high = r_responsive,high - r_constitutive,high`.

`A_high/anoxic = mean_phase(r_responsive - r_constitutive)`.

`D_norm = (A_high/anoxic - A_high) / mean(r_responsive,high, r_constitutive,high)`.

The metric is analysis-only and never enters physiology.

## Prospectively frozen acceptance gate

All four criteria are required:

1. every stable-high and high/anoxic lineage reaches ordinary division criteria before the existing assay timeout;
2. `mean(D_norm) >= +0.02`;
3. `t = mean(D_norm) / [SD(D_norm)/sqrt(24)] >= 2.07`;
4. at least `17/24` independent seed-level `D_norm` values are positive.

The t criterion deliberately uses the approximate two-sided 5% critical magnitude for 23 degrees of freedom despite the directional hypothesis. The sign guard prevents a small number of extreme values from carrying the conclusion.

No seed, phase length, oxygen level, timeout, proteome capacity, regulatory parameter, threshold or sample size may be changed after output from `19001–19024` is observed.

## Interpretation boundary

A pass establishes the M5 exit-gate claim: inherited generic regulatory architecture can change relative reproductive outcome when environmental structure changes, through explicit molecular physiology and without an explicit fitness function.

It does not establish long-run fixation, spontaneous evolution of this exact circuit, optimality of regulation, or universal benefit of sensing. Those belong to later open-evolution experiments.
