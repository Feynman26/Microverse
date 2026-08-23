# M5-C amplified G×E confirmation — frozen before output

## Rationale

The prior final G×E panel using the low-amplitude validation circuit (`19001–19024`) was rejected: mean `D_norm = +0.009986`, SD `0.150358`, `t = 0.325354`, and 12/24 seeds were positive; at least one lineage also missed the 4800-tick endpoint. The result is retained as a falsification.

The mechanistic problem is signal amplitude, not a missing M5-B mechanism: regulation, ligand allostery, explicit expression, finite proteome competition, material/ATP costs, and RNG stream locality are independently validated. In the weak construct, R03 regulation changes too small a fraction of the realized proteome to produce a robust reproductive interaction.

This new construct therefore changes only the **shared basal allocation** to the regulated R03 protein and its catalytically inactive O2-compatible regulator. Both competitors receive the same promoter values and costs. They still differ in exactly one inherited field: the R03 promoter regulatory motif.

## Construct

Both competitors use promoter code `10000` for:

- locus 3: R03 oxidative-phosphorylation protein `0x369C`;
- locus 10: catalytically inactive O2-compatible regulator `0xCCCC`.

All other loci are unchanged from the prior M5-C compressed network. The constitutive R03 motif remains dormant; the responsive R03 motif exactly matches `0xCCCC`.

No fitness variable, reproduction modifier, genotype-specific proteome capacity, reaction bonus, or named oxygen-response routine is added.

## Environment

The environmental comparison is unchanged from the rejected G×E panel:

- reference: constant `O2 = 6`;
- fluctuating: `O2 = 6` for 40 biological minutes, then `O2 = 0` for 40 minutes, repeated;
- offsets: `0, 200, 400, 600` ticks.

The finite proteome capacity remains the production value `5 × 160 = 800` protein units.

## Hypothesis

Under constant high oxygen, ligand binding strongly inhibits the `0xCCCC` repressor, so responsive R03 expression approaches the constitutive state. During sustained anoxia, R03 cannot catalyse oxidative phosphorylation; the responsive architecture can reduce this now-unproductive high-abundance protein and reallocate a finite proteome while the constitutive architecture continues to express it.

Therefore the prospectively predicted interaction is positive:

`D_norm = (A_high/anoxic - A_high) / mean(r_R,high, r_C,high) > 0`.

## Independent confirmation panel

Seeds are fixed at **20001–20032** and have not been used by any prior M5-C panel.

The timeout is fixed at **7200 ticks (720 biological minutes)**. This changes only observation duration, not physiology, and is selected before output because the previous G×E panel showed that a 4800-tick endpoint can censor otherwise viable slow lineages.

All four gates are required:

1. every lineage reaches ordinary division criteria before timeout;
2. `mean(D_norm) >= +0.05`;
3. `t = mean / [SD/sqrt(32)] >= 2.04`;
4. at least `21/32` independent seed-level values are positive.

No seed, timeout, promoter, environment, phase, proteome capacity, regulatory parameter, threshold, direction, or sample size may be changed after output from `20001–20032` is observed.

A pass establishes the M5 exit-gate claim only: an inherited generic regulatory architecture changes relative reproductive performance when environmental structure changes, through explicit physiology and without explicit fitness. Long-run fixation and spontaneous evolution remain later milestones.
