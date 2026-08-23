# M5-C ROS-responsive confirmation

## Why the R03-targeted circuit was retired

Several independent panels showed that regulating oxidative phosphorylation (R03) changed molecular state but did not create a sufficiently reproducible reproductive interaction under the frozen M5-C thresholds. The strongest early signal did not survive independent confirmation. Adding a finite proteome and then a shared translation-capacity constraint improved the biological model but the independent ribosome-limited panel `22001..22064` still produced only a +1.5625 percentage-point interaction with corrected McNemar z=0.588348 (15 responsive-only versus 11 constitutive-only outcomes). That result remains a falsification and is not reinterpreted as success.

The next validation circuit therefore changes the biological target, not the statistical gate. The production chemistry already contains a causal oxidative-stress path:

`O2 -> R03 flux -> ROS -> damage`

and an ATP-consuming control reaction:

`R09: ROS + ATP -> ADP`.

M5-B already supplies generic O2-conditioned regulation through molecular affinity. The final circuit asks whether an alternative inherited promoter architecture can condition R09 abundance on O2 and thereby alter reproductive outcome without any fitness primitive.

## Circuit

Both competitors have identical loci, coding sequences, basal promoters, neutral markers and regulator abundance. R03, R09 and the O2-compatible regulator use promoter code 10000 in both genomes. The sole inherited difference is the R09 regulatory motif:

- constitutive: dormant motif, no compatible promoter edge;
- responsive: motif matches the O2-compatible repressor.

The regulator is catalytically inactive across the M4 reaction catalog. O2 allosterically inhibits its repressive contribution, so responsive R09 expression is more strongly repressed when O2 is absent and is relieved when O2 is high. No code checks for a named stress response.

## Environment

The assay preserves the previously frozen genotype-by-environment schedule:

- stable-high: O2 = 6 continuously;
- fluctuating: O2 alternates 6 -> 0 every 400 ticks (40 biological minutes);
- glucose, nitrogen and phosphorus are maintained at their existing M5-C reservoir values;
- mutations are disabled by the M5-C config;
- success is reaching the ordinary `ready_to_divide()` state before death;
- death before first division is a genuine reproductive failure.

## Independent panel frozen before observation

Seeds: `23001..23064` inclusive. No seed from an earlier M5-C population, lineage, timescale, GxE, proteome, ribosome or reproductive-success panel is reused.

Each seed produces two stable-high runs and four phase-paired fluctuating comparisons per genotype. The phase offsets are 0, 200, 400 and 600 ticks.

## Acceptance criteria frozen before observation

All six criteria must pass simultaneously:

1. constitutive and responsive stable-high reproductive success are both >= 95%;
2. the absolute stable-high success difference is <= 3 percentage points;
3. the responsive genotype-by-environment interaction is >= +4 percentage points;
4. at least 8 fluctuating pairs are discordant;
5. the continuity-corrected signed McNemar statistic is >= 1.96;
6. responsive-only successes are at least 3x constitutive-only successes.

The thresholds are intentionally unchanged from the failed R03-targeted reproductive-success panel. They are not adjusted after observing `23001..23064`.

## Prior mechanistic gate

Before this confirmation was enabled in CI, a separate ROS-circuit mechanics suite had to pass. It verifies that:

- competitors differ in exactly one inherited field;
- the O2-compatible regulator is not a hidden catalyst;
- production R03 consumes O2 and produces ROS;
- the target protein is the exact production R09 ROS-control catalyst;
- high intracellular O2 relieves repression of responsive R09;
- constitutive R09 has no corresponding regulatory edge.

Only after that gate passed was this confirmation panel frozen and enabled.
