# M5-C finite-proteome confirmation — prospectively frozen protocol

## Why a new mechanism is justified

The original M5-C circuit had explicit ATP/AA/NUC expression costs but no finite total proteome budget. A powered 24-seed test at 20-minute O2 phases produced a decisive null (`mean D_norm = -0.008394`, `|t| = 0.135113`). A subsequent independent 48-seed 10-minute confirmation also rejected the timescale hypothesis (`mean D_norm = +0.009202`, `t = +0.245552`, 25/48 negative), while all lineages reached normal division criteria.

Those results rule out further seed or timescale fishing. The next predeclared mechanism is a generic cellular proteome budget: protein expression must have an opportunity cost even when ATP and amino acids are abundant.

## Production mechanism

`SimConfig.proteome_capacity_reference_units = 5.0` defines one shared capacity of five `expression_reference_protein_count` units per cell. With the canonical reference count of 160, the default capacity is 800 realized protein units.

After ordinary stochastic transcription/translation and protein decay, the cell compares its realized proteome with that capacity. If it overshoots, every protein cohort is multiplied by the same scale factor. No sequence, locus, reaction, regulator or environmental condition receives priority.

Removed protein returns only its modeled amino-acid material to the `AA` pool. ATP already spent during translation is not refunded. Therefore unnecessary expression can waste energy and consume finite protein allocation, while down-regulation can free capacity for other proteins. This is a generic physical trade-off, not a fitness rule.

Initial cells are normalized to the same capacity without recycling because initialization constructs a state rather than executing a degradation event.

## Mechanism gates

Before any evolutionary-selection result is interpreted, CI must demonstrate that:

1. an unconstrained ancestral equilibrium exceeds the finite default budget;
2. initialization respects the exact shared capacity while preserving relative promoter-driven allocation;
3. allocation is independent of genome-array order;
4. overexpressing one locus displaces another protein rather than receiving free capacity;
5. runtime stochastic expression can hit the ceiling;
6. excess protein material is recycled without refunding expression ATP;
7. C/N/P accounting remains conserved;
8. all M0-M5-B regression gates remain green.

## Why prior seeds cannot close the new model

Seeds `1011–4044`, `5155–8488`, `13001–13006`, `14001–14024`, `15001–15004`, and `16001–16048` have already been observed under earlier no-budget models or discovery protocols. They are excluded from final inference after the production mechanism changes.

The first CI execution of the finite-proteome code may still execute an old confirmation test while the branch transitions. Any such output is characterization only and cannot close M5-C.

## Independent final assay — frozen before output

Final seeds: **`17001–17048`**, 48 completely new sequential RNG seeds.

Environment:

- stable: `O2 = 0.5`;
- fluctuating: `O2 = 6 -> 0 -> 6 -> 0 ...`;
- 50/50 duty cycle;
- half-period: **200 ticks = 20 biological minutes**;
- fluctuating starting offsets: `0, 100, 200, 300` ticks;
- glucose, nitrogen and phosphorus maintained as identical uniform reservoirs;
- mutation disabled.

The 20-minute schedule is restored as the canonical M5-C protocol because it was the original predeclared environment and equals the existing protein-turnover time constant. The earlier data show that timescale alone does not create a reproducible effect; any new result must arise after the independently justified finite-proteome mechanism.

For each lineage:

`r = ln(2) / T_division`.

For each seed:

`A_stable = r_responsive,stable - r_constitutive,stable`.

`A_fluctuating = mean_phase(r_responsive - r_constitutive)`.

`D_norm = (A_fluctuating - A_stable) / mean(r_responsive,stable, r_constitutive,stable)`.

The gate is deliberately two-sided. The finite proteome mechanism is expected to create a genuine trade-off, but the direction of selection is not prescribed.

All four conditions are required:

1. every stable and fluctuating lineage reaches ordinary division criteria before the existing 3600-tick timeout;
2. `abs(mean(D_norm)) >= 0.02`;
3. `abs(t) >= 2.01`, where `t = mean(D_norm) / [SD(D_norm)/sqrt(48)]`;
4. at least `31/48` seed-level values have the same sign as the observed mean.

The sign guard corresponds to a one-sided fair-sign probability of about 0.03 and prevents a few extreme seed outcomes from carrying the conclusion.

No seed, phase, offset, timeout, capacity, regulatory parameter, threshold, sample size or direction may be changed after `17001–17048` output is observed. If this assay fails, M5-C remains open and the finite-proteome mechanism itself is retained only if its independent mechanistic and regression gates remain scientifically useful; the failed selection result must not be converted into a pass.