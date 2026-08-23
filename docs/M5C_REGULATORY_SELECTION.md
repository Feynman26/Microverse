# M5-C — Environment-dependent regulatory selection

## Purpose

M5-C closes the remaining M5 research gate: demonstrate that inherited regulatory architecture can experience a different reproductive selection gradient when environmental dynamics change, without assigning cells an explicit fitness value.

The assay uses ordinary Microverse physiology. Cells transcribe, translate, spend ATP/material, catalyse reactions, accumulate biomass and satisfy ordinary division criteria through the production engine. The experimental harness only maintains extracellular resources, imposes an oxygen schedule and measures outcomes.

Large-population fixation and long-run evolutionary inference remain M8 scope. M5-C uses causal first-division lineage rates so finite-population drift cannot be mistaken for or hide the molecular effect under test.

## Competitors and causal isolation

The two mutation-free competitors are isogenic except for one inherited promoter regulatory motif at locus 3, whose protein `0x369C` catalyses R03 oxidative phosphorylation.

- `constitutive`: R03 promoter motif is outside the regulator-binding radius.
- `responsive`: R03 promoter motif exactly matches regulator protein `0xCCCC`.

All basal promoter strengths, coding protein signatures, neutral markers, metabolic reactions and expression costs are otherwise identical.

`0xCCCC` is outside every M4 catalytic radius, so it cannot create a hidden metabolic reaction. O2 is the only modeled ligand inside its M5-B allosteric radius. Its sequence makes it a ligand-inhibited repressor in the generic M5-B grammar:

- low intracellular O2 -> stronger R03 repression;
- high intracellular O2 -> O2 weakens repression -> R03 transcription approaches basal expression.

No named oxygen-response or fitness behavior exists in the engine.

## Environmental control

The reference stable treatment is `O2 = 0.5`.

The fluctuating treatment alternates `O2 = 6` and `O2 = 0` with a 50/50 duty cycle. Glucose, nitrogen and phosphorus remain identical uniform reservoirs.

With production `oxygen_transport_km = 0.6`, stable and fluctuating treatments have the same time-averaged empty-cell saturable oxygen transport opportunity:

`0.5/(0.6+0.5) = 0.5 * [6/(0.6+6) + 0/(0.6+0)] = 0.454545...`

This does not force realized uptake to remain identical after physiology diverges. It removes the simpler artifact of unequal mean membrane transport opportunity by construction.

## Molecular response timescale

Production M5 expression uses:

- mRNA decay `0.25 min^-1` -> characteristic time `4 min`;
- protein decay `0.05 min^-1` -> characteristic time `20 min`.

Mechanistic controls independently establish that the responsive R03 promoter changes transcription with O2, one anoxic interval lowers responsive R03 protein relative to constitutive, and protein abundance remains lagged after O2 returns. Thus ligand sensing is immediate at the regulatory layer while the catalytic phenotype has finite molecular memory through ordinary protein turnover.

## Analysis metric

For a single lineage:

`r = ln(2) / T_division`.

For each RNG seed:

`A_stable = r_responsive,stable - r_constitutive,stable`.

Fluctuating advantage is averaged over four equally spaced starting phases of the square wave:

`A_fluctuating = mean_phase(r_responsive - r_constitutive)`.

The dimensionless environment-selection differential is:

`D_norm = (A_fluctuating - A_stable) / mean(r_responsive,stable, r_constitutive,stable)`.

`max_cells=1` in the lineage harness only prevents daughter replacement after maturity. The founder still runs the production `SimulationEngine`, `CellState`, transport, stochastic expression and metabolism and must satisfy ordinary `ready_to_divide()` criteria.

## Falsification history

### Population discovery — positive hypothesis rejected

Seeds `1011, 2022, 3033, 4044` tested a prospectively positive fluctuating-minus-stable descendant-frequency hypothesis. Paired log-ratio deltas were `0.000000, -0.291139, 0.000000, -0.322535`; mean `-0.153418`, with `0/4` positive. Both treatments reproduced and reached generation 3. The original hypothesis was rejected rather than rewritten.

### Population confirmation — negative hypothesis also rejected

Fresh seeds `5155, 6266, 7377, 8488` tested the subsequent protein-lag hypothesis at a fixed 12-division endpoint. All eight populations reached the endpoint without cap/extinction artifacts, but paired deltas were `0, 0, +0.236389, 0`; mean `+0.059097`. This rejected the negative population hypothesis and showed that N≈16 genotype counts are too quantized and drift-sensitive to serve as the causal M5-C gate.

### Six-seed lineage pilot — effect estimate unstable

Fresh seeds `13001–13006`, with 20-minute O2 phases and four phase offsets, produced:

`D_norm = [-0.029572, -0.082668, +0.163598, +0.408477, -0.055432, +0.479383]`.

Mean `+0.147298`, sample SD `0.246588`; all trajectories divided. The material-effect threshold passed but only 3/6 signs agreed, so the panel was not accepted. It was used only to plan a powered independent panel.

### Powered 20-minute lineage panel — decisive null

A prospectively frozen panel of 24 completely new seeds `14001–14024` tested the same 20-minute half-period. Acceptance required all trajectories to divide, `abs(mean D_norm) >= 0.02`, and `abs(t) >= 2.07` for 23 degrees of freedom.

Observed:

- all trajectories reached ordinary division criteria;
- mean `D_norm = -0.008394`;
- sample SD `0.304364`;
- SE `0.062128`;
- `abs(t) = 0.135113`.

The 20-minute regime therefore failed both the material-effect and statistical gates. This is retained as a genuine null result. No additional seeds or relaxed thresholds are added to that hypothesis.

## Timescale discovery after the powered null

Because the 20-minute environmental half-period equals one protein turnover time constant, a new isolated research spike asked whether the same molecular circuit behaves differently when the environment changes faster or slower. Production biology was not changed.

Discovery seeds `15001–15004` were frozen before output. Candidate half-periods were 10, 20, 40 and 80 biological minutes. A new timescale was eligible for independent confirmation only when `abs(mean D_norm) >= 0.10` and at least 3/4 seeds shared the mean sign; if several qualified, the shortest qualifying new timescale was predeclared as the nominee.

Observed discovery results:

| Half-period | Mean D_norm | Sample SD | Sign agreement |
| --- | ---: | ---: | ---: |
| 10 min | -0.142701 | 0.334143 | 3/4 negative |
| 20 min | -0.150933 | 0.322529 | 3/4 negative |
| 40 min | -0.170857 | 0.401908 | 3/4 negative |
| 80 min | -0.166798 | 0.403385 | 3/4 negative |

At least one trajectory somewhere in the full multi-timescale sweep did not divide before the 3600-tick timeout, so the sweep itself is not a validation result. The predeclared nomination rule nevertheless identifies **10 minutes** as the shortest new qualifying timescale. The apparently negative 20-minute four-seed calibration is explicitly not accepted because the independent powered n=24 experiment already established that 20-minute mean as null; this demonstrates why the four-seed sweep is discovery only.

## Prospectively frozen 10-minute confirmation

The final confirmation uses only the nominated 10-minute half-period and **48 completely new RNG seeds `16001–16048`**. No discovery, pilot, population or powered-null seed is reused.

For 10-minute phases the four starting offsets are `0, 50, 100, 150` ticks. Stable O2 remains `0.5`, fluctuating O2 remains `6/0`, and every other molecular/environmental parameter remains unchanged.

The discovery estimate at 10 minutes was mean `-0.142701`, SD `0.334143`, standardized effect `|d| ≈ 0.427`. A two-sided one-sample t power calculation gives about 82.6% power at n=48 and alpha≈0.05. This sample size was fixed before observing any `16001–16048` output.

The confirmation is directional because the discovery nominated a negative differential. **All** of the following must pass:

1. every stable and fluctuating lineage trajectory reaches ordinary division criteria before the existing 3600-tick timeout;
2. `mean(D_norm) <= -0.02`, preserving the discovered negative direction with at least 2% material magnitude;
3. `t = mean(D_norm) / [SD(D_norm)/sqrt(48)] <= -2.01`, the approximate two-sided 5% critical magnitude for 47 degrees of freedom while preserving the negative sign;
4. at least `31/48` independent seed-level `D_norm` values are negative (a complementary sign-consistency guard so a small number of extreme values cannot alone determine the conclusion).

No seed, phase, offset, timeout, direction, threshold or sample size may be changed after this panel is observed. If it fails, timescale alone is considered insufficient and M5-C remains open; the next justified mechanism is a genuine expression/proteome resource trade-off rather than further timescale fishing.

## Interpretation boundary

A passing 10-minute confirmation would establish that changing only temporal environmental structure changes the relative reproductive rate of otherwise matched inherited regulatory architectures through ordinary expression, protein turnover and metabolism, with no explicit fitness variable. It would not establish fixation probability, spontaneous evolution of this exact circuit or universal benefit/cost of regulation; those remain later open-evolution questions.
