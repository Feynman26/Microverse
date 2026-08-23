# M5-C — Environment-dependent regulatory selection

## Purpose

M5-C closes the remaining M5 research gate: demonstrate that inherited regulatory architecture can experience a different **reproductive selection gradient** when environmental dynamics change, without assigning cells an explicit fitness value.

The assay uses ordinary Microverse physiology. Cells transcribe, translate, consume ATP/material, catalyse reactions, accumulate biomass, satisfy ordinary division criteria and can fail physiologically through the production engine. The harness only maintains extracellular resources, imposes an oxygen schedule and measures outcomes.

M5-C separates two scales:

1. **causal reproductive effect** — paired single-lineage generation rates;
2. **small-population frequency dynamics** — retained as characterization, not used as the final causal gate because tiny populations are strongly affected by stochastic lineage branching.

Large-population evolutionary inference belongs to M8.

## Competitors and causal isolation

The two competitors are isogenic except for one inherited promoter regulatory motif at locus 3, whose protein `0x369C` catalyses R03 oxidative phosphorylation.

- `constitutive`: R03 promoter motif is outside the regulator-binding radius.
- `responsive`: R03 promoter motif exactly matches regulator protein `0xCCCC`.

All basal promoter strengths, coding protein signatures, neutral markers, metabolic reactions and expression costs are otherwise identical.

`0xCCCC` is outside every M4 catalytic activity radius, so it cannot create a hidden metabolic reaction. O2 is the only modeled ligand inside its M5-B allosteric radius. Its sequence makes it a ligand-inhibited repressor in the generic M5-B grammar:

- low intracellular O2 -> stronger R03 repression;
- high intracellular O2 -> O2 weakens the repressor -> R03 transcription moves toward basal expression.

There is no named oxygen-response behavior in the engine.

## Environmental control

Stable treatment: `O2 = 0.5`.

Fluctuating treatment: `O2 = 6 -> 0 -> 6 -> 0 ...`, 50/50 square wave, `PHASE_TICKS = 200`.

At `tick_dt_min = 0.10`, each phase lasts 20 biological minutes.

The schedules are matched for mean empty-cell membrane transport opportunity using production `oxygen_transport_km = 0.6`:

`0.5 / (0.6 + 0.5) = 0.5 * [6 / (0.6 + 6) + 0 / (0.6 + 0)] = 0.454545...`

This does not force realized uptake to remain equal after physiology diverges; it removes the trivial artifact of unequal mean saturable transport opportunity by construction.

Glucose, nitrogen and phosphorus are maintained as identical uniform reservoirs. Mutation is disabled to isolate selection on the two inherited architectures.

## Molecular response timescale

M5-A uses:

- mRNA decay `0.25 min^-1` -> `tau_mRNA = 4 min`;
- protein decay `0.05 min^-1` -> `tau_protein = 20 min`.

The O2 phase duration is therefore exactly one protein turnover time constant. Mechanistic controls independently prove that one anoxic phase leaves less R03 protein in the responsive architecture and that R03 remains below constitutive during the early high-O2 opportunity even though ligand sensing has already changed transcription.

## Falsification history

### Population discovery 1 — positive directional hypothesis rejected

Seeds `1011, 2022, 3033, 4044` tested a predeclared expectation that fluctuation would favor the responsive architecture.

Observed paired `delta = log(R/C)_fluctuating - log(R/C)_stable` values were:

- `0.000000`
- `-0.291139`
- `0.000000`
- `-0.322535`

Mean delta `-0.153418`; positive pairs `0/4`. Both treatments reproduced and reached generation 3. The positive hypothesis was therefore rejected rather than rewritten.

### Population confirmation 2 — negative directional hypothesis also rejected

The protein-lag mechanism motivated a fresh, prospectively negative test with seeds `5155, 6266, 7377, 8488` at a fixed 12-division endpoint.

All eight stable/fluctuating populations reached the endpoint without cap or extinction artifacts, but R/C outcomes were:

- seed 5155: stable `8/8`, fluctuating `8/8`, delta `0.000000`;
- seed 6266: stable `8/8`, fluctuating `8/8`, delta `0.000000`;
- seed 7377: stable `8/8`, fluctuating `9/7`, delta `+0.236389`;
- seed 8488: stable `8/8`, fluctuating `8/8`, delta `0.000000`.

Mean delta `+0.059097`; strict negative pairs `0/4`. This rejected the negative population confirmation. It also showed that N≈16 descendant counts are too quantized/noisy to serve as the causal M5-C gate.

### Lineage pilot — material mean effect, insufficient n for sign consistency

To remove population branching drift, the next assay measured first-division growth rate using the production `SimulationEngine`, `CellState`, expression, transport and metabolism, with `max_cells=1` solely to prevent daughter replacement after ordinary division readiness is reached.

For a lineage:

`r = ln(2) / T_division`.

For each seed:

`A_stable = r_responsive,stable - r_constitutive,stable`

and fluctuating advantage is averaged across phase offsets `0, 100, 200, 300` ticks:

`A_fluctuating = mean_phase(r_responsive - r_constitutive)`.

The normalized environment-selection differential is:

`D_norm = (A_fluctuating - A_stable) / mean(r_responsive,stable, r_constitutive,stable)`.

Pilot seeds `13001–13006` were frozen before execution. All 60 lineage trajectories divided before timeout. Their `D_norm` values were:

- `-0.029572`
- `-0.082668`
- `+0.163598`
- `+0.408477`
- `-0.055432`
- `+0.479383`

Mean `+0.147298`; sample SD `0.246588`; standard error at n=6 about `0.100669`.

The predeclared material-effect criterion (`abs(mean) >= 0.02`) passed, but the predeclared requirement for 5/6 identical signs failed (3/6). The pilot is therefore not accepted as final validation. It is used only for variance/sample-size planning.

## Final powered lineage validation — frozen before output

The pilot implies an observed standardized effect of roughly `0.1473 / 0.2466 ≈ 0.60`. A conventional two-sided normal-approximation planning calculation for ~80% power at alpha ~0.05 gives an order of 22 independent seed-level observations. The final panel therefore uses **24 completely new sequential seeds**:

`14001–14024`.

The pilot seeds and all population seeds are excluded from final inference.

Each seed again averages fluctuating physiology across four equally spaced phase offsets `0, 100, 200, 300`, so the result cannot depend on the arbitrary starting phase of the square wave.

The final criterion is deliberately **two-sided**: M5-C requires that temporal environment changes relative reproductive selection on the regulatory architecture; it does not prescribe which architecture must win.

Acceptance requires all of the following:

1. every stable and fluctuating lineage trajectory reaches ordinary division criteria before the 3600-tick timeout;
2. `abs(mean(D_norm)) >= 0.02` — at least a 2% shift relative to baseline growth rate;
3. `abs(t) >= 2.07`, where `t = mean(D_norm) / [SD(D_norm)/sqrt(24)]`, the predeclared approximate two-sided 5% critical value for 23 degrees of freedom.

No sign, threshold, seed, phase-offset or sample-size changes are permitted after this panel is observed. If it fails, M5-C remains open.

## Interpretation boundary

A passing powered lineage gate establishes that stable versus fluctuating O2 changes the **mean relative reproductive rate** of two otherwise matched inherited regulatory architectures through normal expression/protein/metabolic dynamics, without an explicit fitness variable.

It does not establish fixation probability in a finite population, spontaneous evolution of this particular circuit, or a universal benefit/cost of regulation. Those require larger replicated population experiments and belong to M8.

The failed small-population gates are retained because they demonstrate an important emergent property rather than a defect: genetic drift and discrete lineage branching can obscure a causal physiological selection gradient in tiny populations.
