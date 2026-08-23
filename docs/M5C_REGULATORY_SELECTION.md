# M5-C — Environment-dependent regulatory selection

## Purpose

M5-C closes the remaining M5 research gate: demonstrate that inherited regulatory architecture can experience a different **reproductive selection gradient** when environmental dynamics change, without assigning cells an explicit fitness value.

The assay uses ordinary Microverse physiology. Cells transcribe, translate, consume ATP/material, catalyse reactions, accumulate biomass, satisfy ordinary division criteria and can fail physiologically through the production engine. The experiment harness only maintains extracellular resources, imposes an oxygen schedule and measures outcomes.

M5-C deliberately separates two questions:

1. **causal reproductive effect** — measured with paired single-lineage generation rates, where population drift cannot hide a small physiological difference;
2. **small-population frequency dynamics** — retained as characterization, but not used as the final causal gate because N≈16 populations are dominated by stochastic lineage branching.

Large-population evolutionary inference belongs to M8.

## Competitors

The two competitors are isogenic except for one inherited promoter regulatory motif at locus 3, whose protein `0x369C` catalyses R03 oxidative phosphorylation.

- `constitutive`: R03 promoter motif is outside the regulator-binding radius.
- `responsive`: R03 promoter motif exactly matches assay regulator protein signature `0xCCCC`.

All basal promoter strengths, coding protein signatures, neutral markers, metabolic reactions and expression costs are otherwise identical.

The `0xCCCC` protein is deliberately outside every M4 catalytic activity radius, so it cannot create a hidden metabolic reaction. O2 is the only modeled ligand inside its M5-B allosteric radius.

`0xCCCC` has the high bit set, therefore promoter binding is repressive. Bit 14 is also set, so compatible ligand binding inhibits that regulatory contribution. Consequently:

- low intracellular O2 -> stronger R03 repression;
- high intracellular O2 -> O2 weakens the repressor -> R03 transcription moves toward basal expression.

No oxygen-response behavior is named or special-cased in the engine.

## Environmental control

Stable treatment:

`O2 = 0.5`

Fluctuating treatment:

`O2 = 6 -> 0 -> 6 -> 0 ...`

with a 50/50 square wave and `PHASE_TICKS = 200`.

At production `tick_dt_min = 0.10`, each phase lasts 20 biological minutes.

The schedules are matched for mean empty-cell membrane transport opportunity using production `oxygen_transport_km = 0.6`:

`0.5 / (0.6 + 0.5) = 0.5 * [6 / (0.6 + 6) + 0 / (0.6 + 0)] = 0.454545...`

This does not assert that realized uptake stays identical after physiology diverges. It removes the simpler design artifact that one treatment would have a different mean saturable transport opportunity by construction.

Glucose, nitrogen and phosphorus are maintained as identical uniform reservoirs. Mutation is disabled so the assay isolates selection on the two predefined inherited architectures.

## Molecular response timescale

M5-A currently uses:

- mRNA decay `0.25 min^-1` -> `tau_mRNA = 4 min`;
- protein decay `0.05 min^-1` -> `tau_protein = 20 min`.

The oxygen phase duration is therefore exactly one protein turnover time constant.

This creates a mechanistically important lag: ligand sensing and promoter occupancy can change quickly, while the realized catalytic protein pool cannot instantaneously follow. M5-C tests this directly. Deterministic mean-expression controls prove that one anoxic phase leaves less R03 protein in the responsive architecture and that R03 remains below the constitutive level during the early high-O2 opportunity even though sensing has already changed transcription.

## Population discovery 1 — original positive hypothesis rejected

The first population protocol and directional gate were committed before output was available.

Seeds: `1011, 2022, 3033, 4044`.

Original hypothesis:

`delta = log(R/C)_fluctuating - log(R/C)_stable > 0`

with at least 3/4 positive pairs and mean `delta > 0.05`.

The 3600-tick result rejected that hypothesis:

| Seed | Stable R/C | Stable log ratio | Fluctuating R/C | Fluctuating log ratio | Paired delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1011 | 9 / 11 | -0.191055 | 9 / 11 | -0.191055 | 0.000000 |
| 2022 | 10 / 9 | 0.100083 | 9 / 11 | -0.191055 | -0.291139 |
| 3033 | 8 / 9 | -0.111226 | 8 / 9 | -0.111226 | 0.000000 |
| 4044 | 9 / 8 | 0.111226 | 8 / 10 | -0.211309 | -0.322535 |

Mean stable log ratio `-0.022743`; mean fluctuating log ratio `-0.176161`; mean paired delta `-0.153418`; positive pairs `0/4`.

Both treatments reproduced and reached generation 3. The failed positive gate is preserved as a falsification, not rewritten into a pass.

## Population discovery 2 — negative-direction confirmation also rejected

The first discovery suggested that the 20-minute protein timescale might make the responsive architecture too slow for 20-minute O2 phases. That specific mechanism was independently supported by the R03 protein-lag controls.

A new population confirmation was therefore prospectively declared with unseen seeds `5155, 6266, 7377, 8488`, fixed 12-division endpoints, and a negative expected direction.

All eight stable/fluctuating populations reached the 12-division endpoint without extinction or cap artifacts, but genotype frequencies were:

| Seed | Stable R/C | Fluctuating R/C | Paired delta |
| --- | ---: | ---: | ---: |
| 5155 | 8 / 8 | 8 / 8 | 0.000000 |
| 6266 | 8 / 8 | 8 / 8 | 0.000000 |
| 7377 | 8 / 8 | 9 / 7 | +0.236389 |
| 8488 | 8 / 8 | 8 / 8 | 0.000000 |

Mean paired delta was `+0.059097`; strict negative pairs `0/4`. The negative-direction confirmation therefore failed.

This second failure is methodologically useful: the molecular phenotype is measurable, yet descendant-count frequency at N≈16 is quantized and dominated by stochastic lineage branching. Choosing more favorable seeds or weakening the directional criterion would be invalid.

## Final M5-C causal gate — paired lineage reproductive rate

The final gate was specified **before its six lineage seeds were executed**.

It uses the same production `SimulationEngine`, `CellState`, stochastic expression, metabolic solver, transport and viability rules. The only harness-specific operation is `max_cells = 1`, which prevents the simulation engine from replacing a mature founder with daughters. The harness stops at the first tick when that founder satisfies ordinary `ready_to_divide` criteria.

For a dividing lineage:

`r = ln(2) / T_division`

This is an analysis metric only; it never enters physiology.

For each seed:

`A_stable = r_responsive,stable - r_constitutive,stable`

The fluctuating advantage is averaged over four equally spaced starting phases `0, 100, 200, 300` ticks:

`A_fluctuating = mean_phase(r_responsive - r_constitutive)`

The environment-selection differential is:

`D = A_fluctuating - A_stable`

and is normalized by the seed's mean stable growth rate:

`D_norm = D / mean(r_responsive,stable, r_constitutive,stable)`

### Frozen seeds and acceptance criteria

Seeds: `13001, 13002, 13003, 13004, 13005, 13006`.

These seeds have not been used by either prior population experiment.

Four fluctuating phase offsets per seed prevent the result from depending on always starting a founder in the same O2 half-cycle.

The gate is deliberately **two-sided**. M5-C does not require regulation to be beneficial; it requires environmental dynamics to change the reproductive selection gradient.

Acceptance requires all of the following:

1. every constitutive and responsive lineage reaches ordinary division criteria before the 3600-tick timeout;
2. at least 5/6 seed-level normalized differentials have the same sign as the overall mean differential;
3. `abs(mean(D_norm)) >= 0.02`, i.e. the temporal environment shifts relative reproductive selection by at least 2% of baseline growth rate.

No new seeds, sign changes or threshold changes are permitted after this panel is observed. If this gate fails, M5-C remains unclosed and the result is recorded as another falsification.

## Interpretation boundary

A passing lineage gate would establish that:

- the two inherited regulatory architectures differ only through the intended molecular regulatory edge;
- stable versus fluctuating O2 changes their relative reproductive rate through normal expression/protein/metabolic dynamics;
- the effect does not require an explicit fitness variable;
- regulation can therefore be selected differently under different environmental dynamics at the causal lineage level.

It would **not** establish fixation probability in a finite population, spontaneous evolution of this particular circuit, or a universal benefit/cost of regulation. Those require larger replicated population experiments and belong to M8.

The two failed small-population gates remain valuable evidence that genetic drift and discrete lineage branching can obscure a real physiological selection gradient in tiny populations — itself an emergent property that Microverse should preserve rather than suppress.
