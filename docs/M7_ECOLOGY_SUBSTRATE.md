# M7 ecological innovation substrate

M7 begins only after M6 established that physical position is causally relevant. The goal is now to make ordinary extracellular chemistry rich enough that relationships such as cross-feeding, exploitation, detoxification and antagonism can become consequences of molecular fluxes rather than named behaviors.

This document describes **M7-A**, the extracellular molecular substrate and lysis-accounting layer. It deliberately stops before evolvable secondary transport, secretion or extracellular catalysis.

## M7-A molecular world

The intracellular chemistry already contains nineteen explicit metabolite identities. M7-A gives fourteen of those molecules a physical extracellular representation:

| Intracellular ID | Extracellular field | Initial secondary amount | Role at this layer |
| --- | --- | ---: | --- |
| `G` | `glucose` | existing basal value | basal imported resource |
| `O2` | `oxygen` | existing basal value | basal imported resource |
| `NH4` | `nitrogen` | existing basal value | basal imported resource |
| `P` | `phosphorus` | existing basal value | basal imported resource |
| `C2` | `carbon_c2` | 0 | secondary carbon substrate |
| `C3` | `carbon_c3` | 0 | secondary carbon substrate |
| `W1` | `waste_1` | 0 | diffusible metabolic intermediate |
| `W2` | `waste_2` | 0 | diffusible metabolic intermediate |
| `CO2` | `carbon_dioxide` | 0 | oxidized carbon product |
| `AA` | `amino_acids` | 0 | reusable precursor |
| `LIP` | `lipids` | 0 | reusable precursor |
| `NUC` | `nucleotides` | 0 | reusable precursor |
| `ROS` | `oxidant` | 0 | ordinary damaging chemical substrate for later interaction mechanics |
| `X` | `neutral_x` | 0 | semantically neutral diffusible molecule |

`BIO`, ATP/ADP and NAD/NADH remain non-diffusible internal state. Field existence is a physical property only: it does not imply that a cell can import, export, sense or benefit from that molecule.

The ten secondary fields start at zero by default. This preserves the M0-M6 ancestral environment while allowing metabolism, lysis and later secretion to populate a richer extracellular chemical state.

## Diffusion and numerical semantics

Every extracellular field uses the existing finite-difference diffusion solver with reflecting chamber boundaries. Each secondary diffusion coefficient is configuration data and is checked against the same explicit-solver stability condition as the basal resources.

The world registers the four historical fields first and appends secondary fields in deterministic molecular-ID order. This makes the expanded world checksum stable for a fixed model/configuration and preserves exact replay semantics.

## Lysis material accounting

Before M7, death returned only the four basal imported resources. That discarded metabolically meaningful intermediates, precursors and structural biomass from the ecological world.

M7-A changes death into a complete release of every intracellular metabolite that has an extracellular representation. Structural `BIO` is hydrolyzed using the exact reverse material stoichiometry of biomass reaction R12:

`1 BIO -> 2 AA + 1 LIP + 2 NUC`

Given the digital elemental definitions, both sides contain exactly C12/N4/P2. This is a bookkeeping decomposition, not a special corpse-resource or decomposition behavior.

Material physically stored in expression state is also returned using the same accounting already used by M5:

- protein material returns to `AA` according to `translation_aa_cost_per_event`;
- mRNA material returns to `NUC` according to `transcription_nuc_cost_per_event`.

The death event records the exact positive `released_pools` ledger that was added to extracellular fields. This supports later counterfactual and interaction analysis without adding ecological labels to cell state.

## Explicit non-capabilities of M7-A

M7-A does **not** give the ancestor uptake of W1, W2, C2, C3, AA, LIP, NUC, ROS, X or CO2. It does not add secretion, public-good production, toxin behavior, cooperation, cheating, mutualism or cross-feeding APIs.

A test deliberately places W1 outside an ancestral cell and verifies that the field remains untouched and intracellular W1 remains zero under the existing membrane-transport phase. Thus merely registering a molecule in the chamber cannot create a biological function.

## Validation gate

`tests/m7_extracellular_tests.gd` blocks M7-A on:

1. one-to-one catalog-to-world registration for all fourteen extracellular metabolites;
2. zero default inventory for every secondary field;
3. conservative, nonnegative local diffusion of a secondary metabolite;
4. absence of accidental secondary uptake;
5. exact structural C/N/P conservation across lysis, including BIO, protein and mRNA storage;
6. equality between the death-event release ledger and material actually added to the world;
7. exact same-state/seed replay after lysis.

Passing this gate proves only that the ecological chemical substrate exists and is materially accountable.

## Next dependency: M7-B generic membrane transport

The next increment may add secondary import/export only through realized molecular machinery. The transport mechanism must be generic and sequence-derived, must debit intracellular material exactly, must preserve fair allocation under scarce external substrate and must impose protein/energetic opportunity cost where appropriate.

The design must be characterized before integration so the ancestor is not accidentally granted broad secondary transport. Controlled genotypes may be used to prove capability, but labels such as `secretor`, `cross_feeder`, `cooperator` or `toxin_producer` remain forbidden from authoritative cell state.

Only after generic transport passes its own conservation and replay gates should M7 proceed to controlled W1/W2 cross-feeding, extracellular catalysis, exploitation and damaging-chemistry experiments.
