# M5 — explicit expression, regulation, sensing, and phenotypic noise

## Purpose

M5 inserts a transient molecular phenotype between genotype and catalytic function. A gene no longer contributes catalytic activity merely because its promoter code is high. Living cells must transcribe mRNA, translate protein, retain or degrade those molecules, pay resource costs, regulate expression through molecular binding, sense chemistry through generic allostery, and inherit a stochastic fraction of molecular state at division.

The authoritative causal chain becomes:

`environment -> uptake -> intracellular ligands -> allostery -> regulation -> mRNA -> protein -> catalytic affinity -> reaction flux -> ATP/precursors/BIO -> volume -> division`

Two cells with identical genomes can therefore have different immediate physiology. Conversely, a genotype mutation changes future molecular dynamics without instantaneously rewriting the already existing proteome.

## State separation

### Genome

Inherited discrete information:

- stable `locus_id`;
- `promoter_code` controlling basal transcription propensity;
- `protein_signature` controlling catalytic, regulatory, and allosteric molecular identity;
- `regulatory_signature` controlling which proteins can bind the promoter;
- `neutral_marker` for drift/lineage controls.

### ExpressionState

Transient per-cell state:

- mRNA abundance per locus;
- protein abundance per locus.

Expression state is not genotype identity. It is partitioned at division and can differ between clonal sisters.

## Basal expression kinetics

Ignoring resource limitation and regulation:

`dm_i/dt = k_tx * promoter_i - k_mdeg * m_i`

`dp_i/dt = k_tl * m_i - k_pdeg * p_i`

Default rates make ancestral steady-state protein abundance approximately equal to the former M4 promoter-strength proxy. This preserves continuity while removing that proxy from authoritative living-cell physiology.

## Expression cost and material accounting

Transcription consumes ATP and NUC precursor. Translation consumes ATP and AA precursor. ATP becomes ADP. mRNA/protein degradation returns the represented NUC/AA fraction.

Turnover is solved from one molecular snapshot. Material released by degradation during an interval can be reused for synthesis in that same interval, but all release and demand are computed before any locus is updated. Therefore locus array order cannot confer an energetic advantage.

The expression layer exposes structural totals for material sequestered in mRNA/protein. Combined with metabolite structural totals, M5 can audit represented C/N/P across:

`free AA/NUC <-> mRNA/protein`.

This turns expression into an explicit opportunity cost rather than a hidden abundance penalty.

## Generic regulation

M5 contains no named transcription factors, stress genes, dormancy genes, quorum genes, or behavior flags.

Every expressed protein can bind a promoter if its 16-bit `protein_signature` lies within the configured Hamming radius of a target gene's inherited `regulatory_signature`.

`occupancy = effective_protein_abundance * exp(-decay * HammingDistance)`

Proteins with the high signature bit clear contribute activating occupancy; proteins with it set contribute repressing occupancy. This is a digital-chemistry convention, not a biological label.

For a target promoter:

`net = (activation - repression) / (1 + activation + repression)`

`regulation_factor = clamp(1 + gain * net, min_factor, max_factor)`

Transcription is proportional to:

`promoter_strength * regulation_factor`.

Because `protein_signature` also determines catalytic affinity, a coding mutation can alter catalysis, regulation, sensing, or several simultaneously. That pleiotropy is intentional.

## Generic chemical sensing by allostery

M5 does not add a `sensor` class or a `respond_to_starvation()` primitive.

Every digital metabolite has a fixed 16-bit `ligand_signature`. Any protein sufficiently close to a ligand in signature space can be allosterically modulated by the molecule's intracellular concentration.

For a protein-ligand pair inside the allosteric radius:

`affinity = exp(-allosteric_decay * distance)`

`scaled_ligand = concentration * affinity`

`ligand_occupancy = scaled_ligand / (Km + scaled_ligand)`

The strongest compatible ligand modulates the protein's regulatory effective abundance. Protein bit 14 determines whether binding potentiates or inhibits its regulatory activity.

This means the engine never says that glucose means food, ROS means stress, or X means signal. It exposes concentration and molecular compatibility only. A lineage gains an environment-conditioned response only if its evolved molecular wiring makes that correlation useful.

The initial ligand signatures are deliberately farther than the active allosteric radius from all ancestral proteins. Thus generic sensing is possible from generation zero as physics, but the ancestor is not handed a new environment-response circuit. Mutations can move protein signatures into ligand-binding neighborhoods.

## Ancestral regulatory topology

The ancestor contains a sparse ring of exact protein-to-promoter motif matches. It guarantees the regulatory grammar is exercised and measurable without encoding an adaptive strategy. The ring crosses the high-bit boundary, so some edges activate and some repress.

## Stochastic expression

Each transcription and translation request receives bounded multiplicative noise from the single universe RNG:

- same complete state + same RNG state -> same trajectory;
- different seeds can separate molecular trajectories of the same genotype;
- UI/render state never contributes biological randomness.

Direct isolated-cell compatibility controls that omit the universe RNG run noise-free. Authoritative `SimulationEngine` runs always pass the seeded RNG.

## Division and molecular carryover

Metabolite pools use the cell-volume partition ratio. mRNA and proteins receive seeded per-locus partition perturbations while conserving exact parent totals:

`daughter_A + daughter_B = parent`

for every tracked mRNA/protein abundance.

A daughter DNA mutation occurs after molecular partition. Therefore parental proteins can transiently persist in a daughter whose new genome encodes a changed protein. This carryover is deliberate and creates realistic phenotype lag instead of instantaneous genotype-to-phenotype replacement.

## M5 acceptance gates

M5 is accepted only if tests demonstrate:

1. initial expression preserves the M4 activity scale;
2. zero protein abundance means zero catalytic flux;
3. transcription/translation incur ATP and material commitments;
4. expression turnover conserves represented C/N/P when free and macromolecular pools are counted together;
5. generic activator and repressor binding move expression in opposite directions;
6. ordinary metabolite binding can change regulation and downstream mRNA without a named response API;
7. ancestral proteins do not accidentally receive the new ligand-response capability merely from milestone migration;
8. same-seed molecular noise is exactly reproducible;
9. molecular partition conserves mRNA/protein while allowing clonal sister differences;
10. regulatory motifs are inherited and mutable;
11. full chemistry + expression + sensing + genetics remains replayable from the same seed;
12. all M0-M4 gates remain green.

Only after these gates should M6 make finite cell geometry, crowding, pushing, local colony structure, and spatial ecology causally relevant.
