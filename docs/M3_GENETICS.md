# M3 Genetics: inheritance before adaptation

## Purpose

M3 introduces explicit heritable molecular information while intentionally keeping it phenotypically neutral. This separation is deliberate: inheritance, mutation, ancestry and reproducibility must be validated before genotype is allowed to alter metabolism or regulation.

M4 maps protein signatures onto catalytic function. M5 maps promoter/regulatory information onto expression. M10 later makes mutation rate itself evolvable through replication/repair machinery.

## 1. Discrete genotype representation

Each `Gene` contains four discrete fields:

- `locus_id`: stable locus identity used for ancestry and future duplication tracking;
- `promoter_code`: integer 0-10000, later interpreted as expression potential;
- `protein_signature`: 16-bit molecular signature, later interpreted through protein/reaction affinity;
- `neutral_marker`: inherited marker guaranteed not to alter molecular function, retained as a control for drift and ancestry.

Discrete storage avoids float-serialization ambiguity in exact genotype identity and replay.

## 2. Ancestor

The M3 ancestor has twelve loci. Their promoter codes and signatures are deliberately abstract in this milestone. The values are not named `glucose_gene`, `repair_gene`, etc. because that would prematurely bake physiological semantics into genotype.

M4 may replace/recalibrate the ancestral signatures after characterizing the digital catalytic landscape, but the representation and ancestry semantics should remain stable.

## 3. Exact identity

A genome exposes:

- `canonical_key()`: ordered full-locus representation used to resolve exact identity;
- `fingerprint()`: deterministic compact rolling identity used for indexing/telemetry;
- `exact_equals()`: canonical exact equality;
- `checksum()`: state-reproducibility contribution.

The compact fingerprint is not cryptographic and must never be treated as collision-proof. Exact genotype comparisons use the canonical representation.

## 4. Inheritance

At division:

1. the parent's genome is copied deeply for each daughter;
2. no daughter, sibling, or ancestor shares mutable `Gene` objects;
3. each daughter is passed independently through the mutation engine;
4. the mutated genotype becomes that daughter's birth genotype;
5. mutation events are recorded against the daughter and parent identities.

This prevents a mutation in one descendant from retroactively altering a sibling or ancestor.

## 5. M3 mutation operators

### Promoter-code mutation

Per locus, with configured probability:

- sample a bounded integer step;
- add it to promoter code;
- clamp to the valid 0-10000 range;
- guarantee that a triggered mutation changes the code.

M3 consequence: none on physiology.

M5 consequence: expression dynamics.

### Protein-signature mutation

Per locus, with configured probability:

- choose one of the active signature bits;
- flip that bit;
- retain the 16-bit mask.

M3 consequence: none on physiology.

M4 consequence: altered catalytic/binding affinities and possible pleiotropic side activity.

### Neutral-marker mutation

Per locus, with configured probability:

- replace the marker with a new deterministic-RNG sampled integer;
- guarantee a changed value.

This field remains nonfunctional in future milestones and supplies an explicit neutral control class.

## 6. Mutation rate semantics

M3 rates are external simulation parameters and apply independently per gene per replication event:

- promoter mutation probability;
- signature mutation probability;
- neutral-marker mutation probability.

These are not biological claims about real mutation rates. Their role is to validate the evolutionary infrastructure.

M10 replaces the external fixed-rate assumption with a mechanistic baseline error rate modified by costly replication/repair machinery.

## 7. Event semantics

Every mutation event contains at least:

- stable `mutation_id`;
- mutation type;
- locus ID;
- old value;
- new value;
- daughter cell ID;
- parent cell ID;
- generation;
- parent genotype fingerprint;
- resulting daughter genotype fingerprint;
- simulation tick/time supplied by the generic event logger.

Signature mutations additionally record the flipped bit index.

Mutation events are semantic history, not a compressed state-reconstruction substitute. Exact future persistence will use snapshots plus event/intervention history.

## 8. Randomness

All mutation decisions and values use the simulation-owned `DeterministicRng`.

Required invariant:

`same model + same initial state + same RNG state + same interventions -> same mutations`

No mutation code may call a global/unseeded random source.

## 9. Neutrality gate

M3 is intentionally unable to evolve adaptive physiology from promoter/signature changes. This is a feature, not a missing implementation.

The milestone answers:

- Is molecular information inherited correctly?
- Can it mutate without aliasing ancestors/siblings?
- Are exact mutation histories reproducible?
- Can neutral diversity accumulate?
- Do realized mutation counts match configured stochastic expectations?

Only after those answers are validated does M4 give coding signatures functional meaning.

## 10. Bridge to M4

M4 should not replace `protein_signature` with a direct `enzyme_type` enum. Instead it introduces:

- reaction signatures/motifs;
- a validated protein-signature-to-reaction affinity function;
- catalytic promiscuity;
- explicit protein abundance/cost;
- competing metabolic routes.

A signature mutation will then have downstream consequences by changing affinities, not because the mutation is labelled beneficial/deleterious.

## 11. Bridge to M5

M5 interprets promoter/regulatory information through explicit expression dynamics:

`promoter + regulators -> transcription -> mRNA -> translation -> protein abundance`

Promoter mutations therefore acquire context-dependent phenotypic effects only after a costly expression system exists.

## 12. Structural evolution is intentionally deferred

M3 does not add duplication, deletion or rearrangement. Point-like inheritance must be stable first. Gene duplication/deletion enters after M4 can assign functional consequences and after ancestry semantics for duplicated loci are specified.

When duplication arrives, `locus_id` ancestry must distinguish:

- ancestral locus identity;
- duplicated copy identity;
- origin event;
- later divergence.

## 13. Validation suite

M3 tests cover:

- stable 12-locus ancestor;
- deep-copy exact equality;
- mutation-disabled exact inheritance;
- independent sister genomes;
- no parental mutation through aliasing;
- forced promoter/signature/neutral mutation;
- same-seed exact mutation sequence;
- different-seed divergent mutation history;
- neutral marker producing identical physiology;
- realized neutral mutation frequency within a statistical tolerance;
- mutation events containing ancestry/genotype context;
- complete same-seed simulation event-history identity.
