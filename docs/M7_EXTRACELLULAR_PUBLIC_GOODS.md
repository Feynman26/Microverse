# M7-E — extracellular protein catalysis and public-resource opportunity

M7-E adds the first generic extracellular catalytic machinery. The implementation is deliberately molecular: no cell stores a producer, cooperator, cheater, public-good or beneficiary state.

## Physical mechanism

A realized protein cohort may leave the cell only if its 16-bit sequence contains the compressed `Dxxx` secretion/localization motif. Secretion:

- removes the exact protein molecules from the intracellular cohort;
- spends ATP per secreted protein unit;
- deposits those protein molecules into a sequence-specific extracellular field;
- lets those protein molecules diffuse slowly in the chamber.

DNA alone cannot secrete anything. A coding mutation changes extracellular capability only after expression produces the new sequence.

The ancestral proteome contains no `Dxxx` sequence. An ordinary one-bit mutation from ancestral locus 12 (`C136`) to `D136` creates a secretion-compatible protein. This same controlled protein is distance 4 from the first extracellular catalytic reaction signature, so the capability is accessible without gifting it to the ancestor.

## First extracellular reaction

`E01: 1 LIP -> 2 C2`

The modeled structural carbon is exactly conserved:

- `LIP` contains C4;
- each `C2` contains C2;
- one C4 unit therefore becomes two C2 units.

The catalyst is not consumed by the reaction. Flux is bounded by local substrate availability, sequence affinity, local extracellular protein abundance, kinetic saturation and elapsed time.

## Why this can become a public resource

C2 produced outside the cell has no owner. Any neighboring cell that physically encounters it and expresses compatible M7-B C2 transport protein may import it and pay the ordinary transport ATP cost.

A producer therefore pays both:

1. expression/proteome opportunity cost inherited from M5; and
2. explicit secretion ATP/protein cost in M7-E.

A nonproducer may capture the same extracellular product without paying the secretion cost. In a shared local site this creates a mechanistic exploitation opportunity. Spatial structure can change the result because product concentration is highest near the extracellular catalyst; there is no producer bonus or distance rule.

## Tick ordering

For each simulation tick:

1. existing extracellular chemicals and proteins diffuse;
2. basal and secondary membrane exchange runs;
3. secretion-compatible realized proteins are moved outside and ATP is debited;
4. extracellular protein catalysis transforms local chemical fields;
5. intracellular expression/metabolism, death/division and mechanics proceed.

This ordering means newly generated extracellular product cannot be imported at zero elapsed time. It becomes available to transport on the following tick after ordinary diffusion.

## Validation gate

`tests/m7_public_good_tests.gd` blocks M7-E on:

- zero ancestral protein secretion;
- one-bit accessibility of secretion plus extracellular catalytic activity;
- exact intracellular-protein to extracellular-protein transfer accounting;
- explicit ATP-to-ADP secretion cost;
- bounded/nonnegative E01 flux and exact LIP/C2 carbon stoichiometry;
- catalytic protein conservation during reaction;
- capture of producer-generated C2 by a nonproducer with matching transporter;
- equal shared-site product allocation with a measurable producer secretion cost;
- reduced distant nonproducer capture explained by local C2 concentration;
- exact same-state/seed replay including extracellular protein fields.

Passing this gate proves that costly extracellular digestion and exploitation by a nonproducer are physically possible. It does not claim that cooperation, cheating, stable coexistence or public-good production evolves spontaneously.

## Explicitly deferred

M7-E does not yet add extracellular protein degradation/recycling, detoxification, damaging extracellular chemistry, self-resistance, auxotrophic dependence, or growth-level population experiments on the released product. Those remain later M7 increments.
