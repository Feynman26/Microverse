# M7-B generic membrane transport

M7-A made secondary metabolites physically present outside cells without granting any new membrane capability. M7-B adds the missing causal link: secondary molecules may cross the membrane only when the realized proteome contains a compatible transport protein.

## Sequence-derived recognition

Transport uses the same 16-bit protein sequence space as the rest of Microverse, but a distinct recognition domain. For extracellular metabolite `m`, the transport target is:

`transport_target(m) = ligand_signature(m) XOR 0x4190`

The transform is molecule-agnostic: the same operation is applied to every secondary extracellular metabolite. Affinity decays with Hamming distance and is zero beyond distance 4.

The domain mask was characterized against the M6 ancestor before integration. For each of the ten secondary metabolites (`C2`, `C3`, `W1`, `W2`, `CO2`, `AA`, `LIP`, `NUC`, `ROS`, `X`), the closest ancestral coding sequence is exactly distance 5 from its transport target. Therefore:

- the ancestor receives no secondary transport capability;
- one ordinary coding bit flip can move a protein from distance 5 to distance 4 and make the capability nonzero;
- a DNA mutation does not change transport immediately because activity is calculated from realized protein cohorts, preserving the M5 DNA -> expression -> protein causal chain.

This landscape permits promiscuity when one realized protein happens to lie near multiple molecular targets. Such pleiotropy is a consequence of sequence geometry, not a named multifunction flag.

## Reversible exchange

There is no authoritative `import`, `export`, `secretor`, or `cross_feeder` state. A compatible protein supplies membrane permeability/capacity. Direction is determined by the instantaneous concentration gradient:

- external concentration greater than internal concentration -> positive inward exchange;
- internal concentration greater than external concentration -> negative outward exchange;
- equal concentrations -> no net exchange.

The same protein can therefore import or export the same molecule as conditions change. This is the substrate needed for later evolved secretion/cross-feeding without programming either behavior directly.

## Capacity and cost

For a molecule with nonzero transport activity, the requested magnitude is bounded by:

- realized compatible protein abundance;
- cell volume;
- a saturating function of concentration-gradient magnitude;
- configured maximum transport rate;
- intracellular capacity for import;
- available intracellular amount for export.

Every unit actually moved consumes `secondary_transport_atp_cost_per_unit` ATP. ATP is converted to ADP using the existing adenylate accounting. If one cell proposes multiple simultaneous secondary exchanges but lacks ATP for all of them, every proposal is scaled by the same factor. No molecule is privileged by loop order.

Transport proteins also already carry the M5 opportunity costs of transcription, translation, finite ribosomal capacity, and finite proteome capacity.

## Simultaneous allocation and ecological fairness

Secondary transport is solved from a common pre-exchange snapshot.

1. Every living cell computes sequence-derived activities and signed exchange proposals from the same current intracellular/extracellular state.
2. Each cell's complete proposal vector is proportionally scaled to available ATP.
3. At each lattice site and for each metabolite, simultaneous positive import requests are summed.
4. If requests exceed the pre-exchange extracellular inventory, every requester is multiplied by the same availability ratio.
5. Aggregate imports are removed from the field and aggregate exports are then added.
6. Actual cell pool changes and exact ATP costs are applied.

Exported material is intentionally not available to another cell until the next transport phase. This prevents iteration order from creating zero-time producer-to-consumer transfer.

## Scope boundary

M7-B adds only generic membrane exchange. It does not add:

- a fitness or ecological-role API;
- public-good behavior;
- toxin behavior;
- cross-feeding labels;
- extracellular catalytic reactions;
- active pumping against a concentration gradient;
- regulatory intent or goal-directed secretion.

Active pumping can be considered later only if represented as ordinary molecular machinery with explicit energetic coupling rather than a behavioral command.

## Blocking validation

`tests/m7_membrane_transport_tests.gd` requires:

1. zero secondary transport activity in the ancestral proteome;
2. every secondary target exactly one coding bit outside the active ancestral radius;
3. coding mutation alone cannot change transport before compatible protein is realized;
4. the same transport mechanism supports import and export solely by gradient reversal;
5. molecule conservation across cell + extracellular field;
6. ATP expenditure with ATP+ADP conservation;
7. proportional allocation under scarce shared extracellular substrate;
8. proportional ATP limitation across simultaneous molecules;
9. exact deterministic replay from the same state/seed.

Only after this gate passes should M7-C test controlled producer/consumer chemistry such as W1/W2 cross-feeding and determine whether the existing metabolic network is sufficient or requires additional generic extracellular catalysis.
