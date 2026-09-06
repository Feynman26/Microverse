# P3-B — Dense secondary membrane transport

## Purpose

P3-A reduced the median 1,000-cell tick on the AMD Ryzen 5 4500U from
3,022.54 ms to 2,131.76 ms without changing checksums or events. Reconstructing
the unaffected P0 phase timings places secondary transport at approximately
782 ms/tick, ahead of expression at 648 ms/tick and the now-dense metabolism at
approximately 430 ms/tick. Secondary transport is therefore the next measured
hot path.

The historical M7 implementation evaluates ten extracellular metabolites for
every cell. Each evaluation independently sorts and traverses the complete
protein-cohort state, then the allocator builds several nested Dictionaries for
proposals, activities, scarcity scales and realized exchange.

P3-B changes only that execution strategy. It does not change recognition,
gradients, ATP costs, proportional scarcity, world-update ordering or the
public per-cell transport ledger.

## Dense recognition and allocation

`MembraneTransport.proteome_activities()` evaluates all canonical transport
targets in one sorted traversal of a cell's realized proteome. Each target's
Float64 accumulator still receives cohorts in the exact historical
locus/signature order, so interleaving target accumulators does not reorder any
individual sum. A per-engine derived cache stores the ten pure affinity values
for each encountered protein signature, replacing ten repeated hash lookups per
cohort with one reusable numeric row.

The simulation allocator uses the existing canonical
`SECONDARY_EXTRACELLULAR_IDS` order to index:

- realized protein activities;
- desired and ATP-scaled exchanges;
- extracellular scarcity totals and scales;
- aggregate world imports and exports;
- final per-cell exchange.

Only the hot internal representation changes. Before returning, the allocator
reconstructs the same Dictionary ledger consumed by snapshots, tests and
observers.

## Exact reference mode

`SimulationEngine._allocate_secondary_membrane_transport_legacy_reference()`
preserves the M7 Dictionary algorithm as the shadow oracle.
`SimConfig.secondary_transport_use_dense_allocator` selects the path for paired
tests and benchmarks; new simulations use the dense allocator by default. The
switch is execution-only and is excluded from authoritative snapshot data.

The dense implementation preserves:

- canonical metabolite order;
- sorted locus and protein-signature accumulation order;
- common pre-exchange world/cell snapshot semantics;
- proportional ATP and extracellular scarcity scaling;
- aggregate world imports before exports become available;
- molecule and adenylate conservation;
- returned ledger structure and insertion order.

## Paired benchmark

Run identical states through the dense and frozen M7 allocators:

```bash
godot --headless --path . --script benchmarks/p3b_dense_secondary_transport.gd -- \
  --warmup=10 \
  --ticks=30 \
  --populations=1,16,64,256,1000 \
  --output=p3b-dense-secondary-transport.json
```

The runner alternates which allocator executes first on every measured tick,
checks the complete transport ledger and exact checksum after each tick, and
checks final event-history equality. It exits nonzero if any scientific state
diverges or the report cannot be written.

## Acceptance gate

P3-B is eligible to merge only when:

- dense activities exactly equal scalar M7 activities for every transport
  target;
- ATP scaling, active exchange ledgers, full checksums and event histories are
  exactly equal across both paths;
- all M0–M10 and P0–P3-A regressions pass;
- the paired benchmark completes without a scientific divergence;
- three sequential target-laptop runs show at least 10% lower median tick time
  at 64, 256 and 1,000 cells;
- no 1- or 16-cell median regression exceeds 5%.

If the gate is not met, the dense allocator remains unmerged. Scientific
tolerances and biological work are never reduced to manufacture a speedup.
