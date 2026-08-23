# M7-D reciprocal cross-feeding causal assay

M7-C demonstrated a one-way molecular handoff through the existing `R04 -> W1 -> R05` path. M7-D tests whether two chemically distinct waste/recovery routes can operate simultaneously in opposite directions without introducing a mutualism, partner, cooperation, or exchange behavior API.

This remains a controlled capability experiment using hand-constructed realized proteins.

## Reciprocal chemistry

Lineage A carries:

- R04-compatible catalytic protein: `C3 + NADH + ADP -> W1 + NAD + ATP`;
- W1-compatible membrane protein for reversible W1 exchange;
- W2-compatible membrane protein for reversible W2 exchange;
- R11-compatible catalytic protein: `W2 + O2 + ADP -> C2 + CO2 + ATP + ROS`.

Lineage B carries:

- R10-compatible catalytic protein: `C3 + NADH -> W2 + NAD`;
- W2-compatible membrane protein;
- W1-compatible membrane protein;
- R05-compatible catalytic protein: `W1 + NAD + ADP -> C2 + CO2 + NADH + ATP`.

Thus the complete controlled loop is:

`A R04 -> W1 -> environment -> B -> R05 -> C2 + ATP`

and simultaneously:

`B R10 -> W2 -> environment -> A -> R11 -> C2 + ATP + ROS`.

The two directions use different metabolites and different recovery reactions. No exchange value, partner identity, ecological role, or combined reciprocal-reward term exists.

## Snapshot semantics

As in M7-B/C, the first membrane phase sees a pre-exchange extracellular snapshot. A and B can export W1/W2, but neither lineage can consume the partner's zero-time export in that same phase. The subsequent membrane phase sees those extracellular molecules and permits gradient-driven uptake.

Both W1 and W2 must be conserved across cells plus chamber during the membrane handoff.

## Causal ablations

Reciprocity is not inferred merely because both lineages benefit in one configuration. Two directional ablations isolate causality:

- omit B's R10 production while retaining A's W1 production: B may still recover W1, but A must receive no W2 and R11 must remain zero;
- omit A's R04 production while retaining B's W2 production: A may still recover W2, but B must receive no W1 and R05 must remain zero.

These controls demonstrate two independent material dependencies rather than a hidden pairwise bonus.

## Explicit tradeoff

R11 produces ROS while recovering W2. The reciprocal benefit to lineage A therefore includes an ordinary modeled oxidative side effect. The assay requires ROS to increase along with C2/ATP. This prevents the reciprocal interaction from becoming a cost-free ecological reward.

## Reproducibility

A fixed starting state/seed must reproduce both transport ledgers and the complete post-recovery simulation checksum exactly.

## Interpretation boundary

Passing M7-D means the engine can represent **reciprocal metabolic exchange** using generic molecular mechanisms alone. It does not establish that reciprocal exchange is evolutionarily stable, that both lineages reproduce faster together, or that mutualism evolves spontaneously.

The remaining M7 exit work still includes costly extracellular catalysis/public-good exploitation, detoxification, damaging extracellular chemistry/self-resistance, dependence after biosynthetic loss, growth-level secondary-resource controls, and any required genome-structure capability.
