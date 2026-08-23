# M6-B spatial broad phase

M6-A deliberately starts with an O(N^2) all-pairs contact search as a correctness reference. M6-B adds a uniform spatial hash without changing contact physics.

Each finite disk is inserted into every bucket touched by its axis-aligned bounding box. Consequently, any pair of overlapping disks must share at least one bucket even if a disk straddles bucket boundaries or becomes larger through growth. Candidate pairs are deduplicated and sorted canonically by immutable cell ID.

The mechanical solver can run either broad phase. Both feed the same snapshot-based Jacobi contact calculation and the same pair order for all physically relevant contacts. The indexed mode is accepted only if controlled initial states produce exactly the same final positions and overlap as the O(N^2) reference.

The default bucket width is 2.0 chemical-grid units. It is a performance parameter only: changing it may change candidate count but must not change mechanical results.
