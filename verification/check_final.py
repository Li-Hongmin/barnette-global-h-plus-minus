#!/usr/bin/env python3
"""Direct check of the final layer of Theorem 8.2 against the Hamilton-cycle rows (7.15).

For each source role, enumerate all 2^8 gadget subsets and keep those whose union
with the symbolic source skeleton is a single Hamilton cycle.  Assert that each
role has exactly the row of (7.15).  Then assert, for every signed state X/Y in the
40-state closure and each of the eight final gadget positions j, that some role in
X has bit 0 at j (old edge in, new edge out) and some role in Y has bit 1 at j
(new edge in, old edge out); and that the rows separate all 56 ordered pairs of
final gadget edges.  Exits nonzero on failure.
"""
from itertools import product
from verify_reflection_40_state_h_plus_minus import (
    ROLES, GADGET, EDGE_NAMES, source_skeletons, edge, adjacency, components,
    parse_states, closure, local_patterns, relation, SITES, state_name)

ROWS_715 = {"A": "11111100", "B": "10101011", "C": "01010111",
            "D": "01111010", "E": "11001101", "L": "10110111"}

def final_rows():
    rows = {}
    for role, (old_nodes, old_edges) in source_skeletons().items():
        nodes = old_nodes | {"a", "b", "c", "d"}
        found = []
        for mask in range(1 << len(GADGET)):
            support = set(old_edges) | {edge(*GADGET[i]) for i in range(8) if mask >> i & 1}
            g = adjacency(nodes, support)
            if all(len(g[v]) == 2 for v in nodes) and len(components(nodes, support)) == 1:
                found.append("".join(str(mask >> i & 1) for i in range(8)))
        rows[role] = found
    return rows

def main():
    rows = final_rows()
    for role in ROLES:
        assert rows[role] == [ROWS_715[role]], (role, rows[role])
    pairs = 0
    for r, f in product(range(8), repeat=2):
        if r != f:
            assert any(ROWS_715[t][r] == "1" and ROWS_715[t][f] == "0" for t in ROLES), (r, f)
            pairs += 1
    maps = {n: relation(local_patterns(s)) for n, s in SITES.items()}
    sigma_star = closure(parse_states(), maps)
    assert len(sigma_star) == 40
    mixed = 0
    for X, Y in sorted(sigma_star, key=state_name):
        for j in range(8):
            assert any(ROWS_715[t][j] == "0" for t in X), (state_name((X, Y)), EDGE_NAMES[j], "old-in/new-out")
            assert any(ROWS_715[t][j] == "1" for t in Y), (state_name((X, Y)), EDGE_NAMES[j], "new-in/old-out")
            mixed += 2
    print(f"final_rows_match_7.15=6 new_new_pairs={pairs} old_new_conditions={mixed}")
    print("FINAL_LAYER_7_15_VERIFIED")

if __name__ == "__main__":
    main()
