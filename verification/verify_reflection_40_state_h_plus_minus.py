#!/usr/bin/env python3
"""Pure-standard-library verifier for the reflected CEN/UV H^{+-} interface.

The verifier reconstructs the six symbolic source skeletons, enumerates all
256 subsets of the eight C4-gadget edges, and classifies every resulting
support at each of the standard/reflected CEN and UV sites.  It then checks:

* the original 22-state domain and its reflected failure points;
* the least forward closure of that domain under all four maps (40 states);
* the required 0/1 local separation for every closure state and retained edge;
* all 56 ordered pairs of gadget edges are separated by some source-to-target
  interface pattern at every site (the final Hamilton-cycle rows (7.15) and the
  older-edge/final-edge conditions (7.15a) are checked in check_final.py); and
* the birth state of every retained gadget edge.

No graph package, random search, solver, or project module is imported.
"""

from __future__ import annotations

from collections import deque
from itertools import product


ROLES = ("A", "B", "C", "D", "E", "L")
ROLE_ORDER = {role: index for index, role in enumerate(ROLES)}
GADGET = (
    ("u", "a"),
    ("a", "b"),
    ("b", "v"),
    ("x", "c"),
    ("c", "d"),
    ("d", "y"),
    ("a", "d"),
    ("b", "c"),
)
EDGE_NAMES = ("ua", "ab", "bv", "xc", "cd", "dy", "ad", "bc")
SITES = {
    "CEN": ("b", "c", "d", "a"),
    "UV": ("u", "a", "b", "v"),
    "CEN_ref": ("d", "a", "b", "c"),
    "UV_ref": ("b", "v", "u", "a"),
}


def edge(left: str, right: str) -> frozenset[str]:
    return frozenset((left, right))


def path_skeleton(order: tuple[str, ...], prefix: str) -> tuple[set[str], set[frozenset[str]]]:
    nodes = set(order)
    edges: set[frozenset[str]] = set()
    for index, (left, right) in enumerate(zip(order, order[1:])):
        middle = f"{prefix}{index}"
        nodes.add(middle)
        edges.add(edge(left, middle))
        edges.add(edge(middle, right))
    return nodes, edges


def source_skeletons() -> dict[str, tuple[set[str], set[frozenset[str]]]]:
    result = {
        role: path_skeleton(order, role)
        for role, order in {
            "B": ("u", "x", "y", "v"),
            "C": ("x", "u", "v", "y"),
            "D": ("v", "u", "y", "x"),
            "E": ("u", "v", "x", "y"),
        }.items()
    }
    for role, pairs in {
        "A": (("u", "y"), ("v", "x")),
        "L": (("u", "v"), ("x", "y")),
    }.items():
        nodes = {"u", "v", "x", "y"}
        edges: set[frozenset[str]] = set()
        for index, (left, right) in enumerate(pairs):
            middle = f"{role}{index}"
            nodes.add(middle)
            edges.add(edge(left, middle))
            edges.add(edge(middle, right))
        result[role] = nodes, edges
    return result


def adjacency(nodes: set[str], edges: set[frozenset[str]]) -> dict[str, set[str]]:
    result = {node: set() for node in nodes}
    for item in edges:
        left, right = tuple(item)
        result[left].add(right)
        result[right].add(left)
    return result


def components(nodes: set[str], edges: set[frozenset[str]]) -> list[set[str]]:
    graph = adjacency(nodes, edges)
    unseen = set(nodes)
    result: list[set[str]] = []
    while unseen:
        start = next(iter(unseen))
        unseen.remove(start)
        component = {start}
        queue = [start]
        while queue:
            current = queue.pop()
            for neighbour in graph[current]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    queue.append(neighbour)
        result.append(component)
    return result


def endpoint_components(
    nodes: set[str], edges: set[frozenset[str]]
) -> tuple[bool, list[frozenset[str]]]:
    """Validate path components and return their two degree-one endpoints."""

    graph = adjacency(nodes, edges)
    endpoints: list[frozenset[str]] = []
    for component in components(nodes, edges):
        degrees = {node: len(graph[node]) for node in component}
        if any(degree not in (1, 2) for degree in degrees.values()):
            return False, []
        ends = frozenset(node for node, degree in degrees.items() if degree == 1)
        if len(ends) != 2:
            return False, []
        endpoints.append(ends)
    return True, endpoints


def classify_support(
    nodes: set[str], support: set[frozenset[str]], site: tuple[str, str, str, str]
) -> set[str]:
    """Classify every usable target role at one oriented site."""

    u, v, x, y = site
    p, q = edge(u, v), edge(x, y)
    result: set[str] = set()

    full_graph = adjacency(nodes, support)
    full_connected = len(components(nodes, support)) == 1
    full_cycle = full_connected and all(len(full_graph[node]) == 2 for node in nodes)
    if full_cycle:
        state = (p in support, q in support)
        if state == (True, False):
            result.add("B")
        elif state == (False, True):
            result.add("C")
        elif state == (True, True):
            ok, ends = endpoint_components(nodes, support - {p, q})
            target_pairing = frozenset(
                (frozenset((u, y)), frozenset((v, x)))
            )
            if ok and frozenset(ends) == target_pairing:
                result.add("A")

    cut_support = support - {p, q}
    ok, ends = endpoint_components(nodes, cut_support)
    if ok and len(ends) == 1:
        if ends[0] == frozenset((v, x)):
            result.add("D")
        if ends[0] == frozenset((u, y)):
            result.add("E")
    if ok and len(ends) == 2:
        target_pairing = frozenset(
            (frozenset((u, v)), frozenset((x, y)))
        )
        if frozenset(ends) == target_pairing:
            result.add("L")
    return result


def local_patterns(
    site: tuple[str, str, str, str]
) -> dict[tuple[str, str], set[int]]:
    patterns = {(source, target): set() for source in ROLES for target in ROLES}
    for source, (old_nodes, old_edges) in source_skeletons().items():
        nodes = old_nodes | {"a", "b", "c", "d"}
        for mask in range(1 << len(GADGET)):
            support = set(old_edges)
            support.update(
                edge(left, right)
                for index, (left, right) in enumerate(GADGET)
                if mask & (1 << index)
            )
            for target in classify_support(nodes, support, site):
                patterns[source, target].add(mask)
    return patterns


def relation(patterns: dict[tuple[str, str], set[int]]) -> dict[str, frozenset[str]]:
    """Return target <- source role relation."""

    return {
        target: frozenset(
            source for source in ROLES if patterns[source, target]
        )
        for target in ROLES
    }


def transfer_state(
    state: frozenset[str], target_sources: dict[str, frozenset[str]]
) -> frozenset[str]:
    return frozenset(
        target
        for target, sources in target_sources.items()
        if state & sources
    )


def parse_states() -> set[tuple[frozenset[str], frozenset[str]]]:
    tokens = (
        "ABC/ADEL",
        "ABCD/AEL",
        "ABCE/ADL",
        "ACDEL/AB",
        "ABDEL/AC",
        "ABCDEL/ABCDEL",
        "ABCDEL/ABCEL",
        "ABCDEL/ABDEL",
        "ABCDEL/ABEL",
        "ABCDEL/ADEL",
        "ABCDEL/AEL",
        "ABCEL/ABCDEL",
        "ABCEL/ABDEL",
        "ABCEL/ADEL",
        "ABDEL/ABCDEL",
        "ABDEL/ABCEL",
        "ABDEL/ACDEL",
        "ABDEL/ACEL",
        "ACDEL/ABDEL",
        "ACDEL/ABEL",
        "ADEL/ABCDEL",
        "ADEL/ABCEL",
    )
    result = set()
    for token in tokens:
        positive, negative = token.split("/")
        result.add((frozenset(positive), frozenset(negative)))
    assert len(result) == 22
    return result


def state_name(state: tuple[frozenset[str], frozenset[str]]) -> str:
    positive, negative = state
    return "".join(sorted(positive, key=ROLE_ORDER.get)) + "/" + "".join(
        sorted(negative, key=ROLE_ORDER.get)
    )


def state_sort_key(state: tuple[frozenset[str], frozenset[str]]) -> tuple[object, ...]:
    positive, negative = state
    return (
        len(positive) + len(negative),
        len(positive),
        len(negative),
        state_name(state),
    )


def closure(
    seed: set[tuple[frozenset[str], frozenset[str]]],
    maps: dict[str, dict[str, frozenset[str]]],
) -> set[tuple[frozenset[str], frozenset[str]]]:
    result = set(seed)
    queue = deque(seed)
    while queue:
        positive, negative = queue.popleft()
        for target_sources in maps.values():
            image = (
                transfer_state(positive, target_sources),
                transfer_state(negative, target_sources),
            )
            if image not in result:
                result.add(image)
                queue.append(image)
    return result


def selected_indices(site: tuple[str, str, str, str]) -> tuple[int, ...]:
    selected = {edge(site[0], site[1]), edge(site[2], site[3])}
    return tuple(index for index, item in enumerate(GADGET) if edge(*item) in selected)


def has_mask(
    patterns: dict[tuple[str, str], set[int]],
    sources: frozenset[str],
    edge_index: int,
    bit: int,
) -> bool:
    return any(
        any(((mask >> edge_index) & 1) == bit for mask in patterns[source, target])
        for source in sources
        for target in ROLES
    )


def birth_state(
    patterns: dict[tuple[str, str], set[int]], edge_index: int
) -> tuple[frozenset[str], frozenset[str]]:
    positive = frozenset(
        target
        for target in ROLES
        if any(
            any((mask >> edge_index) & 1 for mask in patterns[source, target])
            for source in ROLES
        )
    )
    negative = frozenset(
        target
        for target in ROLES
        if any(
            any(not ((mask >> edge_index) & 1) for mask in patterns[source, target])
            for source in ROLES
        )
    )
    return positive, negative


def assert_relation(
    actual: dict[str, frozenset[str]], expected: dict[str, str]
) -> None:
    assert actual == {
        target: frozenset(sources)
        for target, sources in expected.items()
    }, (actual, expected)


def main() -> None:
    all_patterns = {name: local_patterns(site) for name, site in SITES.items()}
    maps = {name: relation(patterns) for name, patterns in all_patterns.items()}

    assert_relation(
        maps["CEN"],
        {"A": "BCL", "B": "E", "C": "D", "D": "A", "E": "A", "L": "A"},
    )
    assert_relation(
        maps["UV"],
        {"A": "ABL", "B": "E", "C": "D", "D": "C", "E": "A", "L": "A"},
    )
    assert_relation(
        maps["CEN_ref"],
        {"A": "BCL", "B": "D", "C": "E", "D": "A", "E": "A", "L": "A"},
    )
    assert_relation(
        maps["UV_ref"],
        {"A": "ABL", "B": "D", "C": "E", "D": "A", "E": "C", "L": "A"},
    )

    sigma = parse_states()
    standard_failures = {
        name: [
            state
            for state in sigma
            if (
                transfer_state(state[0], maps[name]),
                transfer_state(state[1], maps[name]),
            )
            not in sigma
        ]
        for name in ("CEN", "UV")
    }
    assert all(not failures for failures in standard_failures.values())

    reflected_failures = {
        name: [
            (
                state,
                (
                    transfer_state(state[0], maps[name]),
                    transfer_state(state[1], maps[name]),
                ),
            )
            for state in sigma
            if (
                transfer_state(state[0], maps[name]),
                transfer_state(state[1], maps[name]),
            )
            not in sigma
        ]
        for name in ("CEN_ref", "UV_ref")
    }
    assert len(reflected_failures["CEN_ref"]) == 10
    assert len(reflected_failures["UV_ref"]) == 21
    cen_min = min(reflected_failures["CEN_ref"], key=lambda item: state_sort_key(item[0]))
    uv_min = min(reflected_failures["UV_ref"], key=lambda item: state_sort_key(item[0]))
    assert state_name(cen_min[0]) == "ADEL/ABCEL"
    assert state_name(cen_min[1]) == "ABCDEL/ACDEL"
    assert state_name(uv_min[0]) == "ABC/ADEL"
    assert state_name(uv_min[1]) == "ADEL/ABCDL"

    all_maps = {"CEN": maps["CEN"], "UV": maps["UV"], "CEN_ref": maps["CEN_ref"], "UV_ref": maps["UV_ref"]}
    sigma_star = closure(sigma, all_maps)
    assert len(sigma_star) == 40
    assert len(sigma_star - sigma) == 18
    assert all(positive and negative for positive, negative in sigma_star)
    assert closure(sigma_star, all_maps) == sigma_star

    # Every source role and every target role is represented in every map.
    # The source-side assertion is the part needed to carry a concrete
    # joint-witness separation through an annulus; the target-side assertion
    # records role availability for the next layer.
    for target_sources in all_maps.values():
        assert all(target_sources.values())
        assert all(
            any(source in sources for sources in target_sources.values())
            for source in ROLES
        )

    separation_counts: dict[str, int] = {}
    item2_counts: dict[str, int] = {}
    joint_pair_counts: dict[str, int] = {}
    for name, patterns in all_patterns.items():
        site = SITES[name]
        retained = tuple(index for index in range(len(GADGET)) if index not in selected_indices(site))
        checked = 0
        for state in sigma_star:
            positive, negative = state
            for edge_index in retained:
                assert has_mask(patterns, positive, edge_index, 0)
                assert has_mask(patterns, negative, edge_index, 1)
                checked += 1
        item2_counts[name] = checked

        # A pair of gadget edges which survives as a retained pair must have
        # one concrete source/target transition with local bits (1,0).
        # This is the joint-witness condition at the birth layer; marginal
        # 0/1 checks above do not imply it by themselves.
        joint_checked = 0
        for required, forbidden in product(retained, repeat=2):
            if required == forbidden:
                continue
            assert any(
                any(
                    ((mask >> required) & 1) == 1
                    and ((mask >> forbidden) & 1) == 0
                    for mask in patterns[source, target]
                )
                for source in ROLES
                for target in ROLES
            )
            joint_checked += 1
        joint_pair_counts[name] = joint_checked

        separated = 0
        for required, forbidden in product(range(len(GADGET)), repeat=2):
            if required == forbidden:
                continue
            if any(
                any(
                    ((mask >> required) & 1) == 1
                    and ((mask >> forbidden) & 1) == 0
                    for mask in patterns[source, target]
                )
                for source in ROLES
                for target in ROLES
            ):
                separated += 1
        assert separated == 56
        separation_counts[name] = separated

    expected_births = {
        "CEN": {
            "ua": "ABDEL/AC",
            "ab": "ABCD/AEL",
            "bv": "ACDEL/AB",
            "xc": "ACDEL/AB",
            "cd": "ABCE/ADL",
            "dy": "ABDEL/AC",
        },
        "CEN_ref": {
            "ua": "ACDEL/AB",
            "ab": "ABCE/ADL",
            "bv": "ABDEL/AC",
            "xc": "ABDEL/AC",
            "cd": "ABCD/AEL",
            "dy": "ACDEL/AB",
        },
        "UV": {
            "ab": "ABCE/ADL",
            "xc": "ACDEL/AB",
            "cd": "ABC/ADEL",
            "dy": "ABDEL/AC",
            "ad": "ACDEL/AB",
            "bc": "ABDEL/AC",
        },
        "UV_ref": {
            "ab": "ABCD/AEL",
            "xc": "ABDEL/AC",
            "cd": "ABC/ADEL",
            "dy": "ACDEL/AB",
            "ad": "ABDEL/AC",
            "bc": "ACDEL/AB",
        },
    }
    birth_records: dict[str, dict[str, str]] = {}
    for name, patterns in all_patterns.items():
        actual: dict[str, str] = {}
        for index in range(len(GADGET)):
            if index in selected_indices(SITES[name]):
                continue
            state = birth_state(patterns, index)
            assert state in sigma_star
            actual[EDGE_NAMES[index]] = state_name(state)
        assert actual == expected_births[name]
        birth_records[name] = actual

    print("original_sigma_states=22")
    print("reflected_sigma_failures=CEN_ref:10,UV_ref:21")
    print("minimal_CEN_ref_failure=ADEL/ABCEL->ABCDEL/ACDEL")
    print("minimal_UV_ref_failure=ABC/ADEL->ADEL/ABCDL")
    print(f"four_map_forward_closure={len(sigma_star)}")
    print(f"new_states_beyond_sigma={len(sigma_star - sigma)}")
    print("item2_checks=" + ",".join(f"{name}:{item2_counts[name]}" for name in SITES))
    print("retained_pair_joint_checks=" + ",".join(f"{name}:{joint_pair_counts[name]}" for name in SITES))
    print("interface_ordered_separation=" + ",".join(f"{name}:{separation_counts[name]}" for name in SITES))
    print("birth_states=" + repr(birth_records))
    print("REFLECTION_40_STATE_H_PLUS_MINUS_VERIFIED")


if __name__ == "__main__":
    main()
