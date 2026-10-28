#!/usr/bin/env python3
"""Exact integer-arithmetic replay of the eight-vertex P2604 bound.

Prerequisites:
  g++ -O3 -std=c++17 enumerate_graphs_8.cpp -o enumerate_graphs_8
  ./enumerate_graphs_8 > graphs8_masks.txt
  python3 build_certificates.py   # once, requires numpy and scipy
  python3 verify_exact.py          # this step: standard library only

The numerical LP algorithm is NOT trusted for verification: every supplied
rational dual is checked against all 28 edge inequalities using integers.
"""
import gzip
import itertools
import json
import os

N = 8
pairs = list(itertools.combinations(range(8), 2))
cuts = []
for k in (2, 3, 4):
    for vertices in itertools.combinations(range(8), k):
        s = sum(1 << i for i in vertices)
        if k == 4 and not (s & 1):
            continue
        cuts.append(s)
assert len(pairs) == 28 and len(cuts) == 119


def distances(g):
    d = [[99] * 8 for _ in range(8)]
    for i in range(8):
        d[i][i] = 0
    for k, (i, j) in enumerate(pairs):
        if (g >> k) & 1:
            d[i][j] = d[j][i] = 1
    for k in range(8):
        for i in range(8):
            for j in range(8):
                d[i][j] = min(d[i][j], d[i][k] + d[k][j])
    assert all(d[i][j] < 99 for i in range(8) for j in range(8))
    return d


def tsp_cost(d):
    """Exact Held-Karp recurrence with designated starting vertex 0."""
    dp = [[999] * 8 for _ in range(256)]
    dp[1][0] = 0
    for subset in range(1, 256, 2):
        for i in range(8):
            base = dp[subset][i]
            if base == 999:
                continue
            for j in range(1, 8):
                if not (subset >> j) & 1:
                    nxt = subset | (1 << j)
                    dp[nxt][j] = min(dp[nxt][j], base + d[i][j])
    return min(dp[255][i] + d[i][0] for i in range(1, 8))


def verify_dual(g, t, D, entries, d):
    """Primal min c.x with (-cut)x<=-2, degrees Ex=2, x>=0.

    Dual variables u<=0 for <= inequalities and v unrestricted for
    equalities obey (-cut)^T u + E^T v <= c. We verify after multiplying
    all dual multipliers by their shared positive denominator D.
    """
    assert isinstance(D, int) and D > 0
    z = [0] * 127
    for i, a in entries:
        assert 0 <= i < 127 and z[i] == 0 and a != 0
        z[i] = a
    u, v = z[:119], z[119:]
    assert max(u) <= 0
    for a, b in pairs:
        crosses = sum(
            u[k] for k, S in enumerate(cuts)
            if bool((S >> a) & 1) != bool((S >> b) & 1)
        )
        left = v[a] + v[b] - crosses
        assert left <= D * d[a][b], (g, a, b, left, D * d[a][b])
    objective = 2 * (sum(v) - sum(u))
    assert 10 * objective >= 9 * t * D, (g, t, objective, D)
    return objective, D


def verify_witness():
    # Three disjoint length-three paths, hubs 6 and 7.
    branches = [(5, 0), (4, 1), (3, 2)]
    edges = set()
    for a, b in branches:
        edges.update({
            tuple(sorted((6, a))), tuple(sorted((a, b))),
            tuple(sorted((b, 7)))
        })
    mask = sum(1 << k for k, p in enumerate(pairs) if p in edges)
    assert mask == 43135568
    augmented = mask | (1 << pairs.index((6, 7)))
    assert augmented == 177353296

    for g in (mask, augmented):
        d = distances(g)
        assert tsp_cost(d) == 10

        # Primal candidate with denominator 6.
        x6 = {p: 0 for p in pairs}
        for a, b in branches:
            x6[tuple(sorted((6, a)))] = 4
            x6[tuple(sorted((a, b)))] = 6
            x6[tuple(sorted((b, 7)))] = 4
        for side in ([5, 4, 3], [0, 1, 2]):
            for a, b in itertools.combinations(side, 2):
                x6[tuple(sorted((a, b)))] = 1
        assert all(
            sum(x for (a, b), x in x6.items() if v in (a, b)) == 12
            for v in range(8)
        )
        for S in cuts:
            assert sum(
                x for (a, b), x in x6.items()
                if bool((S >> a) & 1) != bool((S >> b) & 1)
            ) >= 12
        assert sum(x * d[a][b] for (a, b), x in x6.items()) == 54

        # Hand-checkable dual value 9: six interior degree constraints
        # at 1/2 and three two-vertex branch cuts at 1/2.
        branch_sets = [(1 << a) | (1 << b) for a, b in branches]
        assert branch_sets == [33, 18, 12]
        for a, b in pairs:
            lhs2 = (int(a < 6) + int(b < 6)
                    + sum(bool((S >> a) & 1) != bool((S >> b) & 1)
                          for S in branch_sets))
            assert lhs2 <= 2 * d[a][b]
    return mask


def main():
    root = os.path.dirname(__file__)
    graph_file = os.path.join(root, 'graphs8_masks.txt')
    cert_file = os.path.join(root, 'p2604_duals.json.gz')
    if not os.path.exists(graph_file) or not os.path.exists(cert_file):
        raise SystemExit(
            'Generate graph listing and rational dual archive first; '
            'see README.md reproduction steps.'
        )

    with open(graph_file) as f:
        keys = [int(s) for s in f]
    assert len(keys) == 11117
    assert len(set(keys)) == 11117 and keys == sorted(keys)
    with gzip.open(cert_file, 'rt') as f:
        raw = json.load(f)
    cert = {key: (t, D, entries) for key, t, D, entries in raw}
    assert len(raw) == len(cert) == 4921

    n_ham = n_nonham = 0
    saturated = []
    for g in keys:
        d = distances(g)
        t = tsp_cost(d)
        if t == 8:
            n_ham += 1
            assert g not in cert
        else:
            n_nonham += 1
            assert g in cert
            ct, D, entries = cert[g]
            assert ct == t
            numerator, denom = verify_dual(g, t, D, entries, d)
            if 10 * numerator == 9 * t * denom:
                saturated.append(g)
    assert n_ham == 6196 and n_nonham == 4921
    assert set(cert).issubset(set(keys))

    winner = verify_witness()
    assert winner in saturated
    print('PASS: 11,117 connected graph isomorphism classes')
    print('PASS: 6,196 with Hamiltonian tour of length 8')
    print('PASS: 4,921 rational dual certificates verified exactly')
    print('PASS: every graph has LP >= 9/10 times its exact optimal tour')
    print('PASS: theta(3,3,3) has exact tour cost 10 and LP cost 9')
    print('CERTIFIED MAXIMUM: 10/9')
    print('Certificate-bound equality keys:', saturated)


if __name__ == '__main__':
    main()
