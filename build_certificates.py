#!/usr/bin/env python3
"""Regenerate the exact-verifiable rational dual certificates for P2604.

Requires NumPy and SciPy to find candidate LP duals. The numerical solver's
answers are never trusted: each rationalized dual must pass verify_dual()
using integer arithmetic before being written. After this script finishes,
run verify_exact.py to prove the universal bound using Python stdlib alone.

Re-running may produce different but equally valid feasible rational duals
if the version of SciPy/HiGHS selects a different optimum.
"""
import gzip
import json
import math
import os
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

from verify_exact import cuts, distances, pairs, tsp_cost, verify_dual

root = os.path.dirname(__file__)
with open(os.path.join(root, 'graphs8_masks.txt')) as f:
    keys = [int(x) for x in f]

Au = np.asarray([
    [-int(bool((S >> i) & 1) != bool((S >> j) & 1)) for i, j in pairs]
    for S in cuts
], dtype=float)
E = np.asarray([
    [int(v == i or v == j) for i, j in pairs]
    for v in range(8)
], dtype=float)

result = []
for idx, g in enumerate(keys):
    d = distances(g)
    t = tsp_cost(d)
    if t == 8:
        continue
    c = np.asarray([d[i][j] for i, j in pairs], dtype=float)
    sol = linprog(
        c, A_ub=Au, b_ub=[-2.0] * 119,
        A_eq=E, b_eq=[2.0] * 8, bounds=(0, None), method='highs'
    )
    if not sol.success:
        raise RuntimeError(f'LP failed on graph {g}: {sol.message}')

    fractional = [
        Fraction(float(a)).limit_denominator(1000)
        if abs(float(a)) > 1e-8 else Fraction(0)
        for a in list(sol.ineqlin.marginals) + list(sol.eqlin.marginals)
    ]
    D = math.lcm(*(x.denominator for x in fractional))
    integers = [int(x * D) for x in fractional]
    entries = [[i, a] for i, a in enumerate(integers) if a]
    verify_dual(g, t, D, entries, d)  # exact: reject invalid rationalization
    result.append([g, t, D, entries])

    if idx % 1500 == 0:
        print('processed index', idx, flush=True)

assert len(keys) == 11117
assert len(result) == 4921
dest = os.path.join(root, 'p2604_duals.json.gz')
with gzip.open(dest, 'wt') as f:
    json.dump(result, f, separators=(',', ':'))
print('Wrote', dest, 'with', len(result), 'exactly verified dual certificates')
