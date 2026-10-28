# Exact eight-vertex graphic-metric subtour-LP gap — TheoremDB P2604

**Computer-assisted proof candidate — October 27, 2026**

**Proposed result:** The maximum integrality ratio of the symmetric traveling-salesperson subtour-elimination LP over shortest-path metrics of connected simple graphs on exactly eight vertices is

\[
\boxed{\gamma_8 = 10/9}.
\]

**Research status:** This is AI-assisted original computational research with a finite exact-arithmetic verification method. It is **not independently peer reviewed, institutionally accepted, or confirmed novel**. Please reproduce the computation, inspect the enumeration and dual-feasibility logic, and check recent prior art before citing this as an established solution.

**Problem statement:** [TheoremDB P2604](https://www.theoremdb.org/statements/P2604/).

## Definitions

For a connected simple graph \(H=(V,E)\) on \(|V|=8\), equip the complete graph on \(V\) with shortest-path edge costs \(d_H(i,j)\). Let \(T(H)\) be the shortest Hamiltonian cycle in this metric, and let \(P(H)\) be the optimum of

\[
\min \sum_{\{i,j\}} d_H(i,j)x_{ij}
\]
\[
\sum_{j\ne i} x_{ij}=2\quad (i\in V),
\qquad
\sum_{i\in S,j\notin S} x_{ij}\ge2
\quad(\emptyset\ne S\subsetneq V),
\qquad x\ge0.
\]

The target quantity is \(\gamma_8=\max_H T(H)/P(H)\).

## Extremal lower bound: \(\gamma_8\ge10/9\)

Consider the theta graph \(\Theta(3,3,3)\): two hubs \(u,v\), joined by three internally disjoint paths \(u-a_i-b_i-v\) of length three.

**Optimal tour \(T=10\).** A metric tour expands into a connected closed walk in the original graph. To visit the interior of each branch, one either traverses it between hubs (cost \(\ge3\)) or returns to the same hub (cost \(\ge4\)). A closed walk has an even number of hub-to-hub traversals. With none, cost \(\ge12\); with two, cost \(\ge3+3+4=10\). Two branch traversals and a doubled excursion into the third attain 10.

**Primal LP solution \(P\le9\).** Set \(x_{a_i b_i}=1\), \(x_{ua_i}=x_{b_iv}=2/3\), and \(x_{a_i a_j}=x_{b_i b_j}=1/6\) for distinct branches \(i<j\); all others zero. Its degrees equal two; the verifier checks all subtour cuts exactly. Its total metric cost is \(3+4+2=9\).

**Dual certificate \(P\ge9\).** Give each interior-vertex degree equality multiplier \(1/2\), and each of the three branch-interior cut inequalities for \(\{a_i,b_i\}\) multiplier \(1/2\). The resulting per-edge coefficient never exceeds its metric cost, as checked explicitly; the dual objective is \(6+3=9\). Thus the example has \(T/P=10/9\).

## Universal upper bound: exact finite verification

Any connected eight-vertex graph has a vertex whose removal leaves a connected seven-vertex graph. Exhaustively enumerate all labeled seven-vertex graphs, discard disconnected ones, canonically deduplicate under isomorphism, then attach a new vertex to every nonempty subset of the other seven vertices and canonically deduplicate again. The resulting 11,117 classes represent all connected eight-vertex graphs.

For every class:

1. Determine the 28 shortest-path edge costs by integer Floyd–Warshall.
2. Compute the integer optimal tour cost \(T\) by Held–Karp dynamic programming.
3. If \(T=8\), the subtour LP has \(P=8\): every metric edge costs at least one and the degree equations force \(\sum_e x_e=8\).
4. Otherwise, check a rational dual-feasible solution satisfying \(P\ge9T/10\). The independent verifier tests all 28 dual constraints after clearing denominators, using **integer arithmetic only**. The dual certificates can be regenerated using SciPy's LP solver and then verified without relying on floating-point correctness.

The local computation identified **6,196** classes with \(T=8\) and **4,921** remaining classes, each with a verified rational lower bound. Thus \(T/P\le10/9\) for every eight-vertex graph. With the theta witness, \(\gamma_8=10/9\), provided the enumerator and verifier are correct.

The two equality examples identified by the certificate are edge-mask keys \`43135568\` (theta graph) and \`177353296\` (theta graph plus hub-to-hub edge).

## Reproduction

Required: Python 3, C++17 compiler, NumPy and SciPy **for building certificates**. **Exact verification** uses only Python standard library and the generated certificate archive.

\`\`\`bash
g++ -O3 -std=c++17 enumerate_graphs_8.cpp -o enumerate_graphs_8
./enumerate_graphs_8 > graphs8_masks.txt
python3 -m pip install numpy scipy
python3 build_certificates.py
python3 verify_exact.py
\`\`\`

The C++ program should print class counts 853 (seven vertices) and 11,117 (eight vertices) to stderr. The certificate generator creates \`p2604_duals.json.gz\` from numerical LP solutions **only after checking the extracted rational duals with exact integer arithmetic**. The checker then verifies every class and the extremal graph and prints \`CERTIFIED MAXIMUM: 10/9\`.

**Note:** The complete original precomputed certificate archive and graph listing were produced separately in the research package. This repository contains the generators and checker so that all artifacts can be regenerated; if a precomputed \`p2604_duals.json.gz\` is later added here, the \`build_certificates.py\` step can be skipped. Do not claim a fresh replay was completed solely by publishing these files.

## Limits and invitation for review

- This is only about **shortest-path graph metrics on exactly eight vertices**. It does not settle the general graphic-TSP or metric-TSP integrality-gap conjectures.
- TheoremDB listed P2604 open in its reviewed August 1, 2026 record. That is not proof that the answer was unpublished anywhere as of October 27, 2026.
- Please audit the canonical enumeration, the LP dual sign convention, the source of every exact certificate, and the explicit theta primal and dual.
- Research developed with assistance from OpenAI's ChatGPT. Manuscript and implementation are made available for independent checking, not as a claim of external endorsement.
