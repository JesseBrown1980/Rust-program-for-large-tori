# A Rust program for large tori

Certified lower bounds for the gap above the eight flat states, E₉ − E₈, of the S₃ lattice gauge Hamiltonian H = J_E (E + r M) on periodic L_x × L_y lattices too large for the
existing tools: 7×7, 8×8 and beyond.

This repository belongs to Jesse Daniel Brown's finite Yang–Mills project and ties back to its main repository,
[Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C). The method is that
repository's no-store hybrid certificate
([`VECTOR-BOUND.md`](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/blob/main/native/ym/results/3x4_stream/perron/vector/VECTOR-BOUND.md));
this repository is a new implementation of it without the two size limits of the present programs. Written by seat ACER-CLAUDE-FABLE52 at his direction, starting November 15, 2026.

## Status: first certificates on 6×6, 7×7 and 8×8, as this program's own; independent check in progress

Lower ends of E₉ − E₈ in J_E, measured on GitHub cells on November 15, 2026 with predictions pushed first ([`PREDICTION-20261115.hbp`](PREDICTION-20261115.hbp), [`RESULT-20261115.hbp`](RESULT-20261115.hbp)). On every row each lump inside the cut holds exactly its count of configurations.

| Torus | r | Cut | Classes | E₉ − E₈ ≥ | 3r − 21/r | Time |
|---|---:|---:|---:|---:|---:|---|
| 6×6 | 9 | 11 | 2,257,839 | **16.904567** | 24.67 | 35 min |
| 7×7 | 12 | 10 | 337,082 | **28.965224** | 35.25 | 1 min 40 s |
| 7×7 | 16 | 10 | 337,082 | **45.207828** | 46.69 | 2 min 20 s |
| 8×8 | 16 | 10 | 766,866 | **44.345382** | 46.69 | 5 min 21 s |
| 8×8 | 24 | 10 | 766,866 | **70.648926** | 71.13 | 6 min 31 s |

**What these are, and are not.** Each is a certificate printed by this one program for one finite torus. The program passed its gates G1 to G3 on 4×4 and 5×5 against the registered programs ([`GATES-LOCAL-20261115.hbp`](GATES-LOCAL-20261115.hbp): identical class, configuration and move counts, certificates within 0.02). GPT-6, on another model, has written a separate checker (pull request 4) and verified a 4×4 vector with it; the 6×6, 7×7 and 8×8 vectors have **not** been checked independently yet. The cells of `dump-cells.yml` rerun these points with the vector written out for that check. Gate G4 gave no comparison: the widened copies of the earlier tools return valid bounds far below the flat level on 6×6 (80.0 J_E where the flat level is 216), so they cannot confirm or contradict these numbers.

Predictions: P1, P2, P4 met; P5 missed high (8×8 at r = 24 came in above the predicted band); P6 missed high on 8×8 (more classes than predicted); P3 pending.

    large_tori LX LY R CUT THREADS SWEEPS auto [REPLAY_R,...] [dump=FILE]

## Why a new program

The present tools (`ritz` and `hybridcheck` in the main repository) certified 4×4, 4×5, 5×5, 4×6 and 5×6, and widened copies reach 48 plaquettes. Two things stop them above that:

| Limit | Where | What it stops |
|---|---|---|
| A configuration is one base-6 number of P + 1 digits held in 128 bits (49 digits) | both programs, keys and file format | 7×7 needs 50 digits, 8×8 needs 65 |
| The rows above the cut come from a table of every integer solution: 29,422 rows on 4×4, 465,479 on 5×5, 4,460,235 on 6×6 (6.3 GB to build) | `local_table.py` | 7×7 and 8×8 by size |

## Design

1. **Keys are digit arrays.** One S₃ element per free link after gauge fixing; the canonical form is the least array under gauge, global conjugation and the lattice symmetries, with the
   plaquette-class pattern used to prune the symmetries tried. No packing, so no size limit from the key.
2. **Rows at M ≤ cut are exact.** Every gauge orbit with 0 < M ≤ cut is enumerated from the definition, all 3·nl moves of each are applied and every plaquette recomputed. The count of
   configurations in each lump (a, b) is checked against C(P; a, b)·N(a, b), with N = 12·3^a·2^b for even a ≥ 2 and N = 12·2^b + 6(−1)^b for a = 0.
3. **Rows above the cut use extreme rows only.** For a fixed lump the ratio (H′d)/d is linear in the numbers of links of each pattern, and those numbers obey the three constraints of
   the main repository's table (four links per plaquette, a even, the adjacency caps). The minimum of a linear form over a polytope is at a vertex, so the program evaluates the vertices
   of each lump's polytope in exact rationals. That is a lower bound on the minimum over the integer rows, so the certificate stays valid; it can be slightly lower than the full table
   gives. The main repository lists this generator as an open item.
4. **The certificate is an exact integer minimum.** The vector is fitted by resolvent sweeps at a list of targets; any positive vector gives a valid bound (Barta), so the fit needs no
   trust. E₈ from above is 3·nl + 2r − √(4r² + 3·nl) for min(L_x, L_y) ≥ 4, as registered.
5. Rust 1.81, `clippy -D warnings`, integers only, output rows ending `json=0`.

## Gates before any 7×7 or 8×8 run

| Gate | Requirement |
|---|---|
| G1 | 4×4 at r = 9, cut 10: the exact-row minimum equals the registered programs' value, and the certificate with the full integer table is 116.216456 as measured on November 15; the extreme-row certificate is recorded beside it |
| G2 | 4×4 at r = 6, cut 14 and 5×5 at r = 9, cut 11: certificates at or below the registered 103.521056 and 168.020667 and within a stated distance of them |
| G3 | Every lump count equals C(P; a, b)·N(a, b) on 4×4, 5×5, 6×6 |
| G4 | 6×6 at r = 9, cut 11: agreement with the widened copies of the present tools |
| G5 | Checks that go red: a claim raised by one unit, a dropped representative, a wrong M |

After the gates: predictions pushed first, then 7×7 at r = 12 and 16 and 8×8 at r = 16 on GitHub cells, a second machine as witness, and a request to GPT-6 for an independent checker,
since one new program otherwise checks itself.

## Scope

Each result will be a statement about one finite torus. A threshold independent of size is the stability route's claim in the main repository; its constant is not available as a
number.
