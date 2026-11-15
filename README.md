# A Rust program for large tori

Certified lower bounds for the gap above the eight flat states, E₉ − E₈, of the S₃ lattice gauge Hamiltonian H = J_E (E + r M) on periodic L_x × L_y lattices too large for the
existing tools: 7×7, 8×8 and beyond.

This repository belongs to Jesse Daniel Brown's finite Yang–Mills project and ties back to its main repository,
[Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C). The method is that
repository's no-store hybrid certificate
([`VECTOR-BOUND.md`](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/blob/main/native/ym/results/3x4_stream/perron/vector/VECTOR-BOUND.md));
this repository is a new implementation of it without the two size limits of the present programs. Written by seat ACER-CLAUDE-FABLE52 at his direction, starting November 15, 2026.

## Status: all seven original 4×4–8×8 fitted-vector certificates independently verified

The producer's first large-torus cells were registered and predicted before their runs ([`PREDICTION-20261115.hbp`](PREDICTION-20261115.hbp), [`RESULT-20261115.hbp`](RESULT-20261115.hbp)). The original `HYBRID_BIG` fitted vectors were then frozen with SHA-256 in [dump run 43724776455](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43724776455). **A separately authored Rust 1.81 checker independently reconstructed the complete source of all seven original 4×4, 5×5, 6×6, 7×7 and 8×8 fitted-vector certificates.** The last 6×6 full-source receipt and exact source SHA-256 are filed on [issue #3](https://github.com/JesseBrown1980/Rust-program-for-large-tori/issues/3#issuecomment-6896926875); every certificate's fan-in refuses a claim raised by one integer micro-unit.

Lower ends of `E₉ − E₈` in units of `J_E`, with each original producer result and its independent verdict:

| Torus | r | Cut | Classes | E₉ − E₈ ≥ | 3r − 21/r | Producer time | Independent vector replay |
|---|---:|---:|---:|---:|---:|---|---|
| 6×6 | 9 | 11 | 2,257,839 | **16.904567** | 24.67 | 35 min | **PASS**, [64-shard run 43730041762](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43730041762) |
| 7×7 | 12 | 10 | 337,082 | **28.965224** | 35.25 | 1 min 40 s | **PASS**, [32-shard run 43725543227](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43725543227) |
| 7×7 | 16 | 10 | 337,082 | **45.207828** | 46.69 | 2 min 20 s | **PASS**, [run 43725543227](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43725543227) |
| 8×8 | 16 | 10 | 766,866 | **44.345382** | 46.69 | 5 min 21 s | **PASS**, [66-job run 43726083659](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43726083659) |
| 8×8 | 24 | 10 | 766,866 | **70.648926** | 71.13 | 6 min 31 s | **PASS**, [run 43726083659](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43726083659) |

### Separate full-source vector checker, actual frozen producer data

The [independent verifier merged through PR #4](https://github.com/JesseBrown1980/Rust-program-for-large-tori/pull/4), `src/bin/vectorcheck.rs`, has its own six-element S₃, link geometry, spanning-tree gauge, canonical class, orbit-size computation and exact electric-row enumerator. For high-cut rows it independently enumerates the relaxed rational vertices, taking the **largest actual low-cut orbit vector entry** for any transition into that low lump. Its final fan-in requires every disjoint source shard, the original uncompressed vector SHA-256, complete raw-move counts, exact matching minima, and **REFUSAL of the claim +1 micro-unit**.

| Original vector | Independently matched `E₉` floor (`J_E`) | Exact classes checked | Raw electric moves rebuilt | Verification |
|---|---:|---:|---:|---|
| 4×4, r=9, cut=10 | **116.272217** | 10,269 | 985,824 | [Original vector replay 43725375166](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43725375166) PASS |
| 5×5, r=9, cut=11 | **168.034092** | 473,677 | 71,051,550 | [Independent full check 43726465339](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43726465339) PASS |
| 6×6, r=9, cut=11 | **227.666667** | 2,257,839 | 487,693,224 | [Independent full check 43730041762](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43730041762) PASS |
| 7×7, r=12, cut=10 | **317.469462** | 337,082 | 99,102,108 | [Full check 43725543227](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43725543227) PASS |
| 7×7, r=16, cut=10 | **334.903559** | 337,082 | 99,102,108 | [Full check 43725543227](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43725543227) PASS |
| 8×8, r=16, cut=10 | **422.822056** | 766,866 | 294,476,544 | [Full check 43726083659](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43726083659) PASS |
| 8×8, r=24, cut=10 | **450.803001** | 766,866 | 294,476,544 | [Full check 43726083659](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43726083659) PASS |

These are **source-conditional finite-volume spectral lower bounds** based on the registered Perron/Barta and eight-flat-sector comparisons. The original 6×6 vector is **independently confirmed**, with 64 source partitions, all 2,257,839 classes, all 487,693,224 moves and 3,153 above-cut vertices per shard checked; the frozen uncompressed SHA-256 is `dca2181bf4eee2d78e1d3c43ea6d1e7af28a2d083437f0d54a4f5af993bdc633`. Widened earlier tools still provide only very loose valid bounds (80, 81 and −217 versus flat energy 216), so the older numerical G4 *comparison* remains uninformative even though the separate vector checker passed.

The [independent high-lump oracle PR #1](https://github.com/JesseBrown1980/Rust-program-for-large-tori/pull/1) **merged** as `adcee593`, and [independent checker PR #4](https://github.com/JesseBrown1980/Rust-program-for-large-tori/pull/4) **merged** as `e32b533c`. The G5 fault controls were rebased onto the current producer and [merged through PR #8](https://github.com/JesseBrown1980/Rust-program-for-large-tori/pull/8) as `d0c3f707`; earlier draft [#2](https://github.com/JesseBrown1980/Rust-program-for-large-tori/pull/2) was closed as superseded. All three components have passing hosted Rust/proof checks.

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

## Registered gates and post-run status

| Gate | Outcome |
|---|---|
| G1 | 4×4, r=9, cut=10 source/certificate regression PASS; the original fitted producer vector independently reproduced **116.272217** in full |
| G2 | Local controls 4×4, r=6, cut=14 and 5×5, r=9, cut=11 PASS against registered floors; the original 5×5 vector independently reproduced **168.034092** in full |
| G3 | Exact low-cut weighted configuration counts and 3×link electric-row counts PASS on registered smaller runs and on the independently verified 7×7 and 8×8 vector dumps |
| G4 | 6×6, r=9, cut=11 original vector **independently verified in full** ([run 43730041762](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43730041762)); the earlier widened legacy-tool floors are too loose for a useful direct numerical comparison, so that specific historical G4 comparison remains **inconclusive** |
| G5 | Real 4×4 source fault controls (overclaimed integer floor, deleted representative, altered magnetic count, missing electric move) PASS on [Rust 1.81 run 43730317133](https://github.com/JesseBrown1980/Rust-program-for-large-tori/actions/runs/43730317133); [current-producer PR #8](https://github.com/JesseBrown1980/Rust-program-for-large-tori/pull/8) merged, preserving per-level profiles and vector dumps |

The predictions were frozen before the initial 7×7/8×8 production cells. The independent source-check workflow now reads the producer's `rows.hbp` and `vector_*.hbp.gz` artifacts with pinned SHA and claim. The 7×7 and 8×8 per-level profiles show a repeated tight-lump sequence at the low magnetic levels; it remains an **empirical pattern** from these tested sizes and couplings. The 6×6 original dump, the 7×7 cut-11 variant and the cloud witness should be recorded from their own completed receipts.

## Scope

Each result will be a statement about one finite torus. A threshold independent of size is the stability route's claim in the main repository; its constant is not available as a
number.
