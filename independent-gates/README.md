# GPT-6 independent high-lump extreme-row oracle (15 November 2026)

**Purpose:** Check one mathematical substitution in the proposed large-torus hybrid certificate independently of its future Rust implementation. The original source code in the Yang–Mills repository is not imported or reused. This checker is a gate for the *upper-energy lump treatment*, not a full spectral certificate, low-energy orbit census, or result for 7×7 or 8×8.

## Exact reduction

A nonflat configuration is characterized partly by `(a,b)`, with `a` the even number of transposition plaquettes and `b` the 3-cycle count. Among links joining defects let `x` be total T–T, `y` T–C, `z` C–C, and `u` the same-frame subset of T–T. The seven pattern counts are

```
n_TTs=u; n_TTd=x-u; n_TC=y; n_CC=z
n_eT=4a-2x-y
n_eC=4b-y-2z
n_ee=2 L^2 -4(a+b)+x+y+z    [square case]
```

The exact pattern-count constraints include `x>=0,y>=0,z>=0`,
`x<=adj(a)`, `z<=adj(b)`, `x+y+z<=adj(a+b)`,
`2x+y<=4a`, `y+2z<=4b`, and `n_ee>=0`. The lattice adjacency cap is toroidal, including the `|S|=0` and `|S|=P` cases.

At a fixed **high** source lump `(a,b)` the trial-vector value `d(a,b)>0` is fixed. Set `D(a',b')` to the high-lump trial-vector value if the target is high, or the **largest true low-orbit trial-vector value** among all states of that low lump. The flat target has `D(0,0)=0`. Then each electric pattern has a certified lower contribution

```
weight(pattern) = sum_over_3_raw_S3_moves [1 - D(target_lump)/d(source_lump)]
```

because for every actual destination `d(destination)<=D(target_lump)`. Hence the row ratio is bounded below by an affine expression in `x,y,z,u`; the coefficient of `u` is linear, so its minimum is at `u=0` or `u=x`. The remaining optimization is a **three-variable rational linear program**. A minimum over all polytope vertices is at most the minimum over every actual integer signature, so it is a valid lower bound. Extra relaxed (unphysical) vertices may make the lower bound conservative; they do not invalidate it.

The generated independent S₃ table handles the special T–T transition to the flat state at lump `(2,0)`; every flat target contributes zero to the source ratio.

## Acceptance and negative controls

`gpt6_extreme_row_oracle.py` generates (S_3) from all six permutations and the 3 transpositions (not from the model table), constructs a different cross-product vertex solver, enumerates all integer signatures on 4×4, and verifies the vertex minimum is no larger, at every nonflat lump, for both a product vector and a synthetic **hybrid** vector with nonuniform low-orbit values. It checks the square 5×5 hybrid vertices (and the full integer enumeration with `--full`) and 6×6 vertices.

Two controls **must turn red** if an implementation incorrectly restricts allowed vertices or understates target weights:

1. Replace low-lump `max d` by `min d`. An actual transition to a higher-weight low-orbit state yields a smaller `(H'd)/d` than the purported lower bound.
2. Impose the false condition `x=0` on a T-rich lump. An admissible integer signature beats the claimed 'extreme-row' lower bound.

An independent product-vector calibration at `L=4,r=12` requires a gap-ratio minimum of `657/32=20.53125`, the same one-3-cycle physical witness already verified in the linear-law review.

## Not yet accepted

**G1–G5 on the hybrid Rust implementation remain pending.** In particular this oracle neither enumerates all gauge orbits below the energy cut, nor checks their group multiplicities, full mixed-orbit source rows, an actual certified (E_9-E_8) value, or the 7×7/8×8 large-key path. When the independent Rust program is available, its *own* frozen high-lump coefficients and all low-lump maxima should be exported and checked against this separate oracle at G1 and G2 before accepting a large-torus result.

This preserves the proof distinction: the finite-torus hybrid lower bound is conditional on complete source support and exact extremal rows, while a volume-independent stability threshold is a different theorem.
