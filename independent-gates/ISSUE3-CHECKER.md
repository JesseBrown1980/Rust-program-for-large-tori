# Issue #3: independent checker of HYBRID_BIG dumped trial vectors

This is a second **standalone Rust 1.81 binary**:
\`src/bin/vectorcheck.rs\`. It **does not link to or import
\`src/main.rs\`**. Its S₃ multiplication table, 2D torus, geometric
symmetries, spanning-tree gauge, gauge-fixed canonical class, orbit
size, exact electric transition rows, and above-cut linear optimizer
are independently implemented.

## Frozen producer contract

The producer's \`dump=FILE\` writes

    HYBRID_BIG|Lx=N|Ly=N|r=N|cut=N|classes=N|lumps_above=N|claim_micro=N|key=free_link_digits_of_the_comb_gauge|json=0
    P|digits|a|b|weight|value
    G|a|b|value

Every \`P\` record is a nonflat class with \`0<M=2a+3b<=cut\`,
a canonical digit key, an orbit-size weight and positive integer trial
value. Every \`G\` record is a positive trial value for one complete
above-cut lump. No secret state files are required. The checker
refuses absent/duplicate classes, malformed keys, nonpositive values,
incorrect lump or orbit weights, missing high lumps, and numbers
exceeding the source's TOP=2^86 safe-value bound.

## What the checker independently proves

**1. Low-cut census.** Every \`P\` key decodes into the stipulated
comb gauge, and its plaquette classes are recomputed from four
ordered S₃ link products. Independent gauge fixing, translations,
reflections, axis swap on squares and global conjugation produce
the canonical address and orbit weight. The per-lump sum must equal

    C(P;a,b) * (12*3^a*2^b)              if a>=2, even,
    C(P;0,b) * (12*2^b+6*(-1)^b)         if a=0.

**2. Every low-cut electric row.** All \`3 * (2 P)\` raw moves are
generated from the three S₃ transpositions at each link. Every
destination's plaquettes are recomputed, then a below-cut target
is resolved with the independent canonicalizer against the
dumped \`P\` map, while a higher target uses the exact dumped
\`G\` value. Flat targets have zero trial value. The whole exact
integer row ratio is computed with an **outward** ceiling on
the nonnegative electric sum:

    row_lower_micro =
      [3*(2P)+r*(2a+3b)]*1e6 -
      ceil(1e6*sum_of_all_electric_target_values / P_source).

Any row below \`claim_micro\` is a refusal.

**3. Every above-cut lump.** The independent 3D rational polytope
has all nine proof constraints, including nonnegative counts for
every link pattern. Its vertices are enumerated by a separate
cross-product intersection solver. The seven electric link
patterns are freshly generated from S₃ and all three raw moves.
For a move to a low-cut target *lump*, the checker uses the
**maximum actual P value over that lump**. T–T same/different
outgoing sums are maximized independently at their affine
endpoints. The maximum incoming target sum at the extreme
vertices is rounded **up**, producing a conservative floor.

The value \`TOP\` is used only as a conservative penalty for a
pattern transition into an impossible or missing lump; the checker
independently requires **every actual dumped vector value at
most TOP**. The resulting penalty cannot *increase* a certified
lower bound. Every feasible above-cut lump must be represented
in \`G\`, with no extras.

**4. Source-claim acceptance.** The global minimum of all real
low-cut rows and all high-lump vertex rows must equal the
dumped integer-micro claim. A claimed value increased by even
one micro-unit fails. An independent 4-shard fan-in verifies
all source classes are covered exactly once and checks that
the minimum does not admit the +1 claim.

**5. Planted source faults.** A separate hosted test modifies a
real loaded dump by deleting a class, corrupting its magnetic
lump, or removing one raw electric move. Each action must cause
a typed refusal; no broken record may be promoted.

## Performance and scope

For large dumps, run disjoint low-cut source shards:

    vectorcheck VECTOR.dump --shard=0/4
    vectorcheck VECTOR.dump --shard=1/4
    vectorcheck VECTOR.dump --shard=2/4
    vectorcheck VECTOR.dump --shard=3/4

Each row includes \`scope=SHARD_ONLY\`, not a full certificate.
All four must pass against the **same frozen input**; the
\`independent-gates/vectorcheck_fanin.py\` gate sums exact
class/move counts, compares high-source outputs and asserts
that the full lower bound equals the claim.

**Acceptance progression:** first generate and verify an actual
4×4, r=9, cut10 dump on GitHub. After it passes, feed the
completed 6×6 dump, then 7×7 and 8×8. Larger jobs may need
more source shards to stay within memory/time limits, plus
their own fan-in configuration. Source sharding is a resource
strategy, never permission to ignore rows.

This checks a Barta certificate **conditioned on the original
finite-volume \(S_3\) spectral setup**, not the separate continuum
gap or a size-independent stability constant. A pass on 4×4
is not yet an accepted 7×7 or 8×8 certificate.

## Independent review points in producer

- \`certificate()\` rounds \`ceil(sum*1e6/p)\` outward before
  subtraction and therefore cannot overstate a lower bound
  if \`sum\` and \`p\` are exactly represented and positive.
- The negative-row search at high lumps uses maximal outgoing
  neighbor values (including low-lump maxima), not minimal
  values, and maxima over the relaxed vertex domain.
- \`canon()\` computes orbit weight by separating the number
  of spatial pattern placements from the number of distinct
  canonical gauge-fixed/conjugated link arrays. This checker
  recomputes the same finite orbit from an independently built
  lattice symmetry action and compares its weight per class.
- The \`TOP\` fallback is safe only if **all real source vector
  entries are at most TOP**, there are no signed/overflowed
  arithmetic paths, and the unphysical-pattern cost does not
  reduce the outgoing upper bound. The independent checker
  enforces the first restriction and uses checked i128 arithmetic.

All tests remain separate from the producing Rust source. This
draft is an external checker in the same repository and requires
independent review before the operator merges it.
