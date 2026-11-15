//! Hybrid no-store certificate for the S3 lattice gauge Hamiltonian H = J_E (E + r M) on a periodic Lx x Ly torus,
//! written for tori too large for the 128-bit keys and the full row table of the earlier programs.
//!
//!   floor  -> centre : lambda_1(H') >= min over rows of (H'd)/d for a positive d (Barta), so E9 >= that minimum;
//!   ceiling-> centre : lambda_1(H') <= Rayleigh quotient of the same vector restricted to the exact region;
//!   centre outward   : the exact region is every gauge orbit with 0 < M <= CUT, grown from the eight flat states.
//!
//! d is one value per class (gauge, global conjugation, lattice symmetries) with 0 < M <= CUT, one value per lump
//! (a, b) above CUT, and 0 on the flat states; a, b count plaquettes whose holonomy is a transposition, a 3-cycle.
//! Rows at M <= CUT are exact, from the definition. Rows above CUT use the vertices of each lump's polytope of
//! link-pattern counts (four links per plaquette, the adjacency caps): the ratio is linear in those counts, so its
//! minimum over the polytope is at a vertex and bounds the minimum over the integer rows from below.
//!
//! usage: large_tori LX LY R CUT THREADS SWEEPS TARGETS_MICRO[,..]|auto [REPLAY_R[,..]]
//! Integers only. Keys are digit arrays, so no packing limit. Output rows end in json=0.

use std::collections::HashMap;
use std::env;
use std::thread;

const U: i128 = 1_000_000;
const TOP: i128 = 1 << 86;
const PERMS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

struct Group {
    mul: [[u8; 6]; 6],
    inv: [u8; 6],
    cls: [u8; 6], // 0 identity, 1 transposition, 2 three-cycle
    trans: [u8; 3],
    conj: [[u8; 6]; 6],
}

impl Group {
    fn new() -> Self {
        let mut mul = [[0u8; 6]; 6];
        for a in 0..6 {
            for (b, pb) in PERMS.iter().enumerate() {
                let c = [PERMS[a][pb[0]], PERMS[a][pb[1]], PERMS[a][pb[2]]];
                mul[a][b] = PERMS.iter().position(|p| *p == c).unwrap_or(0) as u8;
            }
        }
        let mut inv = [0u8; 6];
        for a in 0..6 {
            inv[a] = (0..6).find(|&b| mul[a][b] == 0).unwrap_or(0) as u8;
        }
        let mut cls = [0u8; 6];
        for (a, slot) in cls.iter_mut().enumerate() {
            let fixed = (0..3).filter(|&j| PERMS[a][j] == j).count();
            *slot = match fixed {
                3 => 0,
                1 => 1,
                _ => 2,
            };
        }
        let t: Vec<u8> = (0..6u8).filter(|&a| cls[a as usize] == 1).collect();
        let mut conj = [[0u8; 6]; 6];
        for c in 0..6 {
            for x in 0..6 {
                conj[c][x] = mul[mul[c][x] as usize][inv[c] as usize];
            }
        }
        Group { mul, inv, cls, trans: [t[0], t[1], t[2]], conj }
    }
    fn m(&self, a: u8, b: u8) -> u8 {
        self.mul[a as usize][b as usize]
    }
    fn i(&self, a: u8) -> u8 {
        self.inv[a as usize]
    }
}

struct Aut {
    link: Vec<(u32, bool)>, // image link, orientation reversed
    plaq: Vec<u32>,
}

struct Torus {
    lx: usize,
    ly: usize,
    p: usize,
    nl: usize,
    free: Vec<u32>,
    tree: Vec<(u32, u32, u32)>, // (from vertex, to vertex, link), in an order that fixes the gauge from the root
    ends: Vec<(u32, u32)>,
    plaq_of: Vec<[u32; 2]>,
    auts: Vec<Aut>,
}

impl Torus {
    fn vid(&self, x: usize, y: usize) -> usize {
        (x % self.lx) * self.ly + (y % self.ly)
    }
    fn h(&self, x: usize, y: usize) -> usize {
        2 * self.vid(x, y)
    }
    fn v(&self, x: usize, y: usize) -> usize {
        2 * self.vid(x, y) + 1
    }

    fn new(lx: usize, ly: usize) -> Self {
        let mut t = Torus { lx, ly, p: lx * ly, nl: 2 * lx * ly, free: Vec::new(), tree: Vec::new(), ends: Vec::new(), plaq_of: Vec::new(), auts: Vec::new() };
        let mut in_tree = vec![false; t.nl];
        for y in 0..ly - 1 {
            t.tree.push((t.vid(0, y) as u32, t.vid(0, y + 1) as u32, t.v(0, y) as u32));
            in_tree[t.v(0, y)] = true;
        }
        for x in 0..lx - 1 {
            for y in 0..ly {
                t.tree.push((t.vid(x, y) as u32, t.vid(x + 1, y) as u32, t.h(x, y) as u32));
                in_tree[t.h(x, y)] = true;
            }
        }
        t.free = (0..t.nl as u32).filter(|&l| !in_tree[l as usize]).collect();
        t.ends = vec![(0, 0); t.nl];
        t.plaq_of = vec![[0, 0]; t.nl];
        for x in 0..lx {
            for y in 0..ly {
                let (hl, vl, here) = (t.h(x, y), t.v(x, y), t.vid(x, y) as u32);
                let (right, up, below, left) = (t.vid(x + 1, y) as u32, t.vid(x, y + 1) as u32, t.vid(x, y + ly - 1) as u32, t.vid(x + lx - 1, y) as u32);
                t.ends[hl] = (here, right);
                t.ends[vl] = (here, up);
                t.plaq_of[hl] = [here, below];
                t.plaq_of[vl] = [here, left];
            }
        }
        // lattice symmetries: translations, the two reflections, and the transposition of the axes on a square
        let mut link_of: HashMap<(u32, u32), (u32, bool)> = HashMap::new();
        for (l, &(a, b)) in t.ends.iter().enumerate() {
            link_of.insert((a, b), (l as u32, false));
            link_of.insert((b, a), (l as u32, true));
        }
        let swaps: &[bool] = if lx == ly { &[false, true] } else { &[false] };
        for &swap in swaps {
            for sx in [false, true] {
                for sy in [false, true] {
                    for tx in 0..lx {
                        for ty in 0..ly {
                            let map = |x: usize, y: usize| -> (usize, usize) {
                                let (x, y) = if swap { (y, x) } else { (x, y) };
                                let x2 = if sx { (lx - x % lx) % lx } else { x % lx };
                                let y2 = if sy { (ly - y % ly) % ly } else { y % ly };
                                ((x2 + tx) % lx, (y2 + ty) % ly)
                            };
                            let vmap = |v: u32| -> u32 {
                                let (x, y) = (v as usize / ly, v as usize % ly);
                                let (a, b) = map(x, y);
                                t.vid(a, b) as u32
                            };
                            let link: Vec<(u32, bool)> = t.ends.iter().map(|&(a, b)| link_of[&(vmap(a), vmap(b))]).collect();
                            let mut plaq = vec![0u32; t.p];
                            for x in 0..lx {
                                for y in 0..ly {
                                    let corners = [map(x, y), map(x + 1, y), map(x, y + 1), map(x + 1, y + 1)];
                                    let base = corners.iter().find(|&&(cx, cy)| {
                                        corners.contains(&((cx + 1) % lx, cy)) && corners.contains(&(cx, (cy + 1) % ly)) && corners.contains(&((cx + 1) % lx, (cy + 1) % ly))
                                    });
                                    let (bx, by) = *base.expect("a plaquette maps to a plaquette");
                                    plaq[t.vid(x, y)] = t.vid(bx, by) as u32;
                                }
                            }
                            t.auts.push(Aut { link, plaq });
                        }
                    }
                }
            }
        }
        t
    }

    fn hol(&self, g: &Group, u: &[u8], p: usize) -> u8 {
        let (x, y) = (p / self.ly, p % self.ly);
        let a = g.m(u[self.h(x, y)], u[self.v(x + 1, y)]);
        g.m(g.m(a, g.i(u[self.h(x, y + 1)])), g.i(u[self.v(x, y)]))
    }

    /// Digits of the gauge-fixed configuration: the tree links become the identity.
    fn gauge_fix(&self, g: &Group, u: &[u8], gv: &mut [u8], out: &mut Vec<u8>) {
        gv[0] = 0;
        for &(a, b, l) in &self.tree {
            gv[b as usize] = g.m(gv[a as usize], u[l as usize]);
        }
        out.clear();
        for &l in &self.free {
            let (a, b) = self.ends[l as usize];
            out.push(g.m(g.m(gv[a as usize], u[l as usize]), g.i(gv[b as usize])));
        }
    }

    fn decode(&self, key: &[u8], u: &mut [u8]) {
        u.fill(0);
        for (d, &l) in key.iter().zip(&self.free) {
            u[l as usize] = *d;
        }
    }
}

struct Scratch {
    cls: Vec<u8>,
    defects: Vec<(u32, u8)>,
    pat: Vec<(u32, u8)>,
    best_pat: Vec<(u32, u8)>,
    cand: Vec<u32>,
    img: Vec<u8>,
    gv: Vec<u8>,
    digits: Vec<u8>,
    trial: Vec<u8>,
    all: Vec<Vec<u8>>,
}

impl Scratch {
    fn new(t: &Torus) -> Self {
        Scratch {
            cls: vec![0; t.p],
            defects: Vec::new(),
            pat: Vec::new(),
            best_pat: Vec::new(),
            cand: Vec::new(),
            img: vec![0; t.nl],
            gv: vec![0; t.p],
            digits: Vec::new(),
            trial: Vec::new(),
            all: Vec::new(),
        }
    }
}

/// Canonical key of a configuration: the least gauge-fixed digit array over the symmetries that give the least
/// pattern of defects, and over the six global conjugations. With `weigh`, also the number of gauge-fixed
/// configurations in its class.
fn canon(t: &Torus, g: &Group, u: &[u8], sc: &mut Scratch, weigh: bool) -> (Vec<u8>, u64) {
    sc.defects.clear();
    for p in 0..t.p {
        let c = g.cls[t.hol(g, u, p) as usize];
        sc.cls[p] = c;
        if c != 0 {
            sc.defects.push((p as u32, c));
        }
    }
    sc.cand.clear();
    let mut first = true;
    for (i, a) in t.auts.iter().enumerate() {
        sc.pat.clear();
        for &(p, c) in &sc.defects {
            sc.pat.push((a.plaq[p as usize], c));
        }
        sc.pat.sort_unstable();
        if first || sc.pat < sc.best_pat {
            first = false;
            std::mem::swap(&mut sc.pat, &mut sc.best_pat);
            sc.cand.clear();
            sc.cand.push(i as u32);
        } else if sc.pat == sc.best_pat {
            sc.cand.push(i as u32);
        }
    }
    let mut best: Vec<u8> = Vec::new();
    sc.all.clear();
    let cand = std::mem::take(&mut sc.cand);
    for &i in &cand {
        let a = &t.auts[i as usize];
        for (l, &(to, rev)) in a.link.iter().enumerate() {
            sc.img[to as usize] = if rev { g.i(u[l]) } else { u[l] };
        }
        let img = std::mem::take(&mut sc.img);
        let mut digits = std::mem::take(&mut sc.digits);
        t.gauge_fix(g, &img, &mut sc.gv, &mut digits);
        sc.img = img;
        for c in 0..6 {
            sc.trial.clear();
            sc.trial.extend(digits.iter().map(|&d| g.conj[c][d as usize]));
            if best.is_empty() || sc.trial < best {
                best.clone_from(&sc.trial);
            }
            if weigh {
                sc.all.push(sc.trial.clone());
            }
        }
        sc.digits = digits;
    }
    let weight = if weigh {
        sc.all.sort_unstable();
        sc.all.dedup();
        (t.auts.len() / cand.len()) as u64 * sc.all.len() as u64
    } else {
        0
    };
    sc.cand = cand;
    (best, weight)
}

fn lump_of(t: &Torus, g: &Group, u: &[u8]) -> (usize, usize) {
    let (mut a, mut b) = (0, 0);
    for p in 0..t.p {
        match g.cls[t.hol(g, u, p) as usize] {
            1 => a += 1,
            2 => b += 1,
            _ => {}
        }
    }
    (a, b)
}

/// The lump after link l is multiplied by transposition s, given the lump before.
fn moved_lump(t: &Torus, g: &Group, u: &mut [u8], l: usize, s: u8, a: usize, b: usize) -> (usize, usize) {
    let [p, q] = t.plaq_of[l];
    let before = [g.cls[t.hol(g, u, p as usize) as usize], g.cls[t.hol(g, u, q as usize) as usize]];
    let old = u[l];
    u[l] = g.m(s, old);
    let after = [g.cls[t.hol(g, u, p as usize) as usize], g.cls[t.hol(g, u, q as usize) as usize]];
    u[l] = old;
    let (mut a2, mut b2) = (a as i64, b as i64);
    for (c, sign) in before.iter().map(|c| (*c, -1i64)).chain(after.iter().map(|c| (*c, 1i64))) {
        match c {
            1 => a2 += sign,
            2 => b2 += sign,
            _ => {}
        }
    }
    (a2 as usize, b2 as usize)
}

fn ceil2sqrt(n: usize) -> usize {
    let mut v = 0;
    while v * v < 4 * n {
        v += 1;
    }
    v
}

fn adj(s: usize, t: &Torus) -> i64 {
    if s == 0 {
        0
    } else if s == t.p {
        2 * s as i64
    } else {
        2 * s as i64 - ceil2sqrt(s).min(t.lx.min(t.ly)).min(ceil2sqrt(t.p - s)) as i64
    }
}

/// Vertices of the polytope of (x, y, z) = numbers of T-T, T-C, C-C links of a lump, as (det, X, Y, Z), det > 0.
fn lump_vertices(a: usize, b: usize, t: &Torus) -> Vec<[i64; 4]> {
    let (a4, b4, nl) = (4 * a as i64, 4 * b as i64, t.nl as i64);
    let rows: [[i64; 4]; 9] = [
        [-1, 0, 0, 0],
        [0, -1, 0, 0],
        [0, 0, -1, 0],
        [1, 0, 0, adj(a, t)],
        [0, 0, 1, adj(b, t)],
        [1, 1, 1, adj(a + b, t)],
        [2, 1, 0, a4],
        [0, 1, 2, b4],
        [-1, -1, -1, nl - a4 - b4],
    ];
    let mut out: Vec<[i64; 4]> = Vec::new();
    for i in 0..9 {
        for j in i + 1..9 {
            for k in j + 1..9 {
                let (r1, r2, r3) = (rows[i], rows[j], rows[k]);
                let det = r1[0] * (r2[1] * r3[2] - r2[2] * r3[1]) - r1[1] * (r2[0] * r3[2] - r2[2] * r3[0]) + r1[2] * (r2[0] * r3[1] - r2[1] * r3[0]);
                if det == 0 {
                    continue;
                }
                let dx = r1[3] * (r2[1] * r3[2] - r2[2] * r3[1]) - r1[1] * (r2[3] * r3[2] - r2[2] * r3[3]) + r1[2] * (r2[3] * r3[1] - r2[1] * r3[3]);
                let dy = r1[0] * (r2[3] * r3[2] - r2[2] * r3[3]) - r1[3] * (r2[0] * r3[2] - r2[2] * r3[0]) + r1[2] * (r2[0] * r3[3] - r2[3] * r3[0]);
                let dz = r1[0] * (r2[1] * r3[3] - r2[3] * r3[1]) - r1[1] * (r2[0] * r3[3] - r2[3] * r3[0]) + r1[3] * (r2[0] * r3[1] - r2[1] * r3[0]);
                let s = if det < 0 { -1 } else { 1 };
                let v = [det * s, dx * s, dy * s, dz * s];
                if rows.iter().all(|c| c[0] * v[1] + c[1] * v[2] + c[2] * v[3] <= c[3] * v[0]) && !out.contains(&v) {
                    out.push(v);
                }
            }
        }
    }
    out
}

struct State {
    key: Vec<u8>,
    a: u16,
    b: u16,
    weight: u64,
    exact: Vec<u32>,        // targets with 0 < M <= cut, with multiplicity
    above: Vec<(u32, u16)>, // (lump index, moves) above the cut
    to_flat: u16,
}

struct Model {
    side: usize,
    cut: usize,
    states: Vec<State>,
    lumps: Vec<(usize, usize, Vec<[i64; 4]>)>, // lumps above the cut that have a non-empty polytope
    lump_at: Vec<u32>,                         // (a, b) -> index into lumps, or u32::MAX
}

fn cost(a: usize, b: usize) -> usize {
    2 * a + 3 * b
}

fn build(t: &Torus, g: &Group, cut: usize, threads: usize) -> Result<Model, String> {
    let side = t.p + 1;
    // lumps above the cut
    let mut lumps = Vec::new();
    let mut lump_at = vec![u32::MAX; side * side];
    for a in (0..=t.p).step_by(2) {
        for b in 0..=(t.p - a) {
            if cost(a, b) > cut {
                let v = lump_vertices(a, b, t);
                if !v.is_empty() {
                    lump_at[a * side + b] = lumps.len() as u32;
                    lumps.push((a, b, v));
                }
            }
        }
    }
    // the 18 flat gauge-fixed configurations: commuting pairs of winding holonomies
    let mut seeds: Vec<Vec<u8>> = Vec::new();
    let mut flat_configs = 0u64;
    let mut sc = Scratch::new(t);
    let mut u = vec![0u8; t.nl];
    for ha in 0..6u8 {
        for hb in 0..6u8 {
            if g.m(ha, hb) != g.m(hb, ha) {
                continue;
            }
            u.fill(0);
            for y in 0..t.ly {
                u[t.h(t.lx - 1, y)] = ha;
            }
            for x in 0..t.lx {
                u[t.v(x, t.ly - 1)] = hb;
            }
            if lump_of(t, g, &u) != (0, 0) {
                return Err("a seed is not flat".into());
            }
            let (k, w) = canon(t, g, &u, &mut sc, true);
            if !seeds.contains(&k) {
                seeds.push(k);
                flat_configs += w;
            }
        }
    }
    // 8 flat gauge orbits, 18 gauge-fixed configurations; the lattice symmetries may merge some of the 8
    if flat_configs != 18 {
        return Err(format!("{flat_configs} flat configurations in {} classes, expected 18", seeds.len()));
    }
    // classes with 0 < M <= cut, grown outward from the flat classes
    let mut index: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut keys: Vec<Vec<u8>> = Vec::new();
    let mut frontier = seeds;
    while !frontier.is_empty() {
        let chunk = frontier.len().div_ceil(threads.max(1)).max(1);
        let found: Vec<Vec<Vec<u8>>> = thread::scope(|s| {
            let handles: Vec<_> = frontier
                .chunks(chunk)
                .map(|part| {
                    s.spawn(move || {
                        let mut sc = Scratch::new(t);
                        let mut u = vec![0u8; t.nl];
                        let mut out: Vec<Vec<u8>> = Vec::new();
                        for key in part {
                            t.decode(key, &mut u);
                            let (a, b) = lump_of(t, g, &u);
                            for l in 0..t.nl {
                                for &s in &g.trans {
                                    let (a2, b2) = moved_lump(t, g, &mut u, l, s, a, b);
                                    let m2 = cost(a2, b2);
                                    if m2 == 0 || m2 > cut {
                                        continue;
                                    }
                                    let old = u[l];
                                    u[l] = g.m(s, old);
                                    out.push(canon(t, g, &u, &mut sc, false).0);
                                    u[l] = old;
                                }
                            }
                        }
                        out.sort_unstable();
                        out.dedup();
                        out
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap_or_default()).collect()
        });
        let mut next = Vec::new();
        for part in found {
            for k in part {
                if !index.contains_key(&k) {
                    index.insert(k.clone(), keys.len() as u32);
                    keys.push(k.clone());
                    next.push(k);
                }
            }
        }
        frontier = next;
    }
    // rows of every class, from the definition
    let n = keys.len();
    let chunk = n.div_ceil(threads.max(1)).max(1);
    let (index, lump_ref) = (&index, &lump_at);
    let parts: Vec<Result<Vec<State>, String>> = thread::scope(|s| {
        let handles: Vec<_> = keys
            .chunks(chunk)
            .map(|part| {
                s.spawn(move || {
                    let mut sc = Scratch::new(t);
                    let mut u = vec![0u8; t.nl];
                    let mut out = Vec::with_capacity(part.len());
                    for key in part {
                        t.decode(key, &mut u);
                        let (a, b) = lump_of(t, g, &u);
                        let (own, weight) = canon(t, g, &u, &mut sc, true);
                        if own != *key {
                            return Err("a key is not its own canonical form".to_string());
                        }
                        let mut st = State { key: key.clone(), a: a as u16, b: b as u16, weight, exact: Vec::new(), above: Vec::new(), to_flat: 0 };
                        let mut above: HashMap<u32, u16> = HashMap::new();
                        for l in 0..t.nl {
                            for &s in &g.trans {
                                let (a2, b2) = moved_lump(t, g, &mut u, l, s, a, b);
                                let m2 = cost(a2, b2);
                                if m2 == 0 {
                                    st.to_flat += 1;
                                } else if m2 > cut {
                                    let li = lump_ref[a2 * side + b2];
                                    if li == u32::MAX {
                                        return Err(format!("a move reaches lump ({a2},{b2}) whose polytope is empty"));
                                    }
                                    *above.entry(li).or_insert(0) += 1;
                                } else {
                                    let old = u[l];
                                    u[l] = g.m(s, old);
                                    let k = canon(t, g, &u, &mut sc, false).0;
                                    u[l] = old;
                                    match index.get(&k) {
                                        Some(&j) => st.exact.push(j),
                                        None => return Err("a move inside the cut reaches a class that the outward growth did not find".to_string()),
                                    }
                                }
                            }
                        }
                        st.above = above.into_iter().collect();
                        st.above.sort_unstable();
                        out.push(st);
                    }
                    Ok(out)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap_or_else(|_| Err("a row thread stopped".into()))).collect()
    });
    let mut states = Vec::with_capacity(n);
    for part in parts {
        states.extend(part?);
    }
    Ok(Model { side, cut, states, lumps, lump_at })
}

fn binomial(n: usize, k: usize) -> u128 {
    let mut out = 1u128;
    for i in 0..k {
        out = out * (n - i) as u128 / (i + 1) as u128;
    }
    out
}

/// Gauge-fixed configurations of one lump: C(P; a, b) * N(a, b), N = 12 * 3^a * 2^b for even a >= 2 and 12 * 2^b + 6 (-1)^b for a = 0.
fn expected(p: usize, a: usize, b: usize) -> u128 {
    let places = binomial(p, a) * binomial(p - a, b);
    let n = if a == 0 {
        let base = 12 * (1u128 << b);
        if b % 2 == 0 {
            base + 6
        } else {
            base - 6
        }
    } else {
        12 * 3u128.pow(a as u32) * (1u128 << b)
    };
    places * n
}

/// Sum of the neighbour values of each link pattern for a lump above the cut: [ee, eT, eC, TT, TC, CC].
fn pattern_sums(md: &Model, t: &Torus, a: usize, b: usize, g: &[i128], dmax: &[i128]) -> [i128; 6] {
    let val = |da: i64, db: i64| -> i128 {
        let (a2, b2) = (a as i64 + da, b as i64 + db);
        if a2 < 0 || b2 < 0 || (a2 + b2) as usize > t.p {
            return TOP; // a move that cannot occur; the largest value keeps the bound safe
        }
        let (a2, b2) = (a2 as usize, b2 as usize);
        if cost(a2, b2) == 0 {
            0
        } else if cost(a2, b2) <= md.cut {
            dmax[a2 * md.side + b2]
        } else {
            match md.lump_at[a2 * md.side + b2] {
                u32::MAX => TOP,
                li => g[li as usize],
            }
        }
    };
    let same = val(-2, 0) + 2 * val(-2, 2);
    let diff = 2 * val(-2, 1) + val(-2, 2);
    [3 * val(2, 0), val(0, 0) + 2 * val(0, 1), 3 * val(2, -1), same.max(diff), val(0, -1) + 2 * val(0, 0), 3 * val(2, -2)]
}

/// Largest sum of neighbour values over the vertices of a lump, as (numerator, det): the row's ratio is diag - num / (det * g).
fn lump_worst(md: &Model, t: &Torus, li: usize, g: &[i128], dmax: &[i128]) -> (i128, i128) {
    let (a, b, verts) = &md.lumps[li];
    let s = pattern_sums(md, t, *a, *b, g, dmax);
    let (a4, b4, nl) = (4 * *a as i128, 4 * *b as i128, t.nl as i128);
    let mut best: (i128, i128) = (0, 1);
    for v in verts {
        let (det, x, y, z) = (i128::from(v[0]), i128::from(v[1]), i128::from(v[2]), i128::from(v[3]));
        let n = [det * (nl - a4 - b4) + x + y + z, det * a4 - 2 * x - y, det * b4 - y - 2 * z, x, y, z];
        let num: i128 = n.iter().zip(&s).map(|(c, w)| c * w).sum();
        if num * best.1 > best.0 * det {
            best = (num, det);
        }
    }
    best
}

struct Fit {
    psi: Vec<i128>,
    g: Vec<i128>,
    cert: i128,
    min_exact: i128,
    min_lump: i128,
}

fn certificate(md: &Model, t: &Torus, r: i128, psi: &[i128], g: &[i128]) -> (i128, i128, Vec<i128>) {
    let nl3 = 3 * t.nl as i128;
    let mut dmax = vec![0i128; md.side * md.side];
    for (st, p) in md.states.iter().zip(psi) {
        let slot = &mut dmax[st.a as usize * md.side + st.b as usize];
        *slot = (*slot).max(*p);
    }
    // per M: the least of (ratio - r M), so that the same vector can be read at another coupling
    let mut electric = vec![i128::MAX; 3 * t.p + 1];
    let mut min_exact = i128::MAX;
    for (st, p) in md.states.iter().zip(psi) {
        let m = cost(st.a as usize, st.b as usize) as i128;
        let s: i128 = st.exact.iter().map(|&j| psi[j as usize]).sum::<i128>() + st.above.iter().map(|&(li, c)| i128::from(c) * g[li as usize]).sum::<i128>();
        let el = nl3 * U - (s * U + p - 1) / p; // 3 nl - s / p, rounded down
        electric[m as usize] = electric[m as usize].min(el);
        min_exact = min_exact.min(el + r * m * U);
    }
    let mut min_lump = i128::MAX;
    for (li, (a, b, _)) in md.lumps.iter().enumerate() {
        let m = cost(*a, *b) as i128;
        let (num, det) = lump_worst(md, t, li, g, &dmax);
        let den = det * g[li];
        let el = nl3 * U - (num * U + den - 1) / den;
        electric[m as usize] = electric[m as usize].min(el);
        min_lump = min_lump.min(el + r * m * U);
    }
    (min_exact, min_lump, electric)
}

fn fit(md: &Model, t: &Torus, r: i128, target: i128, sweeps: usize, threads: usize) -> Option<Fit> {
    let nl3 = 3 * t.nl as i128;
    let n = md.states.len();
    let diag: Vec<i128> = md.states.iter().map(|s| (r * cost(s.a as usize, s.b as usize) as i128 + nl3) * U - target).collect();
    if diag.iter().any(|d| *d <= 0) {
        return None;
    }
    let lump_diag: Vec<i128> = md.lumps.iter().map(|(a, b, _)| (r * cost(*a, *b) as i128 + nl3) * U - target).collect();
    if lump_diag.iter().any(|d| *d <= 0) {
        return None;
    }
    let mut psi = vec![1i128 << 60; n];
    let mut g = vec![1i128; md.lumps.len()];
    let mut dmax = vec![0i128; md.side * md.side];
    let chunk = n.div_ceil(threads.max(1)).max(1);
    for _ in 0..sweeps {
        dmax.fill(0);
        for (st, p) in md.states.iter().zip(&psi) {
            let slot = &mut dmax[st.a as usize * md.side + st.b as usize];
            *slot = (*slot).max(*p);
        }
        // lumps above the cut, three passes from the lowest cost upward
        for _ in 0..3 {
            for li in 0..md.lumps.len() {
                let (num, det) = lump_worst(md, t, li, &g, &dmax);
                let need = (num * U + det * lump_diag[li] - 1) / (det * lump_diag[li]);
                if need > TOP {
                    return None;
                }
                g[li] = need.max(1);
            }
        }
        let mut next = vec![0i128; n];
        let (psi_ref, g_ref, diag_ref) = (&psi, &g, &diag);
        thread::scope(|s| {
            for (ci, out) in next.chunks_mut(chunk).enumerate() {
                s.spawn(move || {
                    for (k, slot) in out.iter_mut().enumerate() {
                        let i = ci * chunk + k;
                        let st = &md.states[i];
                        let sum: i128 = st.exact.iter().map(|&j| psi_ref[j as usize]).sum::<i128>() + st.above.iter().map(|&(li, c)| i128::from(c) * g_ref[li as usize]).sum::<i128>();
                        *slot = ((psi_ref[i] + sum * U / diag_ref[i]) / 2).max(1);
                    }
                });
            }
        });
        let top = next.iter().copied().max().unwrap_or(1).max(1);
        if top > TOP {
            let shift = 128 - (top / TOP).leading_zeros();
            for x in next.iter_mut() {
                *x = (*x >> shift).max(1);
            }
            for x in g.iter_mut() {
                *x = (*x >> shift).max(1);
            }
        } else if top < TOP >> 8 {
            let shift = (TOP / top).leading_zeros().abs_diff(127).saturating_sub(1).min(60);
            for x in next.iter_mut() {
                *x <<= shift;
            }
            for x in g.iter_mut() {
                *x = (*x << shift).min(TOP);
            }
        }
        psi = next;
    }
    let (min_exact, min_lump, _) = certificate(md, t, r, &psi, &g);
    Some(Fit { psi, g, cert: min_exact.min(min_lump), min_exact, min_lump })
}

/// Rayleigh quotient of the vector restricted to the exact region, rounded up: an upper end for lambda_1(H').
fn rayleigh(md: &Model, t: &Torus, r: i128, psi: &[i128]) -> i128 {
    let nl3 = 3 * t.nl as i128;
    let top = psi.iter().copied().max().unwrap_or(1);
    let shift = (128 - top.leading_zeros()).saturating_sub(34);
    let phi: Vec<i128> = psi.iter().map(|p| (*p >> shift).max(1)).collect();
    let (mut num, mut den) = (0i128, 0i128);
    for (i, st) in md.states.iter().enumerate() {
        let d = r * cost(st.a as usize, st.b as usize) as i128 + nl3;
        let s: i128 = st.exact.iter().map(|&j| phi[j as usize]).sum();
        num += i128::from(st.weight) * phi[i] * (d * phi[i] - s);
        den += i128::from(st.weight) * phi[i] * phi[i];
    }
    let (q, rem) = (num.div_euclid(den), num.rem_euclid(den));
    q * U + ((rem as u128 * U as u128) / den as u128) as i128 + 1
}

fn isqrt(n: i128) -> i128 {
    // integer Newton from above
    let mut x = 1i128 << ((128 - n.leading_zeros()) / 2 + 1);
    loop {
        let y = (x + n / x) / 2;
        if y >= x {
            return x;
        }
        x = y;
    }
}

fn run(args: &[String]) -> Result<(), String> {
    let num = |s: &String| s.parse::<i128>().map_err(|_| format!("not a number: {s}"));
    let (lx, ly, r, cut) = (num(&args[1])? as usize, num(&args[2])? as usize, num(&args[3])?, num(&args[4])? as usize);
    let (threads, sweeps) = (num(&args[5])? as usize, num(&args[6])? as usize);
    let auto = args[7] == "auto";
    let targets: Vec<i128> = if auto {
        // a first guess of lambda_1(H'), 3nl - 3nl/(4r) + 3r - 21/r; it only places the targets
        let nl = 2 * (lx * ly) as i128;
        let guess = 3 * nl * U - 3 * nl * U / (4 * r) + 3 * r * U - 21 * U / r;
        [20_000_000i128, 16_000_000, 12_000_000, 9_000_000, 7_000_000, 5_500_000, 4_500_000, 3_500_000, 2_800_000, 2_200_000, 1_700_000, 1_300_000, 1_000_000, 750_000, 500_000, 300_000, 150_000]
            .iter()
            .map(|off| guess - off)
            .collect()
    } else {
        args[7].split(',').map(|x| x.parse::<i128>().map_err(|_| "bad target".to_string())).collect::<Result<_, _>>()?
    };
    if lx < 3 || ly < 3 || cut < 4 {
        return Err("each side at least 3 and CUT at least 4".into());
    }
    let g = Group::new();
    let t = Torus::new(lx, ly);
    let md = build(&t, &g, cut, threads)?;
    // every lump inside the cut must hold exactly its count of configurations
    let mut found: HashMap<(usize, usize), u128> = HashMap::new();
    for st in &md.states {
        *found.entry((st.a as usize, st.b as usize)).or_insert(0) += u128::from(st.weight);
    }
    let (mut bad, mut total) = (0u32, 0u128);
    for a in (0..=t.p).step_by(2) {
        for b in 0..=(t.p - a) {
            if cost(a, b) == 0 || cost(a, b) > cut {
                continue;
            }
            let (want, got) = (expected(t.p, a, b), found.get(&(a, b)).copied().unwrap_or(0));
            total += got;
            if want != got {
                bad += 1;
                println!("BIG_LUMP_BAD|a={a}|b={b}|configs_found={got}|configs_expected={want}|json=0");
            }
        }
    }
    let moves: usize = md.states.iter().map(|s| s.exact.len() + s.above.iter().map(|x| x.1 as usize).sum::<usize>() + s.to_flat as usize).sum();
    println!(
        "BIG_SPACE|Lx={lx}|Ly={ly}|plaquettes={}|links={}|symmetries={}|cut={cut}|classes={}|configs_in_the_cut={total}|lumps_with_a_wrong_count={bad}|moves={moves}|moves_per_class_is_3nl={}|lumps_above={}|vertices_above={}|json=0",
        t.p,
        t.nl,
        t.auts.len(),
        md.states.len(),
        i32::from(moves == md.states.len() * 3 * t.nl),
        md.lumps.len(),
        md.lumps.iter().map(|l| l.2.len()).sum::<usize>()
    );
    if bad != 0 {
        return Err("the exact region is incomplete at this cut".into());
    }
    let mut best: Option<Fit> = None;
    let mut best_target = targets[0];
    let mut tried: Vec<i128> = Vec::new();
    let mut queue = targets.clone();
    let mut refined = !auto;
    loop {
        for &target in &queue {
            tried.push(target);
            match fit(&md, &t, r, target, sweeps, threads) {
                Some(f) => {
                    println!("BIG_FIT|target_micro={target}|cert_micro={}|min_exact_rows={}|min_vertex_rows={}|json=0", f.cert, f.min_exact, f.min_lump);
                    if best.as_ref().map_or(true, |b| f.cert > b.cert) {
                        best = Some(f);
                        best_target = target;
                    }
                }
                None => println!("BIG_FIT|target_micro={target}|cert_micro=NONE|reason=target_above_a_row_or_values_out_of_range|json=0"),
            }
        }
        if refined {
            break;
        }
        // second pass: eight targets between the best target and the next one above it
        refined = true;
        let above = tried.iter().copied().filter(|x| *x > best_target).min().unwrap_or(best_target + 150_000);
        queue = (1..9).map(|k| best_target + (above - best_target) * k / 9).collect();
    }
    let Some(f) = best else { return Err("no target gave a vector".into()) };
    // the certificate, recomputed from the kept vector
    let (min_exact, min_lump, electric) = certificate(&md, &t, r, &f.psi, &f.g);
    let cert = min_exact.min(min_lump);
    let positive = f.psi.iter().all(|x| *x > 0) && f.g.iter().all(|x| *x > 0);
    let upper = rayleigh(&md, &t, r, &f.psi);
    let nl = t.nl as i128;
    let e8_flat = 3 * nl * U;
    let e8 = if lx.min(ly) >= 4 { (3 * nl + 2 * r) * U - isqrt((4 * r * r + 3 * nl) * U * U) } else { e8_flat };
    println!(
        "BIG_CERT|Lx={lx}|Ly={ly}|r={r}|cut={cut}|vector_positive={}|E9_lower_micro={cert}|min_exact_rows={min_exact}|min_vertex_rows={min_lump}|lambda1_upper_micro={upper}|enclosure_width_micro={}|E8_upper_micro={e8}|gap_lower_micro={}|rows_above_the_cut=polytope_vertices|json=0",
        i32::from(positive),
        upper - cert,
        cert - e8
    );
    if !positive || upper < cert {
        return Err("the vector is not positive, or the upper end is below the certificate".into());
    }
    // the same vector at other couplings: every row is linear in r
    if let Some(list) = args.get(8) {
        for r2 in list.split(',').filter_map(|x| x.parse::<i128>().ok()) {
            let c2 = electric.iter().enumerate().filter(|(_, e)| **e != i128::MAX).map(|(m, e)| e + r2 * m as i128 * U).min().unwrap_or(0);
            let e8b = if lx.min(ly) >= 4 { (3 * nl + 2 * r2) * U - isqrt((4 * r2 * r2 + 3 * nl) * U * U) } else { e8_flat };
            println!("BIG_REPLAY|vector_fitted_at_r={r}|read_at_r={r2}|E9_lower_micro={c2}|E8_upper_micro={e8b}|gap_lower_micro={}|json=0", c2 - e8b);
        }
    }
    let _ = &md.states[0].key;
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 8 {
        println!("BIG_REFUSE|reason=usage_large_tori_LX_LY_R_CUT_THREADS_SWEEPS_TARGETS_[REPLAY_R]|json=0");
        std::process::exit(2);
    }
    if let Err(why) = run(&args) {
        println!("BIG_REFUSE|reason={}|json=0", why.replace(' ', "_"));
        std::process::exit(1);
    }
}


#[cfg(test)]
mod g5_fault_controls {
    use super::*;

    // The same objective/census rules used by production, but invoked
    // on deliberately modified records to prove they turn red.
    fn audit_exact_region(t: &Torus, group: &Group, model: &Model) -> Result<(), &'static str> {
        let mut counts: HashMap<(usize, usize), u128> = HashMap::new();
        let mut decoded = vec![0u8; t.nl];
        for state in &model.states {
            t.decode(&state.key, &mut decoded);
            let (a, b) = lump_of(t, group, &decoded);
            if a != usize::from(state.a) || b != usize::from(state.b) {
                return Err("WRONG_MAGNETIC_LUMP");
            }
            let raw_moves = state.exact.len()
                + state.above.iter().map(|(_, n)| usize::from(*n)).sum::<usize>()
                + usize::from(state.to_flat);
            if raw_moves != 3 * t.nl {
                return Err("INCOMPLETE_3NL_ELECTRIC_ROW");
            }
            *counts.entry((a, b)).or_insert(0) += u128::from(state.weight);
        }
        for a in (0..=t.p).step_by(2) {
            for b in 0..=t.p - a {
                let mass = cost(a, b);
                if mass == 0 || mass > model.cut {
                    continue;
                }
                if counts.get(&(a, b)).copied().unwrap_or(0) != expected(t.p, a, b) {
                    return Err("MISSING_ORBIT_REPRESENTATIVE");
                }
            }
        }
        Ok(())
    }

    fn check_reported_lower(
        claimed: i128,
        exact_source_minimum: i128,
        above_cut_minimum: i128,
    ) -> Result<(), &'static str> {
        // Claiming even ONE micro-unit above the exact integer minima
        // must be rejected. A smaller floor is conservative and valid.
        if claimed > exact_source_minimum.min(above_cut_minimum) {
            return Err("UNJUSTIFIED_PLUS_ONE_LOWER_CLAIM");
        }
        Ok(())
    }

    #[test]
    fn g5_exact_small_torus_rejects_sabotage() {
        let torus = Torus::new(4, 4);
        let group = Group::new();
        let mut model = build(&torus, &group, 4, 1)
            .expect("4x4 cut4 must construct a complete original orbit region");
        assert!(!model.states.is_empty());
        audit_exact_region(&torus, &group, &model)
            .expect("original unmodified small-torus census and raw row count");

        // Negative control 1: a reported lower floor inflated by exactly
        // one micro-unit. It MUST refuse, without rounding slack.
        let actual_vector = vec![1i128 << 40; model.states.len()];
        let actual_above_weights = vec![1i128 << 35; model.lumps.len()];
        let (source, above, _) = certificate(
            &model,
            &torus,
            9,
            &actual_vector,
            &actual_above_weights,
        );
        let measured_floor = source.min(above);
        assert!(measured_floor != i128::MIN && measured_floor != i128::MAX);
        assert!(check_reported_lower(measured_floor, source, above).is_ok());
        assert_eq!(
            check_reported_lower(measured_floor + 1, source, above),
            Err("UNJUSTIFIED_PLUS_ONE_LOWER_CLAIM")
        );

        // Negative control 2: remove one genuine positive-weight orbit
        // representative. The independent per-lump combinatorial count
        // MUST fail even though every retained row is unchanged.
        let removed = model.states.pop().expect("at least one actual class");
        assert_eq!(
            audit_exact_region(&torus, &group, &model),
            Err("MISSING_ORBIT_REPRESENTATIVE")
        );
        model.states.push(removed);
        assert!(audit_exact_region(&torus, &group, &model).is_ok());

        // Negative control 3: corrupt the stored magnetic (a,b) count,
        // leaving the canonical source key and electric moves untouched.
        let correct_a = model.states[0].a;
        model.states[0].a = correct_a + 2;
        assert_eq!(
            audit_exact_region(&torus, &group, &model),
            Err("WRONG_MAGNETIC_LUMP")
        );
        model.states[0].a = correct_a;

        // Negative control 4: remove one raw electric move from a
        // source row without changing its other destinations.
        assert!(audit_exact_region(&torus, &group, &model).is_ok());
        let original_flat = model.states[0].to_flat;
        let mut popped = None;
        if original_flat > 0 {
            model.states[0].to_flat -= 1;
        } else if !model.states[0].exact.is_empty() {
            popped = model.states[0].exact.pop();
        } else {
            let (dest, mult) = model.states[0]
                .above
                .first_mut()
                .expect("every nonflat orbit must have electric moves");
            let _ = dest;
            assert!(*mult > 0);
            *mult -= 1;
        }
        assert_eq!(
            audit_exact_region(&torus, &group, &model),
            Err("INCOMPLETE_3NL_ELECTRIC_ROW")
        );
        // Restore only where the source row was modified.
        if original_flat > 0 {
            model.states[0].to_flat = original_flat;
        } else if let Some(index) = popped {
            model.states[0].exact.push(index);
        } else {
            model.states[0].above[0].1 += 1;
        }
        assert!(audit_exact_region(&torus, &group, &model).is_ok());
    }
}
