//! Standalone independent Barta verifier for HYBRID_BIG vector dumps.
//! No inclusion or import of src/main.rs, no use of its canonicalizer.
//! Uses its own S3 permutation multiplication, torus geometry, symmetry
//! action, spanning-tree gauge, orbit size and exact integer polytope.
//! Rust 1.81, std only. Deliberately fail closed on missing inputs.
//!
//! cargo run --release --bin vectorcheck -- VECTOR.dump [--shard=N/T]
//!   [--fault=plus1|drop|magnetic|move]
//! A shard checks only N modulo T of the low rows and ALL high rows;
//! only the fan-in of every shard is a complete source certificate.

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;

const SCALE: i128 = 1_000_000;
const TOP: i128 = 1i128 << 86;
const PERM: [[u8; 3]; 6] = [
    [0, 1, 2], [0, 2, 1], [1, 0, 2],
    [1, 2, 0], [2, 0, 1], [2, 1, 0],
];

type Check<T> = Result<T, String>;

fn demand(ok: bool, why: &str) -> Check<()> {
    if ok { Ok(()) } else { Err(why.to_owned()) }
}
fn vadd(a: i128, b: i128) -> Check<i128> {
    a.checked_add(b).ok_or_else(|| "I128_ADDITION_OVERFLOW".to_owned())
}
fn vmul(a: i128, b: i128) -> Check<i128> {
    a.checked_mul(b).ok_or_else(|| "I128_MULTIPLICATION_OVERFLOW".to_owned())
}

struct S3 {
    prod: [[u8; 6]; 6],
    inv: [u8; 6],
    label: [u8; 6],
    tr: [u8; 3],
}
impl S3 {
    fn new() -> Check<Self> {
        let mut prod = [[0u8; 6]; 6];
        let mut inv = [0u8; 6];
        let mut label = [0u8; 6];
        let mut trans = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                let p = [
                    PERM[i][PERM[j][0] as usize],
                    PERM[i][PERM[j][1] as usize],
                    PERM[i][PERM[j][2] as usize],
                ];
                prod[i][j] = PERM.iter().position(|v| *v == p)
                    .ok_or("BROKEN_PERMUTATION_GROUP")? as u8;
            }
            let fixes = (0..3).filter(|&k| PERM[i][k] == k as u8).count();
            label[i] = match fixes { 3 => 0, 1 => 1, 0 => 2, _ => {
                return Err("UNEXPECTED_S3_CYCLE_TYPE".to_owned());
            }};
            if label[i] == 1 { trans.push(i as u8); }
        }
        for i in 0..6 {
            inv[i] = (0..6).find(|&j| prod[i][j] == 0 && prod[j][i] == 0)
                .ok_or("S3_INVERSE_MISSING")? as u8;
        }
        demand(trans.len() == 3, "S3_TRANS_COUNT")?;
        Ok(Self { prod, inv, label, tr: [trans[0], trans[1], trans[2]] })
    }
    fn mul(&self, a: u8, b: u8) -> u8 {
        self.prod[a as usize][b as usize]
    }
    fn inverse(&self, a: u8) -> u8 {
        self.inv[a as usize]
    }
    fn conj(&self, c: u8, x: u8) -> u8 {
        self.mul(self.mul(c, x), self.inverse(c))
    }
}

#[derive(Clone)]
struct Symmetry {
    edge: Vec<(usize, bool)>,
    face: Vec<usize>,
}
struct SquareTorus {
    lx: usize,
    ly: usize,
    p: usize,
    edges: usize,
    free: Vec<usize>,
    tree: Vec<(usize, usize, usize)>,
    endpoints: Vec<(usize, usize)>,
    syms: Vec<Symmetry>,
}
impl SquareTorus {
    fn vertex(&self, x: usize, y: usize) -> usize {
        (x % self.lx) * self.ly + y % self.ly
    }
    fn horizontal(&self, x: usize, y: usize) -> usize {
        2 * self.vertex(x, y)
    }
    fn vertical(&self, x: usize, y: usize) -> usize {
        2 * self.vertex(x, y) + 1
    }
    fn corner_set(&self, x: usize, y: usize) -> [usize; 4] {
        let mut corners = [
            self.vertex(x, y),
            self.vertex(x + 1, y),
            self.vertex(x, y + 1),
            self.vertex(x + 1, y + 1),
        ];
        corners.sort_unstable();
        corners
    }
    fn new(lx: usize, ly: usize) -> Check<Self> {
        demand(lx >= 4 && ly >= 4 && lx <= 16 && ly <= 16, "SIDE_UNSUPPORTED")?;
        let mut t = Self {
            lx, ly, p: lx * ly, edges: 2 * lx * ly,
            free: Vec::new(), tree: Vec::new(),
            endpoints: Vec::new(), syms: Vec::new(),
        };
        let mut tree_mask = vec![false; t.edges];
        // This is a specification of the dump's comb gauge, not a
        // call into the producer's tree builder.
        for y in 0..ly - 1 {
            let e = t.vertical(0, y);
            t.tree.push((t.vertex(0, y), t.vertex(0, y + 1), e));
            tree_mask[e] = true;
        }
        for x in 0..lx - 1 {
            for y in 0..ly {
                let e = t.horizontal(x, y);
                t.tree.push((t.vertex(x, y), t.vertex(x + 1, y), e));
                tree_mask[e] = true;
            }
        }
        t.free = (0..t.edges).filter(|&e| !tree_mask[e]).collect();
        demand(t.free.len() == t.p + 1, "COMB_GAUGE_NOT_P_PLUS_1")?;
        t.endpoints = vec![(0, 0); t.edges];
        for x in 0..lx {
            for y in 0..ly {
                t.endpoints[t.horizontal(x,y)] =
                    (t.vertex(x,y),t.vertex(x+1,y));
                t.endpoints[t.vertical(x,y)] =
                    (t.vertex(x,y),t.vertex(x,y+1));
            }
        }
        let mut directed = HashMap::new();
        for (edge, &(a,b)) in t.endpoints.iter().enumerate() {
            demand(directed.insert((a,b),(edge,false)).is_none(),
                   "EDGE_COLLISION")?;
            demand(directed.insert((b,a),(edge,true)).is_none(),
                   "REVERSE_EDGE_COLLISION")?;
        }
        let mut face_lookup = HashMap::new();
        for x in 0..lx {
            for y in 0..ly {
                demand(face_lookup.insert(t.corner_set(x,y),t.vertex(x,y)).is_none(),
                       "FACE_CORNER_COLLISION")?;
            }
        }
        for swap in [false,true] {
            if swap && lx != ly { continue; }
            for xreflect in [false,true] {
                for yreflect in [false,true] {
                    for dx in 0..lx {
                        for dy in 0..ly {
                            let f = |v: usize| -> usize {
                                let x = v / ly;
                                let y = v % ly;
                                let (a,b) = if swap {(y,x)} else {(x,y)};
                                let ax = if xreflect {(lx-a)%lx} else {a};
                                let by = if yreflect {(ly-b)%ly} else {b};
                                t.vertex(ax+dx,by+dy)
                            };
                            let mut edge_map = Vec::with_capacity(t.edges);
                            for &(a,b) in &t.endpoints {
                                edge_map.push(*directed.get(&(f(a),f(b)))
                                    .ok_or("IMAGE_OF_EDGE_NOT_PRESENT")?);
                            }
                            let mut face_map = Vec::with_capacity(t.p);
                            for x in 0..lx {
                                for y in 0..ly {
                                    let corners = t.corner_set(x,y);
                                    let mut changed = [
                                        f(corners[0]),f(corners[1]),
                                        f(corners[2]),f(corners[3]),
                                    ];
                                    changed.sort_unstable();
                                    face_map.push(*face_lookup.get(&changed)
                                        .ok_or("IMAGE_OF_FACE_NOT_PRESENT")?);
                                }
                            }
                            // face order is x*ly+y, matching vertex().
                            t.syms.push(Symmetry{edge:edge_map,face:face_map});
                        }
                    }
                }
            }
        }
        demand(t.syms.len()==(if lx==ly {8} else {4})*t.p,
               "SPATIAL_SYMMETRY_COUNT")?;
        Ok(t)
    }
    fn plaquette(&self,g: &S3,u: &[u8],x: usize,y: usize) -> u8 {
        let a = g.mul(u[self.horizontal(x,y)], u[self.vertical(x+1,y)]);
        let b = g.mul(a, g.inverse(u[self.horizontal(x,y+1)]));
        g.mul(b, g.inverse(u[self.vertical(x,y)]))
    }
    fn classes(&self,g: &S3,u: &[u8]) -> Vec<u8> {
        (0..self.p).map(|face|g.label[
            self.plaquette(g,u,face/self.ly,face%self.ly) as usize
        ]).collect()
    }
    fn from_digits(&self,digits: &[u8]) -> Check<Vec<u8>> {
        demand(digits.len()==self.free.len(),"WRONG_FREE_LINK_DIGIT_COUNT")?;
        let mut u=vec![0;self.edges];
        for (&edge,&digit) in self.free.iter().zip(digits) {
            demand(digit<6,"NON_S3_DIGIT")?;
            u[edge]=digit;
        }
        Ok(u)
    }
    fn fix_gauge(&self,g: &S3,u: &[u8]) -> Vec<u8> {
        let mut vertices=vec![0u8;self.p];
        for &(parent,child,edge) in &self.tree {
            vertices[child]=g.mul(vertices[parent],u[edge]);
        }
        self.free.iter().map(|&edge|{
            let (a,b)=self.endpoints[edge];
            g.mul(g.mul(vertices[a],u[edge]),g.inverse(vertices[b]))
        }).collect()
    }
    fn canonical_and_weight(&self,g: &S3,u: &[u8]) -> Check<(Vec<u8>,u64)> {
        let cls = self.classes(g,u);
        let mut least_pattern: Option<Vec<(usize,u8)>> = None;
        let mut matched = Vec::new();
        for (si,sym) in self.syms.iter().enumerate() {
            let mut transformed:Vec<(usize,u8)> = cls.iter().enumerate()
                .filter(|&(_,kind)| *kind!=0)
                .map(|(idx,&kind)|(sym.face[idx],kind)).collect();
            transformed.sort_unstable();
            match least_pattern.as_ref() {
                None => {
                    least_pattern=Some(transformed);
                    matched.clear();
                    matched.push(si);
                }
                Some(best) if transformed < *best => {
                    least_pattern=Some(transformed);
                    matched.clear();
                    matched.push(si);
                }
                Some(best) if transformed == *best => matched.push(si),
                _ => {}
            }
        }
        demand(!matched.is_empty() && self.syms.len()%matched.len()==0,
               "INCOMPLETE_CANONICAL_PATTERN_ORBIT")?;
        let mut unique = HashSet::new();
        let mut least:Option<Vec<u8>>=None;
        for &si in &matched {
            let mut image=vec![0u8;self.edges];
            for (edge,&(other,reversed)) in self.syms[si].edge.iter().enumerate() {
                image[other] = if reversed {g.inverse(u[edge])} else {u[edge]};
            }
            let fixed=self.fix_gauge(g,&image);
            for c in 0..6u8 {
                let key:Vec<u8> = fixed.iter().map(|&x|g.conj(c,x)).collect();
                if least.as_ref().map_or(true,|v|key < *v) {
                    least=Some(key.clone());
                }
                unique.insert(key);
            }
        }
        let weight=(self.syms.len()/matched.len()) as u64
            *unique.len() as u64;
        demand(weight>0,"ZERO_ORBIT_SIZE")?;
        Ok((least.ok_or("CANON_KEY_MISSING")?,weight))
    }
}

struct ClassRow {
    key: Vec<u8>,
    a:usize,b:usize, weight:u64, value:i128,
}
struct Input {
    lx:usize,ly:usize,r:i128,cut:usize,claim:i128,
    classes:Vec<ClassRow>,
    high:HashMap<(usize,usize),i128>,
}
fn header_value<'a>(line:&'a str,name:&str)->Check<&'a str> {
    line.split('|').find_map(|x|x.strip_prefix(&format!("{name}=")))
        .ok_or_else(||format!("MISSING_HEADER_{name}"))
}
fn parse_num<T:std::str::FromStr>(s:&str,tag:&str)->Check<T> {
    s.parse().map_err(|_|format!("BAD_INTEGER_{tag}"))
}
fn load(path:&str)->Check<Input> {
    let body=fs::read_to_string(path).map_err(|e|e.to_string())?;
    let mut lines=body.lines();
    let hdr=lines.next().ok_or("EMPTY_VECTOR_DUMP")?;
    demand(hdr.starts_with("HYBRID_BIG|"),"MISSING_HYBRID_BIG_HEADER")?;
    demand(header_value(hdr,"key")?=="free_link_digits_of_the_comb_gauge",
           "UNRECOGNIZED_KEY_GAUGE")?;
    demand(header_value(hdr,"json")?=="0","DUMP_HEADER_JSON_MARKER")?;
    let lx=parse_num(header_value(hdr,"Lx")?,"LX")?;
    let ly=parse_num(header_value(hdr,"Ly")?,"LY")?;
    let r=parse_num(header_value(hdr,"r")?,"R")?;
    let cut=parse_num(header_value(hdr,"cut")?,"CUT")?;
    let claim=parse_num(header_value(hdr,"claim_micro")?,"CLAIM")?;
    let want_p:usize=parse_num(header_value(hdr,"classes")?,"CLASSES")?;
    let want_g:usize=parse_num(header_value(hdr,"lumps_above")?,"LUMPS")?;
    demand(want_p>0 && want_g>0 && want_p<30_000_000 && want_g<100_000,
           "IMPLAUSIBLE_DUMP_CENSUS")?;
    let mut classes=Vec::with_capacity(want_p);
    let mut high=HashMap::new();
    let mut started_g=false;
    for line in lines {
        let parts:Vec<_>=line.split('|').collect();
        match parts.first().copied() {
            Some("P") => {
                demand(!started_g && parts.len()==6,"BAD_P_RECORD")?;
                demand(parts[1].bytes().all(|b|(b'0'..=b'5').contains(&b)),
                       "BAD_P_LINK_DIGIT")?;
                let key=parts[1].bytes().map(|b|b-b'0').collect();
                let a=parse_num(parts[2],"P_A")?;
                let b=parse_num(parts[3],"P_B")?;
                let weight=parse_num(parts[4],"P_WEIGHT")?;
                let value=parse_num(parts[5],"P_VALUE")?;
                classes.push(ClassRow{key,a,b,weight,value});
                demand(classes.len()<=want_p,"EXTRA_P_RECORD")?;
            }
            Some("G") => {
                started_g=true;
                demand(parts.len()==4,"BAD_G_RECORD")?;
                let a=parse_num(parts[1],"G_A")?;
                let b=parse_num(parts[2],"G_B")?;
                let value=parse_num(parts[3],"G_VALUE")?;
                demand(high.insert((a,b),value).is_none(),
                       "DUPLICATED_ABOVE_LUMP")?;
                demand(high.len()<=want_g,"EXTRA_G_RECORD")?;
            }
            _ => return Err("UNKNOWN_VECTOR_DUMP_RECORD".to_owned()),
        }
    }
    demand(classes.len()==want_p && high.len()==want_g,
           "INCOMPLETE_VECTOR_DUMP")?;
    demand(r>=0 && cut>=4 && claim>0,"INVALID_CLAIM_PARAMETERS")?;
    Ok(Input{lx,ly,r,cut,claim,classes,high})
}
fn nchoose(n:usize,k:usize)->Check<u128> {
    let k=k.min(n-k);
    let mut value=1u128;
    for i in 1..=k {
        value=value.checked_mul((n+1-i) as u128)
            .ok_or("BINOMIAL_OVERFLOW")?/(i as u128);
    }
    Ok(value)
}
fn expected_lump(p:usize,a:usize,b:usize)->Check<u128> {
    let comb=nchoose(p,a)?.checked_mul(nchoose(p-a,b)?)
        .ok_or("EXPECTED_LUMP_BINOMIAL_OVERFLOW")?;
    let base=if a==0 {
        let two=2u128.checked_pow(b as u32).ok_or("LUMP_POWER_OVERFLOW")?;
        let x=12u128.checked_mul(two).ok_or("LUMP_POWER_OVERFLOW")?;
        if b%2==0 {x+6}else{x-6}
    } else {
        let three=3u128.checked_pow(a as u32).ok_or("LUMP_POWER_OVERFLOW")?;
        let two=2u128.checked_pow(b as u32).ok_or("LUMP_POWER_OVERFLOW")?;
        12u128.checked_mul(three).and_then(|v|v.checked_mul(two))
            .ok_or("LUMP_POWER_OVERFLOW")?
    };
    comb.checked_mul(base).ok_or("EXPECTED_LUMP_OVERFLOW".to_owned())
}
fn mcount(a:usize,b:usize)->usize{2*a+3*b}
fn ceil_2sqrt(n:usize)->usize {
    let mut low=0usize;
    let mut high=2*n+1;
    while low+1<high {
        let mid=(low+high)/2;
        if mid*mid>=4*n { high=mid; } else {low=mid;}
    }
    high
}
fn adjacency(s:usize,t:&SquareTorus)->i128 {
    if s==0 {return 0;}
    if s==t.p{return (2*s) as i128;}
    (2*s-ceil_2sqrt(s).min(t.lx.min(t.ly))
        .min(ceil_2sqrt(t.p-s))) as i128
}
fn cross(a:[i128;3],b:[i128;3])->[i128;3]{
    [a[1]*b[2]-a[2]*b[1],
     a[2]*b[0]-a[0]*b[2],
     a[0]*b[1]-a[1]*b[0]]
}
fn dot(a:[i128;3],b:[i128;3])->i128{
    a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
}
fn extreme_vertices(t:&SquareTorus,a:usize,b:usize)->Check<Vec<[i128;4]>> {
    let na=4*a as i128;
    let nb=4*b as i128;
    let nl=t.edges as i128;
    let row:[[i128;4];9]=[
        [-1,0,0,0],[0,-1,0,0],[0,0,-1,0],
        [1,0,0,adjacency(a,t)],[0,0,1,adjacency(b,t)],
        [1,1,1,adjacency(a+b,t)],[2,1,0,na],
        [0,1,2,nb],[-1,-1,-1,nl-na-nb],
    ];
    let mut unique=HashSet::new();
    for i in 0..9 {
        for j in i+1..9 {
            for k in j+1..9 {
                let n1=[row[i][0],row[i][1],row[i][2]];
                let n2=[row[j][0],row[j][1],row[j][2]];
                let n3=[row[k][0],row[k][1],row[k][2]];
                let b1=row[i][3];let b2=row[j][3];let b3=row[k][3];
                let c23=cross(n2,n3);
                let c31=cross(n3,n1);
                let c12=cross(n1,n2);
                let d=dot(n1,c23);
                if d==0 {continue;}
                let mut xyz=[b1*c23[0]+b2*c31[0]+b3*c12[0],
                             b1*c23[1]+b2*c31[1]+b3*c12[1],
                             b1*c23[2]+b2*c31[2]+b3*c12[2]];
                let den=d.abs();
                if d<0 { for v in &mut xyz{*v = -*v;} }
                let valid=row.iter().all(|c|{
                    c[0]*xyz[0]+c[1]*xyz[1]+c[2]*xyz[2]<=c[3]*den
                });
                if valid{unique.insert([den,xyz[0],xyz[1],xyz[2]]);}
            }
        }
    }
    Ok(unique.into_iter().collect())
}
fn high_values(t:&SquareTorus,input:&Input, maxima:&HashMap<(usize,usize),i128>)
    -> Check<(i128,usize,usize)> {
    let group=S3::new()?;
    let mut expected_high=0;
    let mut min=i128::MAX;
    let mut vertex_count=0usize;
    let mut high_rows=0usize;
    for a in (0..=t.p).step_by(2) {
        for b in 0..=t.p-a {
            if mcount(a,b)<=input.cut {continue;}
            let vertices=extreme_vertices(t,a,b)?;
            if vertices.is_empty() {
                demand(!input.high.contains_key(&(a,b)),
                       "G_LUMP_WITHOUT_FEASIBLE_VERTEX")?;
                continue;
            }
            expected_high+=1;
            let value=*input.high.get(&(a,b)).ok_or_else(||
                format!("MISSING_G_LUMP_{a}_{b}"))?;
            demand(value>0 && value<=TOP,"G_VALUE_OUTSIDE_SAFE_TOP")?;
            let tr=group.tr[0];
            let ct=(0..6u8).find(|&q|group.label[q as usize]==2)
                .ok_or("MISSING_THREE_CYCLE")?;
            let tt_other=group.tr[1];
            let patterns=[
                (0u8,0u8),(0u8,tr),(0u8,ct),
                (tr,tr),(tr,tt_other),(tr,ct),(ct,group.inverse(ct))
            ];
            let target=|aa:i128,bb:i128|->Check<i128> {
                if aa<0||bb<0||aa+bb>t.p as i128||aa%2!=0 {
                    return Ok(TOP);
                }
                let dest=(aa as usize,bb as usize);
                if mcount(dest.0,dest.1)==0 {return Ok(0);}
                if mcount(dest.0,dest.1)<=input.cut {
                    Ok(*maxima.get(&dest).ok_or(
                        format!("MISSING_LOW_LUMP_MAX_{}_{}",dest.0,dest.1))?)
                } else {
                    Ok(*input.high.get(&dest).unwrap_or(&TOP))
                }
            };
            let mut weights=[0i128;7];
            for (i,&(x,y)) in patterns.iter().enumerate() {
                let olda=(group.label[x as usize]==1) as i128
                         +(group.label[y as usize]==1) as i128;
                let oldb=(group.label[x as usize]==2) as i128
                         +(group.label[y as usize]==2) as i128;
                let mut sum=0i128;
                for &tr in &group.tr {
                    let xx=group.mul(tr,x);
                    let yy=group.mul(tr,y);
                    let da=(group.label[xx as usize]==1) as i128
                           +(group.label[yy as usize]==1) as i128-olda;
                    let db=(group.label[xx as usize]==2) as i128
                           +(group.label[yy as usize]==2) as i128-oldb;
                    sum=vadd(sum,target(a as i128+da,b as i128+db)?)?;
                }
                weights[i]=sum;
            }
            let sums=[weights[0],weights[1],weights[2],
                      weights[3].max(weights[4]),weights[5],weights[6]];
            let na=4*a as i128;
            let nb=4*b as i128;
            let nl=t.edges as i128;
            for [den,x,y,z] in vertices {
                let counts=[den*(nl-na-nb)+x+y+z,
                    den*na-2*x-y,den*nb-y-2*z,x,y,z];
                demand(counts.iter().all(|&c|c>=0),"NEGATIVE_PATTERN_MULTIPLICITY")?;
                let mut sum=0i128;
                for (count,weight) in counts.into_iter().zip(sums) {
                    sum=vadd(sum,vmul(count,weight)?)?;
                }
                let denominator=vmul(den,value)?;
                demand(denominator>0,"ZERO_VERTEX_DENOMINATOR")?;
                let electric_floor=(3*nl)*SCALE - div_ceil(vmul(sum,SCALE)?,denominator)?;
                let ratio=vadd(electric_floor,input.r*mcount(a,b) as i128*SCALE)?;
                min=min.min(ratio);
                vertex_count+=1;
                high_rows+=1;
            }
        }
    }
    demand(expected_high==input.high.len(),"EXTRA_G_LUMP")?;
    demand(min!=i128::MAX,"NO_ABOVE_CUT_ROWS")?;
    Ok((min,expected_high,vertex_count.max(high_rows)))
}
fn div_ceil(a:i128,b:i128)->Check<i128> {
    demand(b>0&&a>=0,"INVALID_POSITIVE_CEIL")?;
    Ok(a/b+i128::from(a%b!=0))
}
fn check_input(input:&Input,t:&SquareTorus)->Check<HashMap<Vec<u8>,usize>> {
    let mut map=HashMap::with_capacity(input.classes.len());
    for (i,row) in input.classes.iter().enumerate() {
        demand(row.key.len()==t.free.len(),"BAD_KEY_WIDTH")?;
        demand(row.a<=t.p&&row.b<=t.p&&row.a+row.b<=t.p&&row.a%2==0,
               "INVALID_P_LUMP")?;
        demand(mcount(row.a,row.b)>0&&mcount(row.a,row.b)<=input.cut,
               "P_OUTSIDE_EXACT_REGION")?;
        demand(row.value>0 && row.value<=TOP,"P_VALUE_OUTSIDE_SAFE_TOP")?;
        demand(row.weight>0,"ZERO_P_WEIGHT")?;
        demand(map.insert(row.key.clone(),i).is_none(),"DUPLICATE_P_KEY")?;
    }
    Ok(map)
}
fn check_count(input:&Input,t:&SquareTorus)->Check<()> {
    let mut counts=HashMap::<(usize,usize),u128>::new();
    for row in &input.classes {
        let v=counts.entry((row.a,row.b)).or_default();
        *v=v.checked_add(row.weight as u128)
            .ok_or("ORBIT_CENSUS_OVERFLOW")?;
    }
    for a in (0..=t.p).step_by(2) {
        for b in 0..=t.p-a {
            let m=mcount(a,b);
            if m==0||m>input.cut{continue;}
            let actual=counts.get(&(a,b)).copied().unwrap_or(0);
            let wanted=expected_lump(t.p,a,b)?;
            demand(actual==wanted,&format!(
                "LOW_LUMP_COUNT_FAIL_a{a}_b{b}_got{actual}_want{wanted}"))?;
        }
    }
    Ok(())
}
fn main_check(path:&str,shard:(usize,usize),fault:Option<&str>)->Check<()> {
    let mut input=load(path)?;
    let torus=SquareTorus::new(input.lx,input.ly)?;
    let group=S3::new()?;
    if fault==Some("drop") {
        input.classes.remove(0);
    }
    if fault==Some("magnetic") {
        input.classes[0].a+=2;
    }
    let lookup=check_input(&input,&torus)?;
    check_count(&input,&torus)?;
    let mut max_low=HashMap::<(usize,usize),i128>::new();
    for row in &input.classes {
        let e=max_low.entry((row.a,row.b)).or_insert(0);
        *e=(*e).max(row.value);
    }
    let (high_min,gcount,high_vertices)=high_values(&torus,&input,&max_low)?;
    let (mut low_min,mut checked,mut moves)=(i128::MAX,0usize,0usize);
    for (idx,row) in input.classes.iter().enumerate() {
        if idx%shard.1!=shard.0 {continue;}
        let mut u=torus.from_digits(&row.key)?;
        let cls=torus.classes(&group,&u);
        let a=cls.iter().filter(|&&k|k==1).count();
        let b=cls.iter().filter(|&&k|k==2).count();
        demand(a==row.a&&b==row.b,"DECODED_MAGNETIC_LUMP_MISMATCH")?;
        let (own,weight)=torus.canonical_and_weight(&group,&u)?;
        demand(own==row.key,"NONCANONICAL_P_KEY")?;
        demand(weight==row.weight,"CLASS_ORBIT_WEIGHT_MISMATCH")?;
        let mut destination_sum=0i128;
        let mut row_moves=0usize;
        for edge in 0..torus.edges {
            let old=u[edge];
            for &tr in &group.tr {
                u[edge]=group.mul(tr,old);
                let changed=torus.classes(&group,&u);
                let na=changed.iter().filter(|&&k|k==1).count();
                let nb=changed.iter().filter(|&&k|k==2).count();
                let m=mcount(na,nb);
                let v=if m==0 {
                    0
                }else if m<=input.cut {
                    let (key,_)=torus.canonical_and_weight(&group,&u)?;
                    let j=*lookup.get(&key).ok_or("ELECTRIC_TARGET_P_MISSING")?;
                    let target=&input.classes[j];
                    demand((na,nb)==(target.a,target.b),
                           "ELECTRIC_TARGET_LUMP_DISAGREES")?;
                    target.value
                }else{
                    *input.high.get(&(na,nb)).ok_or("ELECTRIC_G_LUMP_MISSING")?
                };
                destination_sum=vadd(destination_sum,v)?;
                row_moves+=1;
                u[edge]=old;
            }
        }
        if fault==Some("move") && checked==0 {row_moves-=1;}
        demand(row_moves==3*torus.edges,"INCOMPLETE_3NL_RAW_MOVES")?;
        let energy=(3*torus.edges as i128
            +input.r*mcount(a,b) as i128)*SCALE;
        let floored=energy-div_ceil(vmul(SCALE,destination_sum)?,row.value)?;
        low_min=low_min.min(floored);
        checked+=1;
        moves+=row_moves;
    }
    demand(checked>0,"EMPTY_SOURCE_SHARD")?;
    let full=shard.1==1;
    let minimum=low_min.min(high_min);
    let effective_claim=if fault==Some("plus1"){input.claim+1}else{input.claim};
    demand(minimum>=effective_claim,"CLAIM_EXCEEDS_INDEPENDENT_ROW_MINIMUM")?;
    if full {
        demand(minimum==input.claim,"DUMP_CLAIM_NOT_THE_INTEGER_ROW_MINIMUM")?;
    }
    println!("VECTORCHECK|Lx={}|Ly={}|r={}|cut={}|shard={}/{}|classes_total={}|classes_checked={checked}|raw_moves_checked={moves}|lumps_above={gcount}|vertices_above={high_vertices}|low_min_micro={low_min}|high_min_micro={high_min}|claim_micro={}|min_checked_micro={minimum}|max_targets=ACTUAL_LOW_ORBIT_MAX|S3=INDEPENDENT|canonical=INDEPENDENT|result=PASS|scope={}|json=0",
        input.lx,input.ly,input.r,input.cut,shard.0,shard.1,input.classes.len(),input.claim,
        if full {"COMPLETE"} else {"SHARD_ONLY"});
    Ok(())
}
fn run()->Check<()> {
    let args:Vec<String>=env::args().collect();
    demand(args.len()>=2,"USAGE_VECTORCHECK_DUMP_SHARD_OPTION")?;
    let mut shard=(0usize,1usize);
    let mut fault=None;
    for arg in &args[2..] {
        if let Some(v)=arg.strip_prefix("--shard=") {
            let (a,b)=v.split_once('/').ok_or("BAD_SHARD_SYNTAX")?;
            shard=(parse_num(a,"SHARD_INDEX")?,parse_num(b,"SHARD_TOTAL")?);
        }else if let Some(v)=arg.strip_prefix("--fault="){
            demand(["plus1","drop","magnetic","move"].contains(&v),
                   "UNKNOWN_FAULT_MODE")?;
            fault=Some(v);
        }else{return Err("UNKNOWN_CHECKER_ARGUMENT".to_owned());}
    }
    demand(shard.1>0&&shard.0<shard.1,"INVALID_SHARD_INDEX")?;
    main_check(&args[1],shard,fault)
}
fn main(){
    if let Err(err)=run(){
        println!("VECTORCHECK_REFUSE|reason={}|json=0",err.replace(' ',"_"));
        std::process::exit(1);
    }
}
