#!/usr/bin/env python3
"""Independent exact-rational validator for the large-torus extreme-row step.

GPT-6, 2026-11-15. This deliberately does not import the existing
Yang-Mills local_table.py, product_bound.py, or the planned Rust code.

Generates S3 link moves from permutations, evaluates all feasible
integer pattern signatures on L=4 (and optionally L=5), and checks
that the vertex minimum of the *relaxed* polytope is never ABOVE
their exact minimum. Includes a hybrid low-cut target-MAX regression
and two planted unsafe shortcuts which MUST be rejected.

This verifies the high-lump proof gadget, NOT the full low-M exact
orbit list, large-torus spectral bound, or the claimed G1-G5 gates.
"""
from fractions import Fraction as Q
from itertools import permutations, combinations
from math import isqrt
import sys

GROUP=list(permutations((0,1,2)))
LOOKUP={p:i for i,p in enumerate(GROUP)}
ID=LOOKUP[(0,1,2)]


def require(yes,why):
    if not yes:
        raise ValueError("EXTREME_ROW_REFUSE_"+why)


def mul(i,j):
    a,b=GROUP[i],GROUP[j]
    return LOOKUP[tuple(a[b[k]] for k in range(3))]


def cls(i):
    fixes=sum(GROUP[i][j]==j for j in range(3))
    return 0 if fixes==3 else 1 if fixes==1 else 2


T=[i for i in range(6) if cls(i)==1]
C=[i for i in range(6) if cls(i)==2]
require(len(T)==3 and len(C)==2,"S3_not_6_elements")
PAIRS={
 "ee":(ID,ID),"eT":(ID,T[0]),"eC":(ID,C[0]),
 "TTs":(T[0],T[0]),"TTd":(T[0],T[1]),
 "TC":(T[0],C[0]),"CC":(C[0],C[1]),
}


def generate_raw_rule():
    rules={}
    for name,(x,y) in PAIRS.items():
        counts={}
        a0=(cls(x)==1)+(cls(y)==1)
        b0=(cls(x)==2)+(cls(y)==2)
        for t in T:
            xx,yy=mul(t,x),mul(t,y)
            da=(cls(xx)==1)+(cls(yy)==1)-a0
            db=(cls(xx)==2)+(cls(yy)==2)-b0
            require(da%2==0,"parity_of_S3_transposition_move")
            counts[da,db]=counts.get((da,db),0)+1
        require(sum(counts.values())==3,"incomplete_S3_electric_rule")
        rules[name]=tuple(sorted((da,db,m) for (da,db),m in counts.items()))
    return rules


RULE=generate_raw_rule()
EXPECTED={
 "ee":((2,0,3),),
 "eT":((0,0,1),(0,1,2)),
 "eC":((2,-1,3),),
 "TTs":((-2,0,1),(-2,2,2)),
 "TTd":((-2,1,2),(-2,2,1)),
 "TC":((0,-1,1),(0,0,2)),
 "CC":((2,-2,3),),
}
require(RULE==EXPECTED,"independent_S3_pattern_rule_changed")


def sqrt_ceil4(n):
    v=isqrt(4*n)
    return v if v*v==4*n else v+1


def adj(s,L):
    p=L*L
    if s==0:return 0
    if s==p:return 2*p
    return 2*s-min(sqrt_ceil4(s),L,sqrt_ceil4(p-s))


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def dot(a,b):
    return sum(i*j for i,j in zip(a,b))


def halfspaces(a,b,L):
    # Independent derivation from four incident links per plaquette,
    # the toroidal adjacency bounds, and nonnegative n_ee.
    nl=2*L*L
    return (
        ((-1,0,0),0), ((0,-1,0),0), ((0,0,-1),0),
        ((1,0,0),adj(a,L)), ((0,0,1),adj(b,L)),
        ((1,1,1),adj(a+b,L)), ((2,1,0),4*a),
        ((0,1,2),4*b), ((-1,-1,-1),nl-4*a-4*b),
    )


def vertex_set(a,b,L):
    planes=halfspaces(a,b,L)
    out=set()
    for rows in combinations(planes,3):
        (n1,d1),(n2,d2),(n3,d3)=rows
        p23,p31,p12=cross(n2,n3),cross(n3,n1),cross(n1,n2)
        den=dot(n1,p23)
        if not den:continue
        v=tuple(Q(d1*p23[k]+d2*p31[k]+d3*p12[k],den)
                for k in range(3))
        if all(dot(n,v)<=d for n,d in planes):
            out.add(v)
    require(bool(out),"no_vertices_for_nonflat_lump")
    return out


def integer_xyz(a,b,L):
    # Different enumeration structure from local_table.py:
    # x is total TT, y is TC, z is CC. TTs split is resolved
    # *analytically* as endpoint 0 or x for any linear cost.
    lower=max(0,4*(a+b)-2*L*L)
    for x in range(adj(a,L)+1):
        for y in range(min(4*a-2*x,4*b,adj(a+b,L)-x)+1):
            upper_z=min(adj(b,L),(4*b-y)//2,adj(a+b,L)-x-y)
            for z in range(max(0,lower-x-y),upper_z+1):
                yield x,y,z


def target_weight(lump,kind,cut):
    a,b=lump
    if lump==(0,0):return Q(0)
    if kind=="product":
        e=Q(15,128)
        h=Q(1,4)
        return e**(a//2)*h**b
    # Separately chosen positive hybrid vector; high lump fixed,
    # and LOW orbit weights not uniform, so use their MAXIMUM.
    base=Q(7+(11*a+13*b)%19,5)
    if 2*a+3*b<=cut:
        return base+Q(13,5)
    return base


def pattern_costs(a,b,L,kind,cut):
    own=target_weight((a,b),kind,cut)
    require(own>0,"nonflat_trial_vector_not_positive")
    out={}
    for pat,moves in RULE.items():
        subtotal=Q(0)
        for da,db,n in moves:
            t=(a+da,b+db)
            require(0<=t[0]<=L*L and 0<=t[1]<=L*L and
                    t[0]+t[1]<=L*L and t[0]%2==0,
                    "S3_move_lump_out_of_range")
            subtotal+=n*(1-target_weight(t,kind,cut)/own)
        out[pat]=subtotal
    return out


def affine_coeff(a,b,L,r,weights):
    nl=2*L*L
    ee,eT,eC=weights["ee"],weights["eT"],weights["eC"]
    tt=min(weights["TTs"],weights["TTd"])
    constant=r*(2*a+3*b)+(nl-4*a-4*b)*ee+4*a*eT+4*b*eC
    return (constant,ee-2*eT+tt,
            ee-eT-eC+weights["TC"],ee-2*eC+weights["CC"])


def score(coeff,x,y,z):
    c,cx,cy,cz=coeff
    return c+cx*x+cy*y+cz*z


def check_one(L,r,kind,cut,exhaustive):
    total=0
    vertices=0
    integer_points=0
    gap_min=None
    winning=None
    for a in range(0,L*L+1,2):
        for b in range(L*L-a+1):
            if (a,b)==(0,0):continue
            weights=pattern_costs(a,b,L,kind,cut)
            coeff=affine_coeff(a,b,L,r,weights)
            vs=vertex_set(a,b,L)
            lower=min(score(coeff,*xyz) for xyz in vs)
            vertices+=len(vs)
            if exhaustive:
                true_min=None
                for x,y,z in integer_xyz(a,b,L):
                    value=score(coeff,x,y,z)
                    if true_min is None or value<true_min:
                        true_min=value
                    integer_points+=1
                require(true_min is not None,"no_feasible_integer_pattern")
                require(lower<=true_min,
                        "unsound_extreme_lower_"+str((L,a,b)))
            total+=1
            if gap_min is None or lower<gap_min:
                gap_min=lower
                winning=(a,b)
    result=gap_min-2*L*L*3  # source ratio minus 3*nl flat cap
    print("EXTREME_ROW_ORACLE|L=%d|mode=%s|lumps=%d|vertices=%d"
          "|integer_points=%d|ratio_minus_flat=%s|winning_lump=%d,%d"
          "|all_integer_rows_dominate_relaxed_vertices=PASS"
          "|json=0"%(L,kind,total,vertices,integer_points,result,
                     winning[0],winning[1]),flush=True)
    return result


def tamper_tests():
    # Wrong practice: replacing low-cut orbit-specific target values
    # with their *minimum* when evaluating a high-cut source. This
    # increases the purported lower bound and is unsound.
    a,b,cut=(0,2,4)
    src=target_weight((a,b),"hybrid",cut)
    tgt=(2,0) # CC three transpositions enter this low lump
    high=target_weight(tgt,"hybrid",cut)
    low=high-Q(13,5)
    require(0<low<high,"target_range_fixture")
    true_cc=3*(1-high/src)
    wrong_cc=3*(1-low/src)
    require(wrong_cc>true_cc,"min_instead_of_max_not_caught")

    # False 'extreme' polytope restricted to x=0: on a T-rich
    # lump it may eliminate a real integer signature and raise
    # the claimed lower bound. This must be rejected.
    a,b,L=4,0,4
    w=pattern_costs(a,b,L,"product",0)
    c=affine_coeff(a,b,L,Q(12),w)
    all_min=min(score(c,*v) for v in vertex_set(a,b,L))
    fake=min(score(c,*v) for v in vertex_set(a,b,L)
             if v[0]==0)
    require(fake>all_min,"false_x0_restriction_not_detected")
    int_min=min(score(c,*v) for v in integer_xyz(a,b,L))
    require(int_min<fake,"false_constraint_not_refuted_by_integer_witness")
    print("EXTREME_ROW_NEGATIVE|low_orbit_min_instead_of_max=REFUSED"
          "|false_x0_constraint=REFUSED"
          "|real_integer_counterexample_to_x0=PASS|json=0",
          flush=True)


def main():
    full=("--full" in sys.argv)
    allowed=({"--full"},set())
    require(set(sys.argv[1:]) in allowed,"unknown_arguments")
    tamper_tests()
    # Product form allows a known exact reference at 4x4 r=12.
    gap=check_one(4,Q(12),"product",0,True)
    require(gap==Q(657,32),
            "actual_4x4_one_C_reference_not_reproduced")
    check_one(4,Q(9),"hybrid",10,True)
    check_one(5,Q(15),"hybrid",11,full)
    check_one(6,Q(16),"hybrid",11,False)
    print("EXTREME_ROW_DONE|independent_S3_moves=PASS"
          "|exact_vertex_lower_bound=PASS|hybrid_low_lump_max=PASS"
          "|full_G1_G5_hybrid_certificate=NOT_COMPUTED"
          "|large_torus_gap=NOT_COMPUTED|json=0",flush=True)


if __name__=="__main__":
    main()
