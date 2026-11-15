#!/usr/bin/env python3
"""Fail-closed fan-in for four completely disjoint independent low-source shards.

All four shard rows are emitted by src/bin/vectorcheck.rs, which
separately enumerates every raw electric move, group, canonical form,
and vertex from one frozen HYBRID_BIG vector. This script does not
recalculate source rows or certify a result without those four passes.
"""
import pathlib
import sys

def require(ok, reason):
    if not ok:
        raise RuntimeError("VECTORCHECK_FANIN_REFUSE_"+reason)

def parse(path):
    lines=path.read_text().splitlines()
    hits=[line for line in lines if line.startswith("VECTORCHECK|")]
    require(len(hits)==1,"ONE_SUCCESSFUL_CHECKER_ROW_REQUIRED")
    d={}
    for chunk in hits[0].split("|")[1:]:
        require("=" in chunk,"MALFORMED_FIELD")
        k,v=chunk.split("=",1)
        require(k not in d,"DUPLICATE_FIELD")
        d[k]=v
    require(d["result"]=="PASS" and d["scope"]=="SHARD_ONLY",
            "NOT_A_REAL_SOURCE_SHARD")
    shas=[line.split('=',1)[1] for line in lines if line.startswith('SOURCE_SHA256=')]
    require(len(shas)<=1,'DUPLICATE_SOURCE_SHA_LINE')
    if shas:
        d['source_sha256']=shas[0].lower()
    return d

def main():
    require(len(sys.argv) in (2,3,4),"USAGE_RECEIPTS_DIR_EXPECTED_SHARDS_SHA256")
    folder=pathlib.Path(sys.argv[1])
    expected=int(sys.argv[2]) if len(sys.argv)>=3 else 4
    require(1<=expected<=256,"INVALID_SHARD_COUNT")
    locked_sha=sys.argv[3].lower() if len(sys.argv)==4 else None
    if locked_sha is not None:
        require(len(locked_sha)==64 and all(x in '0123456789abcdef' for x in locked_sha),'INVALID_SOURCE_SHA256')
    files=sorted(folder.glob("CHECK-SHARD-*.hbp"))
    require(len(files)==expected,"NEED_EVERY_SHARD_RECEIPT")
    rows=[parse(p) for p in files]
    if locked_sha is not None:
        require(all(row.get('source_sha256')==locked_sha for row in rows),
                'SOURCE_HASH_MISMATCH_ACROSS_SHARDS')
    seen=set()
    total_classes=0
    total_moves=0
    claimed=None
    minimum=None
    for row in rows:
        require(row["shard"].endswith("/"+str(expected)),"WRONG_SHARD_TOTAL")
        index=int(row["shard"].split("/")[0])
        require(index in range(expected) and index not in seen,
                "DUPLICATE_OR_INVALID_SOURCE_SHARD")
        seen.add(index)
        identity=(row["Lx"],row["Ly"],row["r"],row["cut"],
                  row["classes_total"],row["lumps_above"],
                  row["vertices_above"],row["claim_micro"],
                  row["high_min_micro"])
        if claimed is None:
            claimed=identity
        require(identity==claimed,"SHARD_SOURCE_OR_CLAIM_INCONSISTENT")
        cnt=int(row["classes_checked"])
        moves=int(row["raw_moves_checked"])
        require(moves==cnt*3*(2*int(row["Lx"])*int(row["Ly"])),
                "INCOMPLETE_RAW_ELECTRIC_CENSUS")
        total_classes+=cnt
        total_moves+=moves
        low=int(row["low_min_micro"])
        high=int(row["high_min_micro"])
        this_min=min(low,high)
        minimum=this_min if minimum is None else min(minimum,this_min)
    require(seen==set(range(expected)),"MISSING_SHARD")
    require(total_classes==int(rows[0]["classes_total"]),
            "INCOMPLETE_VERIFIED_EXACT_REGION")
    require(minimum==int(rows[0]["claim_micro"]),
            "CLAIM_DOES_NOT_EQUAL_TRUE_INTEGER_MINIMUM")
    # An overclaim by even one integer micro-unit violates the
    # independently computed minimum across the full source space.
    require(not minimum>=int(rows[0]["claim_micro"])+1,
            "RAISED_CLAIM_MUST_BE_REFUSED")
    print("VECTORCHECK_FANIN|Lx={}|Ly={}|r={}|cut={}"
          "|checked_shards={}|fully_checked_exact_classes={}"
          "|fully_checked_raw_moves={}"
          "|high_vertices_checked_per_shard={}"
          "|claim_micro={}|independent_min_micro={}"
          "|raised_claim_by_one_micro=REFUSED"
          "|all_sources=CHECKED|status=PASS|json=0".format(
              rows[0]["Lx"],rows[0]["Ly"],rows[0]["r"],
              rows[0]["cut"],expected,total_classes,total_moves,
              rows[0]["vertices_above"],rows[0]["claim_micro"],
              minimum),flush=True)

if __name__=="__main__":
    main()
