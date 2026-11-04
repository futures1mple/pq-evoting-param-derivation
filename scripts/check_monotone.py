#!/usr/bin/env python3
"""Check the monotonicity of the stored hardness estimates (hypothesis (H4) of the paper).

  M-LWE: for fixed (d, L, model) the estimate is nonincreasing in q;
  M-SIS: for fixed (scheme, d, L, beta, model) nondecreasing in q, and for fixed (d, L, q, model)
         nonincreasing in beta.
All files of results/estimates/ and results/estimates_evolve/ are used, including the sub-bit moduli.
Prints the number of entries, of compared pairs and of violations.

usage: python3 scripts/check_monotone.py
"""
import glob
import json
import math
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def load(pattern):
    out = []
    for p in sorted(glob.glob(os.path.join(ROOT, "results", pattern))):
        with open(p) as f:
            out += [dict(json.loads(x), _file=os.path.relpath(p, ROOT)) for x in f if x.strip()]
    return out


def chains(groups, tol=1e-9):
    pairs = bad = 0
    viol = []
    for k, pts in groups.items():
        pts.sort()
        for (x1, v1), (x2, v2) in zip(pts, pts[1:]):
            if x1 == x2:
                continue
            pairs += 1
            if v2 > v1 + tol:
                bad += 1
                viol.append((k, x1, v1, x2, v2))
    return pairs, bad, viol


def main():
    lwe = [r for r in load("estimates/lwe*.jsonl") if r.get("min") is not None]
    sis = [r for r in load("estimates/sis*.jsonl") + load("estimates_evolve/sis.jsonl") if r.get("min") is not None]
    nulls = len([r for r in load("estimates/*.jsonl") + load("estimates_evolve/sis.jsonl")
                 if r.get("kind") in ("lwe", "sis") and r.get("min") is None])
    g = {}
    for r in lwe:
        g.setdefault((r["d"], r["L"], r["model"]), {})[r["q"]] = r["min"]      # duplicates: same value
    p1, b1, v1 = chains({k: list(v.items()) for k, v in g.items()})
    gq, gb = {}, {}
    for r in sis:
        lb = round(math.log2(r["beta_SIS"]), 9)
        gq.setdefault((r.get("scheme", "T-EVOLVE"), r["d"], r["L"], lb, r["model"]), {})[r["q"]] = -r["min"]
        gb.setdefault((r["d"], r["L"], r["q"], r["model"]), {})[lb] = r["min"]
    p2, b2, v2 = chains({k: list(v.items()) for k, v in gq.items()})
    p3, b3, v3 = chains({k: list(v.items()) for k, v in gb.items()})
    print("M-LWE entries %d (files: %s)" % (len(lwe), ", ".join(sorted({r["_file"] for r in lwe}))))
    print("M-SIS entries %d (files: %s)" % (len(sis), ", ".join(sorted({r["_file"] for r in sis}))))
    print("entries without a value: %d" % nulls)
    print("M-LWE nonincreasing in q: %d pairs, %d violations" % (p1, b1))
    print("M-SIS nondecreasing in q: %d pairs, %d violations" % (p2, b2))
    print("M-SIS nonincreasing in beta: %d pairs, %d violations" % (p3, b3))
    for v in (v1 + v2 + v3)[:20]:
        print("  violation", v)


if __name__ == "__main__":
    main()
