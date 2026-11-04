#!/usr/bin/env python3
"""The method applied to EVOLVE (del Pino, Lyubashevsky, Neven, Seiler, ACM CCS 2017).

EVOLVE's own parameter analysis is used as published:
  * Section 5 lists the constraints on the proof bounds (n = ring degree = 256, d = module rank,
    sigma = 1, N_A authorities, bucket size l = 30):
        B_OR      = 2 N_A sqrt(n(2d+1)) sigma          B'_OR     = 44 sqrt(60 n(2d+1)) B_OR
        B_Amo,1   = 2 sqrt(n(2d+1)) sigma              B'_Amo,1  = 2684 n sqrt(2d+1) B_Amo,1
        B_Amo,2   = 2 (l+1) sqrt(n(2d+1)) sigma        B'_Amo,2  = 2684 n sqrt(2d+1) B_Amo,2
    and the consistency conditions  2 B'_OR <= B_r,  2 sqrt(60) N_A B'_Amo,1 <= B_r,
    (l+1) B'_Amo,1 <= B_r,  B'_Amo,2 <= B_r  (the square root in (4) covers 60 only);
  * Section 3.1: the commitment is binding under M-SIS_{q, n, d, 2d+1, 2 B_r}.
So beta_SIS = 2 B_r with B_r the largest of the four left-hand sides. The bounds do not involve the
number of voters N_V (the aggregation tree of Section 4.2 keeps B_r independent of N_V; each
amortized proof covers at least about 2800 commitments, Section 3.3), so the method has N_A, d and
the target as inputs. Commitment and hiding instance have the same shape as T-EVOLVE with L = 1
(width 2d+1, M-LWE with n = dN, m = (d+1)N), so the M-LWE table of estimate_grid.py is reused; this
script adds the entries for d = 5 and the M-SIS table for EVOLVE's bound.

  estimate  M-LWE for d = 5, b = 34..37 and d = 6..12, b = 35 (added to results/estimates/lwe.jsonl) and M-SIS for
            N_A in NA_LIST, d = 5..12, b from the smallest q_b > beta_SIS to +6 bits
            (results/estimates_evolve/sis.jsonl); needs SageMath and the estimator
  select    selection for lambda in {128, 192} and the cost model sets of select_params.py
            (results/evolve_selection.{json,md}); no estimator needed

usage: python3 scripts/evolve_chain.py {estimate,select} [--jobs J]
"""
import argparse
import json
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from estimate_grid import prime_below, run, N  # noqa: E402

NA_LIST = (2, 3, 4, 5, 7, 9, 15)
D_RANGE = range(5, 13)
LAMBDAS = (128, 192)
SETS = {"MATZOV": ("MATZOV",), "MATZOV+MATZOV_q": ("MATZOV", "MATZOV_q")}
MODELS = ("MATZOV", "MATZOV_q", "ADPS16", "ADPS16_q")
L_BUCKET, SIGMA = 30, 1.0


def evolve_bounds(d, NA, l=L_BUCKET, sigma=SIGMA):
    w = math.sqrt(N * (2 * d + 1))
    B_OR = 2 * NA * w * sigma
    B_ORp = 44 * math.sqrt(60 * N * (2 * d + 1)) * B_OR
    B_A1 = 2 * w * sigma
    B_A1p = 2684 * N * math.sqrt(2 * d + 1) * B_A1
    B_A2 = 2 * (l + 1) * w * sigma
    B_A2p = 2684 * N * math.sqrt(2 * d + 1) * B_A2
    cons = {"(3) 2 B'_OR": 2 * B_ORp, "(4) 2 sqrt(60) N_A B'_Amo,1": 2 * math.sqrt(60) * NA * B_A1p,
            "(7) (l+1) B'_Amo,1": (l + 1) * B_A1p, "(8) B'_Amo,2": B_A2p}
    dom = max(cons, key=cons.get)
    B_r = cons[dom]
    return {"d": d, "N_A": NA, "B_r": B_r, "beta_SIS": 2 * B_r, "log2beta": math.log2(2 * B_r),
            "dominant": dom, "constraints_log2": {k: math.log2(v) for k, v in cons.items()}}


def first_b(beta):
    b = int(math.log2(beta)) + 1
    while prime_below(b) <= beta:
        b += 1
    return b


def estimate(jobs):
    est = os.path.join(ROOT, "results", "estimates")
    # EVOLVE needs smaller moduli than T-EVOLVE: the smallest admissible bit length is 34 for d = 5
    # and 35 for d >= 6 (N_A = 2); the table of estimate_grid.py starts at b = 36
    lwe = [{"kind": "lwe", "d": 5, "L": 1, "b": b, "q": prime_below(b), "model": m}
           for b in range(34, 38) for m in MODELS] + \
          [{"kind": "lwe", "d": d, "L": 1, "b": 35, "q": prime_below(35), "model": m}
           for d in range(6, 13) for m in MODELS]
    run(lwe, os.path.join(est, "lwe.jsonl"), jobs)
    out = os.path.join(ROOT, "results", "estimates_evolve")
    os.makedirs(out, exist_ok=True)
    sis = []
    for NA in NA_LIST:
        for d in D_RANGE:
            e = evolve_bounds(d, NA)
            b0 = first_b(e["beta_SIS"])
            for b in range(b0, min(50, b0 + 6) + 1):
                for m in MODELS:
                    # NV = 0 and n = N_A mark EVOLVE entries; L = 1 (one message row)
                    sis.append({"kind": "sis", "scheme": "EVOLVE", "NV": 0, "n": NA, "L": 1, "d": d, "b": b,
                                "q": prime_below(b), "beta_SIS": e["beta_SIS"], "model": m})
    run(sis, os.path.join(out, "sis.jsonl"), jobs)


def load(path):
    with open(path) as f:
        return [json.loads(line) for line in f]


def select():
    H = {(r["d"], r["b"], r["model"]): r["min"] for r in load(os.path.join(ROOT, "results", "estimates", "lwe.jsonl"))}
    S = {(r["n"], r["d"], r["b"], r["model"]): r["min"]
         for r in load(os.path.join(ROOT, "results", "estimates_evolve", "sis.jsonl"))}
    out = {"scheme": "EVOLVE (CCS 2017), constraints of its Section 5, binding at 2 B_r (Section 3.1)",
           "bounds": [evolve_bounds(d, NA) for NA in NA_LIST for d in D_RANGE], "targets": []}
    md = ["# EVOLVE through the method", "",
          "Constraints of EVOLVE Section 5, binding at 2 B_r (Section 3.1); independent of the number of voters.",
          "Commitment = (d+1) ring elements mod q (one share of one voter, per authority); KiB = 1024 bytes.", ""]
    md += ["| N_A | log2 beta_SIS (d = 6) | dominant constraint |", "|---|---|---|"]
    for NA in NA_LIST:
        e = evolve_bounds(6, NA)
        md.append("| %d | %.2f | %s |" % (NA, e["log2beta"], e["dominant"]))
    md.append("")
    for lam in LAMBDAS:
        for sname, ms in SETS.items():
            rows = []
            for NA in NA_LIST:
                adm = [(d, b) for d in D_RANGE for b in range(34, 51)
                       if all((d, b, m) in H and (NA, d, b, m) in S for m in ms)
                       and all(H[(d, b, m)] >= lam and S[(NA, d, b, m)] >= lam for m in ms)]
                if not adm:
                    rows.append({"N_A": NA, "d": None})
                    continue
                d, b = min(adm)
                e = evolve_bounds(d, NA)
                rows.append({"N_A": NA, "d": d, "b": b, "q": prime_below(b), "log2beta": e["log2beta"],
                             "b_lo": first_b(e["beta_SIS"]), "binding_limits": b > first_b(e["beta_SIS"]),
                             "hiding": {m: H[(d, b, m)] for m in MODELS},
                             "binding": {m: S[(NA, d, b, m)] for m in MODELS},
                             "commitment_bytes": (d + 1) * N * b / 8})
            out["targets"].append({"lambda": lam, "models": list(ms), "selection": rows})
            md += ["## lambda = %d, cost models: %s" % (lam, ", ".join(ms)), "",
                   "| N_A | d | log2 q | log2 beta_SIS | hiding | binding | commitment (KiB) |",
                   "|---|---|---|---|---|---|---|"]
            for r in rows:
                if r["d"] is None:
                    md.append("| %d | none | | | | | |" % r["N_A"])
                    continue
                md.append("| %d | %d | %d | %.2f | %s | %s | %.2f |" % (
                    r["N_A"], r["d"], r["b"], r["log2beta"], " / ".join("%.1f" % r["hiding"][m] for m in ms),
                    " / ".join("%.1f" % r["binding"][m] for m in ms), r["commitment_bytes"] / 1024))
            md.append("")
    with open(os.path.join(ROOT, "results", "evolve_selection.json"), "w") as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(ROOT, "results", "evolve_selection.md"), "w") as f:
        f.write("\n".join(md) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=("estimate", "select"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    if a.what == "estimate":
        estimate(a.jobs)
    else:
        select()


if __name__ == "__main__":
    main()
