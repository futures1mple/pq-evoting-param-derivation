#!/usr/bin/env python3
"""Superseded by derive.py; kept for reference (selection of the first manuscript version).

Parameter selection restricted to the moduli q_b, the largest prime q = 17 (mod 32) below 2^b
(one modulus per bit length). derive.py searches all primes q = 17 (mod 32) and is used for the paper.

Reads results/estimates/{lwe,sis}.jsonl (written by estimate_grid.py) and, for every target level
lambda and set of cost models,
  1. computes b_hi(d): the largest bit length b <= 50 at which the M-LWE estimate is >= lambda in
     every model of the set (by (H4) the admissible bit lengths for hiding are b <= b_hi(d));
  2. for every configuration (N_V, n, L) of the tables selects the admissible pair (d, b) with the
     smallest d and, for it, the smallest b (admissible: q_b > beta_SIS, hiding and binding >= lambda
     in every model), and checks that no other admissible pair gives a smaller ballot;
  3. computes the largest electorate N_V*(d) served by each rank d (Corollary 8): the largest N_V
     with beta_SIS(d, N_V) < q_{b_hi(d)}, by bisection with param_report; if SageMath and the
     estimator are available (--binding), the M-SIS estimate at (N_V*(d), d, b_hi(d)) is added, to
     confirm that binding does not limit the threshold.

Writes results/selection.json and results/selection.md.

usage: python3 scripts/select_params.py [--binding]
"""
import argparse
import json
import math
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from estimate_grid import prime_below, BIN, EST, N  # noqa: E402

LAMBDAS = (128, 192)
SETS = {"MATZOV": ("MATZOV",), "MATZOV+MATZOV_q": ("MATZOV", "MATZOV_q")}
ALL_MODELS = ("MATZOV", "MATZOV_q", "ADPS16", "ADPS16_q")
N_LO, N_HI = 100, 10 ** 8


def report(nv, n, L, d, q):
    args = [BIN, str(nv), str(n), "1", str(L), "-1", str(d), str(q), "1", "64"]
    return json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)


def load(path):
    with open(path) as f:
        return [json.loads(line) for line in f]


_size = {}


def ballot_bytes(n, L, d, q):
    """Ballot size; independent of N_V (the ballot side of the chain does not involve N_V)."""
    k = (n, L, d, q)
    if k not in _size:
        _size[k] = report(10 ** 4, n, L, d, q)["ballot_bytes"]
    return _size[k]


_beta = {}


def beta_sis(nv, L, d):
    k = (nv, L, d)
    if k not in _beta:
        _beta[k] = report(nv, 4, L, d, prime_below(45))["beta_SIS"]
    return _beta[k]


def threshold(L, d, q):
    """Largest N_V in [N_LO, N_HI] with beta_SIS(d, N_V) < q (beta_SIS nondecreasing in N_V, Lemma 3)."""
    if beta_sis(N_LO, L, d) >= q:
        return 0
    if beta_sis(N_HI, L, d) < q:
        return N_HI
    lo, hi = N_LO, N_HI
    while hi - lo > 1 and hi / lo > 1.002:
        mid = min(hi - 1, max(lo + 1, int(round(math.sqrt(lo * hi)))))
        if beta_sis(mid, L, d) < q:
            lo = mid
        else:
            hi = mid
    return lo


def sis_estimate(d, L, q, beta, models):
    sys.path.insert(0, EST)
    from sage.all import log, oo
    from estimator import SIS, RC
    from estimator.reduction import MATZOV
    m = {"MATZOV": RC.MATZOV, "MATZOV_q": MATZOV(nn="quantum")}
    out = {}
    for name in models:
        res = SIS.estimate(SIS.Parameters(n=d * N, q=q, length_bound=beta, m=(2 * d + L) * N, norm=2),
                           red_cost_model=m[name], quiet=True)
        out[name] = min(float(log(v["rop"], 2)) for v in res.values() if v["rop"] != oo)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binding", action="store_true", help="estimate M-SIS at the thresholds (needs Sage)")
    a = ap.parse_args()
    est = os.path.join(ROOT, "results", "estimates")
    lwe, sis = load(os.path.join(est, "lwe.jsonl")), load(os.path.join(est, "sis.jsonl"))
    H = {(r["d"], r["b"], r["model"]): r["min"] for r in lwe}
    B = {(r["NV"], r["n"], r["L"], r["d"], r["b"], r["model"]): r["min"] for r in sis}
    beta_tab = {(r["NV"], r["n"], r["L"], r["d"]): r["beta_SIS"] for r in sis}
    configs = sorted({k[:3] for k in B})
    ds = sorted({r["d"] for r in lwe})
    bs = sorted({r["b"] for r in lwe})
    out = {"tables": "results/estimates (lattice estimator 53da598)", "targets": []}
    md = ["# Parameter selection", "",
          "From `results/estimates/{lwe,sis}.jsonl` by `scripts/select_params.py`. Ballot: T-EVOLVE, free weight, signed second challenge; sizes in KiB (1024 bytes).", ""]
    for lam in LAMBDAS:
        for sname, ms in SETS.items():
            t = {"lambda": lam, "models": list(ms)}
            bhi = {d: max([b for b in bs if all(H[(d, b, m)] >= lam for m in ms)], default=None) for d in ds}
            t["b_hi"] = bhi
            rows = []
            for (nv, n, L) in configs:
                adm = [(d, b) for d in ds for b in bs
                       if all((nv, n, L, d, b, m) in B for m in ms)
                       and all(H[(d, b, m)] >= lam and B[(nv, n, L, d, b, m)] >= lam for m in ms)]
                if not adm:
                    rows.append({"NV": nv, "n": n, "L": L, "d": None})
                    continue
                d, b = min(adm)
                q = prime_below(b)
                sizes = {p: ballot_bytes(n, L, p[0], prime_below(p[1])) for p in adm if p[0] <= d + 1}
                smallest = min(sizes, key=sizes.get)
                b_lo = min(bb for bb in bs if (nv, n, L, d, bb, ms[0]) in B)
                rows.append({"NV": nv, "n": n, "L": L, "d": d, "b": b, "q": q, "beta_SIS": beta_tab[(nv, n, L, d)],
                             "b_lo": b_lo, "binding_limits": b > b_lo,
                             "hiding": {m: H[(d, b, m)] for m in ALL_MODELS},
                             "binding": {m: B[(nv, n, L, d, b, m)] for m in ALL_MODELS},
                             "ballot_bytes": sizes[(d, b)], "min_size_pair": list(smallest),
                             "lexicographic_is_min_size": smallest == (d, b)})
            t["selection"] = rows
            thr = []
            for L in (1, 2):
                for d in ds:
                    if bhi[d] is None:
                        continue
                    q = prime_below(bhi[d])
                    nstar = threshold(L, d, q)
                    rec = {"L": L, "d": d, "b_hi": bhi[d], "q": q, "N_V_star": nstar,
                           "beyond_range": nstar == N_HI}
                    if a.binding and 0 < nstar < N_HI:
                        rec["binding_at_threshold"] = sis_estimate(d, L, q, beta_sis(nstar, L, d), ms)
                    thr.append(rec)
                    print("lambda=%d %s L=%d d=%d b_hi=%d N_V*=%s" % (lam, sname, L, d, bhi[d],
                          (">= 1e8" if nstar == N_HI else "%.4g" % nstar)), flush=True)
                    if nstar == N_HI:
                        break            # every larger rank also serves the whole range
            t["thresholds"] = thr
            out["targets"].append(t)
            md += ["## lambda = %d, cost models: %s" % (lam, ", ".join(ms)), "",
                   "b_hi(d): " + ", ".join("%d: %s" % (d, bhi[d]) for d in ds), "",
                   "| N_V | n | L | d | log2 q | log2 beta_SIS | hiding | binding | ballot (KiB) |",
                   "|---|---|---|---|---|---|---|---|---|"]
            for r in rows:
                if r["d"] is None:
                    md.append("| %d | %d | %d | none | | | | | |" % (r["NV"], r["n"], r["L"]))
                    continue
                md.append("| %d | %d | %d | %d | %d | %.2f | %s | %s | %.1f |" % (
                    r["NV"], r["n"], r["L"], r["d"], r["b"], math.log2(r["beta_SIS"]),
                    " / ".join("%.1f" % r["hiding"][m] for m in ms), " / ".join("%.1f" % r["binding"][m] for m in ms),
                    r["ballot_bytes"] / 1024))
            md += ["", "Largest electorate per rank (L = 1 / L = 2):", ""]
            for d in ds:
                xs = [x for x in thr if x["d"] == d]
                if not xs:
                    continue
                md.append("- d = %d (b_hi = %d): %s" % (d, xs[0]["b_hi"], " / ".join(
                    ">= 10^8" if x["beyond_range"] else ("none" if x["N_V_star"] == 0 else "%.3g" % x["N_V_star"]) for x in xs)))
            md.append("")
    with open(os.path.join(ROOT, "results", "selection.json"), "w") as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(ROOT, "results", "selection.md"), "w") as f:
        f.write("\n".join(md) + "\n")


if __name__ == "__main__":
    main()
