#!/usr/bin/env python3
"""Binding at the rank thresholds (needs Sage for the estimates).

The thresholds N_V*(d) of results/derivation.json are computed from the hiding condition
q_lo(d; N_V) <= q_hi(d) alone. At N_V*(d) the bound beta_SIS is largest for rank d, so the M-SIS
estimate at (d, q_lo, beta_SIS(N_V*)) is the smallest binding value of rank d up to its threshold.
This script estimates it for every finite threshold (L = 1, n = 4) in every model of the target, and
also checks the feasibility conditions S5(i) and S5(iv) there.
Writes results/estimates/sis_thr.jsonl and results/threshold_binding.md.

usage: python3 scripts/threshold_binding.py [--jobs J]
"""
import argparse
import json
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from estimate_grid import run  # noqa: E402
from derive import TEVOLVE, TARGETS, prime_above, EST  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    with open(os.path.join(ROOT, "results", "derivation.json")) as f:
        der = json.load(f)
    pts, jobs = [], []
    seen = set()
    for tname, t in der["targets"].items():
        lam_b, mb = TARGETS[tname][2], TARGETS[tname][3]
        for x in t["thresholds"]:
            ns = x["N_V_star"]
            if x["L"] != 1 or not ns or ns >= 10 ** 8:
                continue
            cfg, d = (ns, 4, 1), x["d"]
            beta = TEVOLVE.beta(cfg, d)
            q = prime_above(beta)
            ok, why = TEVOLVE.feasible(cfg, d)
            pts.append({"target": tname, "d": d, "N_V": ns, "q": q, "beta_SIS": beta, "lam_b": lam_b, "models": mb,
                        "feasible": ok, "log2_budget": TEVOLVE.budget(cfg, d),
                        "log2_sigma_max": math.log2(max(TEVOLVE.report(cfg, d)[k] for k in ("sigma_J", "sigma1", "sigma2")))})
            for m in mb:
                if (d, q, m) in seen:
                    continue
                seen.add((d, q, m))
                jobs.append({"kind": "sis", "scheme": "T-EVOLVE", "NV": ns, "n": 4, "L": 1, "d": d,
                             "b": round(math.log2(q), 6), "q": q, "beta_SIS": beta, "model": m})
    path = os.path.join(EST, "sis_thr.jsonl")
    run(jobs, path, a.jobs)
    with open(path) as f:
        est = [json.loads(x) for x in f]
    md = ["# Binding at the rank thresholds (L = 1, n = 4)", "", "`scripts/threshold_binding.py`", "",
          "| target | d | N_V* | log2 q | log2 beta | binding (min over models) | margin | S5 feasible | log2 budget | log2 max sigma |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for p in pts:
        v = {r["model"]: r["min"] for r in est if r["d"] == p["d"] and r["q"] == p["q"] and r["model"] in p["models"]
             and abs(math.log2(r["beta_SIS"]) - math.log2(p["beta_SIS"])) < 1e-9}
        b = min(v.values()) if len(v) == len(p["models"]) else None
        p["binding"] = b
        md.append("| %s | %d | %d | %.3f | %.3f | %s | %s | %s | %.2f | %.2f |" % (
            p["target"], p["d"], p["N_V"], math.log2(p["q"]), math.log2(p["beta_SIS"]),
            "-" if b is None else "%.1f" % b, "-" if b is None else "%.1f" % (b - p["lam_b"]),
            p["feasible"], p["log2_budget"], p["log2_sigma_max"]))
    with open(os.path.join(ROOT, "results", "threshold_binding.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(ROOT, "results", "threshold_binding.json"), "w") as f:
        json.dump(pts, f, indent=1, default=list)
    print("\n".join(md))


if __name__ == "__main__":
    main()
