#!/usr/bin/env python3
"""Hardness tables for the parameter search (needs SageMath/passagemath and the lattice estimator).

The search of the method needs, for every candidate module rank d and bit length b of the modulus
q_b (largest prime q = 17 mod 32 below 2^b):
  * the M-LWE estimate (hiding), which depends on (d, L, q_b) only, and
  * the M-SIS estimate (binding), which depends on (d, L, q_b, beta_SIS), where beta_SIS comes from
    param_report for the electorate N_V and the number of authorities n.
This script computes both tables under four cost models and writes every attack of every estimate
as one JSON line, so that the selection (any target level, any set of cost models) and the
monotonicity checks can be done afterwards without the estimator.

  lwe  M-LWE table: d in D_RANGE, b in B_LWE, L = 1 (secret and error D_1, m = (d+L)N samples;
       arora-gb and bkw skipped, as in t-evolve-cryptoalg tools/quantum_models.py)
  sis  M-SIS table: for every (N_V, n, L) in CONFIGS and d in D_RANGE, b from the smallest bit length
       with q_b > beta_SIS up to 50 (Euclidean norm, m = (2d+L)N columns)
  all  both, SIS first

Cost models: MATZOV (classical), MATZOV_q (quantum nearest-neighbour cost of AGPS20), ADPS16
(core-SVP 2^(0.292 beta)), ADPS16_q (2^(0.265 beta)).

The runs can be interrupted and restarted: finished entries (same key) are skipped.

usage: python3 scripts/estimate_grid.py {lwe,sis,all} [--jobs J] [--quick]
       --quick: a few entries only, to test the installation (written to results/estimates_quick/)
Environment: ESTIMATOR_DIR (default deps/lattice-estimator).
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BIN = os.path.join(ROOT, "deps", "t-evolve-national-scale", "build", "param_report")
EST = os.environ.get("ESTIMATOR_DIR", os.path.join(ROOT, "deps", "lattice-estimator"))
N = 256
MODELS = ("MATZOV", "MATZOV_q", "ADPS16", "ADPS16_q")
D_RANGE = range(6, 13)
B_LWE = range(36, 51)
B_MAX = 50                                   # q < 2^50 in the implementation
CONFIGS = [(nv, 4, L) for nv in (10 ** 4, 10 ** 5, 10 ** 6, 6480000, 10 ** 7) for L in (1, 2)] + \
          [(10 ** 6, 9, 1)]                  # n = 9: check that n does not matter


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_below(b):
    """Largest prime q = 17 (mod 32) below 2^b."""
    q = 2 ** b - 15
    while not is_prime(q):
        q -= 32
    return q


def param_report(nv, n, L, d):
    q = prime_below(45)                      # any valid q: the bounds do not depend on q
    args = [BIN, str(nv), str(n), "1", str(L), "-1", str(d), str(q), "1", "64"]
    return json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)


def estimate(job):
    """Runs in a fresh worker process: one estimate under one cost model."""
    sys.path.insert(0, EST)
    from sage.all import log, oo
    from estimator import LWE, SIS, ND, RC
    from estimator.reduction import ADPS16, MATZOV
    models = {"MATZOV": RC.MATZOV, "MATZOV_q": MATZOV(nn="quantum"), "ADPS16": RC.ADPS16,
              "ADPS16_q": ADPS16(mode="quantum")}
    t0 = time.time()
    d, L, q, model = job["d"], job["L"], job["q"], job["model"]
    if job["kind"] == "lwe":
        p = LWE.Parameters(n=d * N, q=q, Xs=ND.DiscreteGaussian(1.0), Xe=ND.DiscreteGaussian(1.0), m=(d + L) * N)
        res = LWE.estimate(p, red_cost_model=models[model], deny_list=("arora-gb", "bkw"), quiet=True)
    else:
        p = SIS.Parameters(n=d * N, q=q, length_bound=job["beta_SIS"], m=(2 * d + L) * N, norm=2)
        res = SIS.estimate(p, red_cost_model=models[model], quiet=True)
    attacks = {k: (float(log(v["rop"], 2)) if v["rop"] != oo else None) for k, v in res.items()}
    blocks = {k: int(v["beta"]) for k, v in res.items() if "beta" in v and v["rop"] != oo}
    finite = [v for v in attacks.values() if v is not None]
    out = dict(job, attacks=attacks, block_sizes=blocks, min=min(finite) if finite else None,
               seconds=round(time.time() - t0, 1), estimator="53da598")
    return out


def key(rec):
    if rec["kind"] == "lwe":
        return ("lwe", rec["d"], rec["L"], rec["b"], rec["model"])
    return ("sis", rec["NV"], rec["n"], rec["L"], rec["d"], rec["b"], rec["model"])


def done_keys(path):
    keys = set()
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                try:
                    keys.add(key(json.loads(line)))
                except (ValueError, KeyError):
                    pass
    return keys


def lwe_jobs(quick):
    ds, bs, models = ((7,), (45,), ("MATZOV",)) if quick else (D_RANGE, B_LWE, MODELS)
    return [{"kind": "lwe", "d": d, "L": 1, "b": b, "q": prime_below(b), "model": m}
            for d in ds for b in bs for m in models]


def sis_jobs(quick, outdir):
    configs, ds, models = ((CONFIGS[0],), (7,), ("MATZOV",)) if quick else (CONFIGS, D_RANGE, MODELS)
    jobs = []
    with open(os.path.join(outdir, "bounds.jsonl"), "w") as fb:
        for (nv, n, L) in configs:
            for d in ds:
                r = param_report(nv, n, L, d)
                fb.write(json.dumps(dict(r, NV=nv)) + "\n")
                fb.flush()
                if r["sigma1"] <= 0 or r["sigma2"] <= 0 or r["sigma_J"] <= 0:
                    continue                 # no sampler plan: infeasible, nothing to estimate
                beta = r["beta_SIS"]
                b = int(math.log2(beta)) + 1
                while prime_below(b) <= beta:
                    b += 1
                for bb in range(b, (b + 1 if quick else B_MAX + 1)):
                    q = prime_below(bb)
                    for m in models:
                        jobs.append({"kind": "sis", "NV": nv, "n": n, "L": L, "d": d, "b": bb, "q": q,
                                     "beta_SIS": beta, "model": m})
    return jobs


def run(jobs, path, nproc):
    done = done_keys(path)
    todo = [j for j in jobs if key(j) not in done]
    print("%s: %d entries, %d already done, %d to run on %d processes" % (os.path.basename(path), len(jobs),
          len(jobs) - len(todo), len(todo), nproc), flush=True)
    t0 = time.time()
    ctx = mp.get_context("spawn")
    with ctx.Pool(nproc, maxtasksperchild=1) as pool, open(path, "a") as f:
        for i, rec in enumerate(pool.imap_unordered(estimate, todo), 1):
            f.write(json.dumps(rec) + "\n")
            f.flush()
            print("[%d/%d, %.0f min] %s min=%s (%.0f s)" % (i, len(todo), (time.time() - t0) / 60, key(rec),
                  None if rec["min"] is None else round(rec["min"], 2), rec["seconds"]), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=("lwe", "sis", "all"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    outdir = os.path.join(ROOT, "results", "estimates_quick" if a.quick else "estimates")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "run_info.jsonl"), "a") as f:
        f.write(json.dumps({"start": time.strftime("%Y-%m-%d %H:%M:%S"), "what": a.what, "jobs": a.jobs,
                            "quick": a.quick, "python": sys.version.split()[0], "cpus": os.cpu_count(),
                            "estimator_dir": EST}) + "\n")
    if a.what in ("sis", "all"):
        run(sis_jobs(a.quick, outdir), os.path.join(outdir, "sis.jsonl"), a.jobs)
    if a.what in ("lwe", "all"):
        run(lwe_jobs(a.quick), os.path.join(outdir, "lwe.jsonl"), a.jobs)
    print("finished", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    main()
