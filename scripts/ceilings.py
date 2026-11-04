#!/usr/bin/env python3
"""Hiding ceilings at sub-bit resolution (needs SageMath and the estimator).

The integer table results/estimates/lwe.jsonl gives the M-LWE estimate at q_b (largest prime
q = 17 mod 32 below 2^b). The modulus need not be of that form: any prime q = 17 (mod 32) is
admissible, and a modulus slightly above 2^b costs the same number of bits per coefficient as q_{b+1}.
For every rank d, cost model and target lambda in {128, 192} whose crossing lies inside the table,
this script interpolates the crossing x* in log2 q between the two neighbouring table entries and
estimates M-LWE at the largest prime q = 17 (mod 32) below 2^(x* - 0.05) and below 2^(x* + 0.05).
With the monotonicity of the estimate in q (hypothesis (H4) of the paper, Lemma 5; checked by check_monotone.py) the ceiling then lies
between the two values, at 0.1-bit resolution.

Writes results/estimates/lwe_fine.jsonl (same format as lwe.jsonl, with the field "x" = log2 of the
power of two below which q was taken; "b" holds x rounded to 3 decimals).

With --refine, the script then closes every remaining gap: for each (d, model, lambda), if the largest
computed modulus that reaches lambda and the smallest larger computed modulus that fails are more than
0.1 bit apart (L = 1, all stored M-LWE files), it estimates at the largest prime q = 17 (mod 32) below
the midpoint in log2 q, and repeats until every gap is at most 0.1 bit.

With --refine-dual-free, the same is done for the MATZOV estimate without the dual attacks.

usage: python3 scripts/ceilings.py [--jobs J] [--refine | --refine-dual-free]
"""
import argparse
import json
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from estimate_grid import is_prime, run  # noqa: E402

MODELS = ("MATZOV", "MATZOV_q", "ADPS16", "ADPS16_q")
LAMBDAS = (128, 192)
STEP = 0.05


def prime_below_real(x):
    """Largest prime q = 17 (mod 32) below 2^x (x real)."""
    q = int(math.floor(2 ** x))
    q -= (q - 17) % 32
    if q >= 2 ** x:
        q -= 32
    while not is_prime(q):
        q -= 32
    return q


def crossings(lwe):
    H = {}
    for r in lwe:
        H.setdefault((r["d"], r["model"]), {})[r["b"]] = r["min"]
    out = []
    for (d, m), tab in sorted(H.items()):
        bs = sorted(tab)
        for lam in LAMBDAS:
            for b1, b2 in zip(bs, bs[1:]):
                if b2 == b1 + 1 and tab[b1] >= lam > tab[b2]:
                    x = b1 + (tab[b1] - lam) / (tab[b1] - tab[b2])
                    out.append({"d": d, "model": m, "lambda": lam, "b_int": b1, "x_interp": x})
    return out


def load_l1(est):
    out = []
    for name in ("lwe.jsonl", "lwe_fine.jsonl", "lwe_sel.jsonl", "lwe_naive.jsonl"):
        path = os.path.join(est, name)
        if os.path.exists(path):
            with open(path) as f:
                out += [json.loads(line) for line in f if line.strip()]
    return [r for r in out if r["L"] == 1 and r["min"] is not None]


def no_dual(r):
    return min(v for k, v in r["attacks"].items() if v is not None and not k.startswith("dual"))


def gaps(rows, width=0.1, value=None, models=None):
    """(d, model, lambda, log2 q_pass, log2 q_fail) for every gap wider than `width` bits."""
    H = {}
    for r in rows:
        if models is None or r["model"] in models:
            H.setdefault((r["d"], r["model"]), {})[r["q"]] = r["min"] if value is None else value(r)
    out = []
    for (d, m), tab in sorted(H.items()):
        for lam in LAMBDAS:
            ok = [q for q, v in tab.items() if v >= lam]
            if not ok:
                continue
            qp = max(ok)
            bad = [q for q, v in tab.items() if v < lam and q > qp]
            if bad and math.log2(min(bad)) - math.log2(qp) > width + 1e-3:
                out.append((d, m, lam, math.log2(qp), math.log2(min(bad))))
    return out


def refine(est, nproc, dual_free=False):
    while True:
        g = gaps(load_l1(est), value=no_dual if dual_free else None, models=("MATZOV",) if dual_free else None)
        jobs, seen = [], set()
        for d, m, lam, lp, lf in g:
            q = prime_below_real((lp + lf) / 2)
            if (d, q, m) in seen or q <= 2 ** lp:
                continue
            seen.add((d, q, m))
            x = round(math.log2(q), 3)
            jobs.append({"kind": "lwe", "d": d, "L": 1, "b": x, "x": x, "q": q, "model": m})
        print("gaps wider than 0.1 bit: %d" % len(g), flush=True)
        if not jobs:
            return
        run(jobs, os.path.join(est, "lwe_fine.jsonl"), nproc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--refine", action="store_true")
    ap.add_argument("--refine-dual-free", action="store_true",
                    help="the same for the MATZOV estimate without the dual attacks (sensitivity check)")
    a = ap.parse_args()
    est = os.path.join(ROOT, "results", "estimates")
    if a.refine or a.refine_dual_free:
        refine(est, a.jobs, dual_free=a.refine_dual_free)
        return
    with open(os.path.join(est, "lwe.jsonl")) as f:
        lwe = [json.loads(line) for line in f]
    jobs, seen = [], set()
    for c in crossings(lwe):
        for x in (c["x_interp"] - STEP, c["x_interp"] + STEP):
            x = round(x, 3)
            if x <= c["b_int"] or x >= c["b_int"] + 1:
                x = round(min(max(x, c["b_int"] + 0.01), c["b_int"] + 0.99), 3)
            q = prime_below_real(x)
            if (c["d"], q, c["model"]) in seen:
                continue
            seen.add((c["d"], q, c["model"]))
            jobs.append({"kind": "lwe", "d": c["d"], "L": 1, "b": x, "x": x, "q": q, "model": c["model"]})
    with open(os.path.join(est, "crossings.json"), "w") as f:
        json.dump(crossings(lwe), f, indent=1)
    run(jobs, os.path.join(est, "lwe_fine.jsonl"), a.jobs)


if __name__ == "__main__":
    main()
