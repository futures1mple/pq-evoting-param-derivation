#!/usr/bin/env python3
"""Structural checks of the bound chain (no estimator needed).

For T-EVOLVE (n = 4, t = 3, L = 1, free weight, signed second challenge) this script
  1. sweeps the number of voters N_V over a logarithmic grid and checks that E, sigma_1, sigma_2,
     B_A, beta_BL and beta_SIS never decrease (Lemma 3 of the method);
  2. reports the dominant term of beta_SIS and the slope of log2 beta_SIS per decade
     (square-root law: 0.5 log2 10 = 1.661);
  3. checks that the bounds do not depend on q or on t, and are nondecreasing in d, n and L;
  4. finds, by bisection, the largest N_V for which beta_SIS stays below the largest prime
     q = 17 (mod 32) under 2^b, and the N_V at which sigma_1 needs a third sampler level.

Writes results/sweep_electorate.json and prints a summary. Run time about 10-15 minutes, dominated by
calls above about 3 * 10^6 voters (third sampler level).

usage: python3 scripts/sweep_electorate.py      (after scripts/fetch_deps.sh)
"""
import json
import math
import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BIN = os.path.join(ROOT, "deps", "t-evolve-national-scale", "build", "param_report")
OUT = os.path.join(ROOT, "results", "sweep_electorate.json")
MONO = ("E", "sigma1", "sigma2", "B_A", "beta_BL", "beta_SIS")
TERMS = {"root": "t_root", "tree": "t_tree", "ballot": "t_ballot", "or": "t_or", "BL": "beta_BL"}


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


Q45 = prime_below(45)


def rep(NV, n=4, t=3, L=1, d=7, q=Q45):
    args = [BIN, str(NV), str(n), str(t), str(L), "-1", str(d), str(q), "1", "64"]
    return json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)


def dominant(r):
    return max(TERMS, key=lambda k: r[TERMS[k]])


def largest_below(pred, lo, hi):
    """Largest N_V in [lo, hi] with pred false, for pred monotone (false then true)."""
    if not pred(hi):
        return hi
    while hi - lo > max(1, lo // 2000):
        mid = (lo + hi) // 2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return lo


def main():
    out = {"param_report": "t-evolve-national-scale d0885a3", "config": "n=4 t=3 L=1 w=free signed"}
    grid = sorted(set([int(round(10 ** (2 + i * 0.02))) for i in range(236)] +
                      [int(10 ** x) for x in (6.75, 6.8, 6.9, 7.0, 7.25, 7.5)] + [6480000]))
    rows, viol = {}, []
    prev = None
    for NV in grid:
        r = rep(NV)
        rows[NV] = r
        if prev is not None:
            viol += [(k, prev, NV) for k in MONO if r[k] < rows[prev][k] * (1 - 1e-12)]
        prev = NV
    out["grid"] = [{"NV": NV, "E": r["E"], "log2_sigma1": math.log2(r["sigma1"]),
                    "log2_sigma2": math.log2(r["sigma2"]), "log2_beta_SIS": r["log2beta"],
                    "dominant": dominant(r)} for NV, r in rows.items()]
    out["monotonicity_violations"] = viol
    print("N_V grid: %d points from 1e2 to 10^7.5; monotonicity violations: %d" % (len(grid), len(viol)))
    for NV in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 6480000, 10 ** 7, int(10 ** 7.5)):
        r = rows[NV]
        print("  N_V=%-9d E=%-6d log2 sigma1=%.2f log2 beta_SIS=%.3f dominant=%s"
              % (NV, r["E"], math.log2(r["sigma1"]), r["log2beta"], dominant(r)))
    slope = (rows[int(10 ** 7.5)]["log2beta"] - rows[10 ** 6]["log2beta"]) / 1.5
    out["slope_per_decade_1e6_to_1e7.5"] = slope
    print("slope of log2 beta_SIS per decade, 1e6..10^7.5: %.3f (square-root law 1.661)" % slope)

    qs = [prime_below(b) for b in (30, 38, 45, 49)]
    out["q_independent"] = len({rep(10 ** 5, q=q)["beta_SIS"] for q in qs}) == 1
    out["t_independent"] = len({rep(10 ** 5, n=9, t=t)["beta_SIS"] for t in range(1, 10)}) == 1
    for name, rng, kw in (("d", range(1, 13), "d"), ("n", range(2, 16), "n"), ("L", range(1, 11), "L")):
        vals = [rep(10 ** 5, **({kw: x, "t": 2} if kw == "n" else {kw: x}))["beta_SIS"] for x in rng]
        out["monotone_in_" + name] = all(a <= b for a, b in zip(vals, vals[1:]))
        out["log2_beta_range_" + name] = [math.log2(vals[0]), math.log2(vals[-1])]
    print("independent of q: %s; independent of t: %s; nondecreasing in d, n, L: %s, %s, %s"
          % (out["q_independent"], out["t_independent"], out["monotone_in_d"], out["monotone_in_n"],
             out["monotone_in_L"]))

    cross = {}
    for d in (7, 8):
        for b in (45, 46, 47):
            q = prime_below(b)
            cross["d%d_q%d" % (d, b)] = largest_below(lambda x: rep(x, d=d, q=Q45)["beta_SIS"] >= q,
                                                       10 ** 4, 10 ** 8)
            print("d=%d: beta_SIS < q_%d up to N_V ~ %.3g" % (d, b, cross["d%d_q%d" % (d, b)]))
    cross["sigma1_above_2^27.5_d7"] = largest_below(lambda x: math.log2(rep(x)["sigma1"]) > 27.5, 10 ** 6, 10 ** 7)
    print("sigma_1 above 2^27.5 (two sampler levels) from N_V ~ %.3g (d=7)" % cross["sigma1_above_2^27.5_d7"])
    out["crossings"] = cross
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
