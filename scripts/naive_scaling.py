#!/usr/bin/env python3
"""Simpler parameter strategies against the derived set at 6.48 million voters (paragraph "Naive scaling" of the paper).

128 bits, classical MATZOV, n = 4, L = 1. The strategies:
  reuse10^6   the set derived for 10^6 voters (d = 7, its modulus) used at 6.48e6;
  samerank    d = 7 with the smallest admissible modulus above beta_SIS(6.48e6, d = 7);
  linear      beta_SIS extrapolated linearly in N_V from 10^6, smallest modulus above it, d = 7;
  national    the derived 6.48e6 set used at 10^7 voters; also the largest electorate it serves, with its
              own modulus and with the largest prime q = 17 (mod 32) of the same bit length;
  oneforall   the derived 6.48e6 set used for 10^4 voters (size comparison).
The M-LWE estimates at the new moduli are computed with the estimator (needs Sage) and stored in
results/estimates/lwe_naive.jsonl; the rest uses param_report. Output: results/naive.{json,md}.

usage: python3 scripts/naive_scaling.py [--jobs J]
"""
import argparse
import json
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from estimate_grid import prime_below, run  # noqa: E402
from derive import TEVOLVE, prime_above, EST, Hardness  # noqa: E402

NAT = 6480000


def beta(nv, d):
    return TEVOLVE.beta((nv, 4, 1), d)


def largest_nv(d, q, lo=10 ** 6, hi=10 ** 8):
    """Largest N_V with beta_SIS(N_V, d) < q (bisection, 0.02 % resolution)."""
    while hi / lo > 1.0002:
        mid = int(round(math.sqrt(lo * hi)))
        if beta(mid, d) < q:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    b6, b_nat7 = beta(10 ** 6, 7), beta(NAT, 7)
    q6 = prime_above(b6)                                   # modulus of the derived 10^6 set
    q_same = prime_above(b_nat7)
    q_lin = prime_above(b6 * NAT / 10 ** 6)
    jobs = [{"kind": "lwe", "d": 7, "L": 1, "b": round(math.log2(q), 6), "q": q, "model": "MATZOV"} for q in (q_same, q_lin)]
    run(jobs, os.path.join(EST, "lwe_naive.jsonl"), a.jobs)
    H = Hardness()
    with open(os.path.join(EST, "lwe_naive.jsonl")) as f:
        hid = {json.loads(x)["q"]: json.loads(x)["min"] for x in f}
    r7 = TEVOLVE.report((NAT, 4, 1), 7)
    bud7 = TEVOLVE.budget((NAT, 4, 1), 7)
    b_nat8 = beta(NAT, 8)
    q_nat = prime_above(b_nat8)
    bits_nat = q_nat.bit_length()
    q_top = prime_below(bits_nat)
    out = {
        "reuse_1e6": {"d": 7, "log2q": math.log2(q6), "log2beta": math.log2(b_nat7), "log2sigma1": math.log2(r7["sigma1"]),
                      "two_level_max_log2": 27.5, "log2_budget_tail_at_187": math.log2(30 * (r7["E"] + 1) + 1) - 187,
                      "E": r7["E"], "log2_budget_derived": bud7},
        "same_rank": {"d": 7, "log2q": math.log2(q_same), "bits": q_same.bit_length(), "hiding_MATZOV": hid[q_same]},
        "linear": {"d": 7, "log2q": math.log2(q_lin), "bits": q_lin.bit_length(), "hiding_MATZOV": hid[q_lin]},
        "national_at_1e7": {"d": 8, "log2q": math.log2(q_nat), "log2beta_1e7": math.log2(beta(10 ** 7, 8)),
                            "max_NV_own_q": largest_nv(8, q_nat, lo=NAT), "max_NV_top_q": largest_nv(8, q_top, lo=NAT),
                            "log2q_top": math.log2(q_top), "hiding_top_MATZOV": H.hiding(8, q_top, "MATZOV")},
        "one_for_all": {"size_nat_KiB": TEVOLVE.size((10 ** 4, 4, 1), 8, bits_nat) / 1024,
                        "size_derived_1e4_KiB": TEVOLVE.size((10 ** 4, 4, 1), 7, prime_above(beta(10 ** 4, 7)).bit_length()) / 1024},
    }
    with open(os.path.join(ROOT, "results", "naive.json"), "w") as f:
        json.dump(out, f, indent=1)
    md = ["# Simpler strategies at 6.48e6 voters (128 bits, classical MATZOV, n = 4, L = 1)", "",
          "`scripts/naive_scaling.py`", "", "```", json.dumps(out, indent=1), "```", ""]
    with open(os.path.join(ROOT, "results", "naive.md"), "w") as f:
        f.write("\n".join(md))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
