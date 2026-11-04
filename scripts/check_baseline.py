#!/usr/bin/env python3
"""Baseline check: recomputes, with the param_report used by the tool (t-evolve-national-scale), every
configuration logged in results/params_chain_signed.jsonl of t-evolve-cryptoalg v1.1 and compares all
bounds and sizes.

The estimator values in the log are not recomputed here (they need SageMath); this script checks
that the part of the chain that precedes the estimator is reproduced exactly.

usage: python3 scripts/check_baseline.py      (after scripts/fetch_deps.sh)
"""
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
LOG = os.path.join(ROOT, "deps", "t-evolve-cryptoalg", "results", "params_chain_signed.jsonl")
BIN = os.path.join(ROOT, "deps", "t-evolve-national-scale", "build", "param_report")
FIELDS = ("E", "sigma_J", "B_J", "sigma1", "sigma2", "B_A", "beta_SIS", "log2beta",
          "ballot_bytes", "auth_per_ballot_bytes")


def run(rec, logq):
    args = [BIN, str(rec["NV"]), str(rec["n"]), str(rec["t"]), str(rec["L"]), str(rec["w"]),
            str(rec["d"]), str(rec["q"]), "1" if rec["mode"] == "signed" else "0", logq]
    return json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)


def main():
    logq, n_ok, n_bad = "64", 0, 0
    with open(LOG) as f:
        for line in f:
            rec = json.loads(line)
            if "# chain" in rec:
                logq = str(rec.get("logQ", "64"))
                continue
            if "beta_SIS" not in rec:
                continue
            new = run(rec, logq)
            bad = [k for k in FIELDS if abs(new[k] - rec[k]) > 1e-6 * max(1.0, abs(rec[k]))]
            tag = "ok " if not bad else "DIFF " + ",".join(bad)
            print(f"NV={rec['NV']:>8} n={rec['n']} t={rec['t']} L={rec['L']} d={rec['d']} "
                  f"log2q={rec['log2q']:.0f} log2beta={new['log2beta']:.3f}  {tag}")
            n_ok += not bad
            n_bad += bool(bad)
    print(f"{n_ok} configurations reproduced, {n_bad} differ")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
