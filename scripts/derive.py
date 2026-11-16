#!/usr/bin/env python3
"""Parameter derivation (step S6 of the method), generic over the protocol.

A protocol module supplies, for a configuration and a module rank d,
  * beta(cfg, d)        the M-SIS bound of steps S1-S4 (independent of q, Lemma 2),
  * feasible(cfg, d)    the feasibility conditions S5(i), (ii), (iv) that do not involve q,
  * size(cfg, d, bits)  the size that the selection minimises (bits = bit length of q),
  * L(cfg)              the number of message slots (shape of the M-LWE / M-SIS instances).
Two modules are provided: TEVOLVE (bounds from param_report of the T-EVOLVE code) and EVOLVE (the
constraints of EVOLVE's Section 5, binding at 2 B_r). Everything else is generic:

  hiding     M-LWE estimates from results/estimates/lwe.jsonl and lwe_fine.jsonl (integer and
             sub-bit moduli); by (H4) the hiding ceiling q_hi(d) is the largest computed q whose
             estimate reaches the target in every model of the set (conservative: the true ceiling
             lies between q_hi and the smallest computed q that fails);
  modulus    any prime q = 17 (mod 32); the selection takes q_lo(d), the smallest such prime above
             beta, and the size depends on the bit length of q only;
  hiding at q  M-LWE at the selected q (results/estimates/lwe_sel.jsonl, L = 1; lwe_L2.jsonl, L = 2);
  binding    M-SIS estimates at (d, q_lo, beta) from results/estimates/sis_fine.jsonl (computed by
             `derive.py estimate`), falling back to the integer tables when q_lo is a table modulus;
  selection  the smallest d with q_lo(d) <= q_hi(d) and binding >= target; Lemma 7 of the paper
             shows that it also minimises the size;
  thresholds for TEVOLVE, the largest N_V served by each rank (bisection with param_report on the
             range 10^2 <= N_V <= 10^8, to a relative precision of 0.02 %; "none" means that the rank
             is not admissible anywhere in that range).

Targets: separate levels and model sets for hiding and binding.

usage: python3 scripts/derive.py estimate [--jobs J]   # needs Sage: M-SIS at the selected moduli, M-LWE with L = 2
       python3 scripts/derive.py select                # results/derivation.{json,md}
"""
import argparse
import json
import math
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from estimate_grid import is_prime, prime_below, run, BIN, N  # noqa: E402
from evolve_chain import evolve_bounds  # noqa: E402

EST = os.path.join(ROOT, "results", "estimates")
D_RANGE = range(5, 13)
Q_MAX_BITS = 50

# targets: name -> (hiding level, hiding models, binding level, binding models)
TARGETS = {
    "128 classical": (128, ("MATZOV",), 128, ("MATZOV",)),
    "128 classical+quantum": (128, ("MATZOV", "MATZOV_q"), 128, ("MATZOV", "MATZOV_q")),
    "128 core-SVP": (128, ("ADPS16",), 128, ("ADPS16",)),
    "128 core-SVP classical+quantum": (128, ("ADPS16", "ADPS16_q"), 128, ("ADPS16", "ADPS16_q")),
    "192 classical": (192, ("MATZOV",), 192, ("MATZOV",)),
    "192 classical+quantum": (192, ("MATZOV", "MATZOV_q"), 192, ("MATZOV", "MATZOV_q")),
    "hiding 192 / binding 128, classical": (192, ("MATZOV",), 128, ("MATZOV",)),
}
MODELS = ("MATZOV", "MATZOV_q", "ADPS16", "ADPS16_q")


def prime_above(beta):
    """Smallest prime q = 17 (mod 32) with q > beta."""
    q = int(math.floor(beta)) + 1
    q += (17 - q) % 32
    while not is_prime(q):
        q += 32
    return q


# ------------------------------------------------------------------ protocol modules

CACHE = os.path.join(ROOT, "deps", "param_report_cache.jsonl")   # outputs of param_report (not tracked)


class TEVOLVE:
    name = "T-EVOLVE"
    CONFIGS = [(nv, 4, L) for nv in (10 ** 4, 10 ** 5, 10 ** 6, 6480000, 10 ** 7) for L in (1, 2)] + [(10 ** 6, 9, 1)]
    _cache = None

    @classmethod
    def report(cls, cfg, d, q=None):
        """Output of param_report (deterministic; cached on disk in deps/)."""
        if cls._cache is None:
            cls._cache = {}
            if os.path.exists(CACHE):
                with open(CACHE) as f:
                    for line in f:
                        r = json.loads(line)
                        cls._cache[tuple(r["key"])] = r["out"]
        nv, n, L = cfg
        q = q or prime_below(45)
        k = (nv, n, L, d, q)
        if k not in cls._cache:
            out = subprocess.run([BIN, str(nv), str(n), "1", str(L), "-1", str(d), str(q), "1", "64"],
                                 capture_output=True, text=True, check=True).stdout
            cls._cache[k] = json.loads(out)
            os.makedirs(os.path.dirname(CACHE), exist_ok=True)
            with open(CACHE, "a") as f:
                f.write(json.dumps({"key": list(k), "out": cls._cache[k]}) + "\n")
        return cls._cache[k]

    @classmethod
    def L(cls, cfg):
        return cfg[2]

    @classmethod
    def beta(cls, cfg, d):
        return cls.report(cfg, d)["beta_SIS"]

    @classmethod
    def budget(cls, cfg, d):
        """Per-attempt statistical error of the aggregation proof (log2), S5(iv)."""
        r = cls.report(cfg, d)
        E, l = r["E"], 30
        mu = 2 * d + cfg[2]
        tau_bits = max(187.0, math.ceil(171.0 + math.log2(l * (E + 1) + 1)))
        tail = math.log2(l * (E + 1) + 1) - tau_bits
        rr, lnM = 309 / 20, 86929 / 157922
        rej = math.log2(2 * math.exp(-rr * rr / 2) / math.exp(lnM))
        levels3 = max(r["sigma1"], r["sigma2"]) > 2 ** 27.5
        samp = math.log2(E * N * mu * (500 + 517)) + (-215 if levels3 else -216)
        fresh = -192
        return math.log2(sum(2 ** x for x in (tail, rej, samp, fresh)))

    @classmethod
    def feasible(cls, cfg, d):
        r = cls.report(cfg, d)
        if min(r["sigma_J"], r["sigma1"], r["sigma2"]) <= 0:
            return False, "no sampler plan"
        if cls.budget(cfg, d) > -170.3:
            return False, "statistical budget"
        return True, ""

    @classmethod
    def size(cls, cfg, d, bits):
        return cls.report(cfg, d, prime_below(bits))["ballot_bytes"]


class EVOLVE:
    name = "EVOLVE"
    CONFIGS = [(na,) for na in (2, 3, 4, 5, 7, 9, 15)]

    @classmethod
    def L(cls, cfg):
        return 1

    @classmethod
    def beta(cls, cfg, d):
        return evolve_bounds(d, cfg[0])["beta_SIS"]

    @classmethod
    def feasible(cls, cfg, d):
        return True, ""

    @classmethod
    def size(cls, cfg, d, bits):
        return (d + 1) * N * bits / 8          # one commitment


# ------------------------------------------------------------------ hardness data

def load(name):
    path = os.path.join(EST, name)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [json.loads(line) for line in f]


class Hardness:
    def __init__(self):
        self.lwe = {}                            # (d, L, model) -> list of (q, min)
        for r in load("lwe.jsonl") + load("lwe_fine.jsonl") + load("lwe_L2.jsonl") + load("lwe_sel.jsonl"):
            self.lwe.setdefault((r["d"], r["L"], r["model"]), []).append((r["q"], r["min"]))
        for v in self.lwe.values():
            v.sort()
        self.sis = {}                            # (d, L, q, beta-rounded, model) -> min
        for r in load("sis.jsonl") + load("sis_fine.jsonl") + \
                [dict(x, L=1) for x in load(os.path.join("..", "estimates_evolve", "sis.jsonl"))]:
            self.sis[(r["d"], r["L"], r["q"], round(math.log2(r["beta_SIS"]), 9), r["model"])] = r["min"]

    def ceiling(self, d, lam, models, L=1):
        """Largest computed q with M-LWE >= lam in every model (None if none in the data)."""
        best = None
        for m in models:
            pts = self.lwe.get((d, 1, m), [])
            ok = [q for q, v in pts if v >= lam]
            if not ok:
                return None
            best = max(ok) if best is None else min(best, max(ok))
        return best

    def hiding(self, d, q, m, L=1):
        for qq, v in self.lwe.get((d, L, m), []):
            if qq == q:
                return v
        return None

    def binding(self, d, L, q, beta, m):
        return self.sis.get((d, L, q, round(math.log2(beta), 9), m))


# ------------------------------------------------------------------ selection

def select_cell(P, cfg, H, target, need):
    lam_h, mh, lam_b, mb = target
    L = P.L(cfg)
    for d in D_RANGE:
        ok, why = P.feasible(cfg, d)
        if not ok:
            continue
        beta = P.beta(cfg, d)
        q_lo = prime_above(beta)
        q_hi = H.ceiling(d, lam_h, mh)
        if q_hi is None or q_lo > q_hi or q_lo >= 2 ** Q_MAX_BITS:
            continue
        bind = {m: H.binding(d, L, q_lo, beta, m) for m in mb}
        if any(v is None for v in bind.values()):
            need.update(("sis", P.name, cfg, d, q_lo, beta, m) for m, v in bind.items() if v is None)
            return {"d": d, "q": q_lo, "pending": True}
        if min(bind.values()) < lam_b:
            continue                              # binding would need a larger q: not seen on the grid
        bits = q_lo.bit_length()
        rec = {"d": d, "q": q_lo, "log2q": math.log2(q_lo), "bits": bits, "log2beta": math.log2(beta),
               "q_hi_log2": math.log2(q_hi), "binding": bind, "size": P.size(cfg, d, bits)}
        if L == 2:
            h2 = {m: H.hiding(d, q_lo, m, L=2) for m in mh}
            if any(v is None for v in h2.values()):
                need.update(("lwe2", d, q_lo, m) for m, v in h2.items() if v is None)
            elif min(h2.values()) < lam_h:
                continue                          # the L = 2 instance misses the target at q_lo: next rank
            rec["hiding_L2"] = h2
        h1 = {m: H.hiding(d, q_lo, m) for m in mh}
        if L == 1 and any(v is None for v in h1.values()):
            need.update(("lwe1", d, q_lo, m) for m, v in h1.items() if v is None)
        rec["hiding_at_q"] = h1
        if hasattr(P, "budget"):
            rec["log2_budget"] = P.budget(cfg, d)
        return rec
    return {"d": None}


def threshold(L, d, q_hi, lo=100, hi=10 ** 8):
    def beta(nv):
        return TEVOLVE.beta((nv, 4, L), d)
    if beta(lo) >= q_hi:
        return 0
    if beta(hi) < q_hi:
        return hi
    while hi - lo > 1 and hi / lo > 1.0002:
        mid = min(hi - 1, max(lo + 1, int(round(math.sqrt(lo * hi)))))
        if beta(mid) < q_hi:
            lo = mid
        else:
            hi = mid
    return lo


def dual_free_ceiling(d, lam, models):
    """Ceiling when the dual and dual-hybrid attacks are excluded (from the stored attack lists)."""
    best = None
    for m in models:
        pts = []
        for r in load("lwe.jsonl") + load("lwe_fine.jsonl") + load("lwe_sel.jsonl") + load("lwe_naive.jsonl"):
            if r["d"] == d and r["model"] == m and r["L"] == 1:
                v = [x for k, x in r["attacks"].items() if x is not None and not k.startswith("dual")]
                pts.append((r["q"], min(v)))
        ok = [q for q, v in pts if v >= lam]
        if not ok:
            return None
        best = max(ok) if best is None else min(best, max(ok))
    return best


def selection(with_thresholds=True):
    H = Hardness()
    need = set()
    out = {"targets": {}}
    for tname, target in TARGETS.items():
        t = {"target": list(target), "ceilings": {}, "cells": {}, "thresholds": []}
        for d in D_RANGE:
            c = H.ceiling(d, target[0], target[1])
            t["ceilings"][d] = None if c is None else math.log2(c)
        for P in (TEVOLVE, EVOLVE):
            for cfg in P.CONFIGS:
                t["cells"]["%s %s" % (P.name, cfg)] = select_cell(P, cfg, H, target, need)
        for L in ((1, 2) if with_thresholds else ()):
            for d in D_RANGE:
                c = H.ceiling(d, target[0], target[1])
                if c is None:
                    continue
                ns = threshold(L, d, c)
                t["thresholds"].append({"L": L, "d": d, "q_hi_log2": math.log2(c), "N_V_star": ns})
                if ns == 10 ** 8:
                    break
        out["targets"][tname] = t
        print("target", tname, "done", flush=True)
    dual = {}
    for lam in (128, 192):
        for d in D_RANGE:
            a, b = H.ceiling(d, lam, ("MATZOV",)), dual_free_ceiling(d, lam, ("MATZOV",))
            if a and b:
                dual["%d %d" % (lam, d)] = [math.log2(a), math.log2(b)]
    out["dual_free_ceilings_MATZOV"] = dual
    return out, need


def write_md(out):
    md = ["# Parameter derivation", "",
          "`scripts/derive.py select`. Modulus: smallest prime q = 17 (mod 32) above beta_SIS; size by bit length.",
          "T-EVOLVE: ballot bytes (n, L as given); EVOLVE: one commitment. q_hi: hiding ceiling (largest computed q).", ""]
    for tname, t in out["targets"].items():
        md += ["## " + tname, "", "Hiding ceilings log2 q_hi(d): " + ", ".join(
            "%s: %s" % (d, "-" if v is None else "%.2f" % v) for d, v in t["ceilings"].items()), "",
            "| cell | d | log2 q | bits | log2 beta | hiding at q (min over models) | binding (min over models) | size (KiB) |",
            "|---|---|---|---|---|---|---|---|"]
        for k, r in t["cells"].items():
            if r.get("d") is None:
                md.append("| %s | none | | | | | | |" % k)
            elif r.get("pending"):
                md.append("| %s | %d (binding pending) | | | | | | |" % (k, r["d"]))
            else:
                flag = ""
                h2 = r.get("hiding_L2") or {}
                if h2 and (None in h2.values() or min(h2.values()) < t["target"][0]):
                    flag = " (L=2 hiding: %s)" % h2
                hv = r.get("hiding_L2") if r.get("hiding_L2") else r.get("hiding_at_q")
                hs = "-" if not hv or None in hv.values() else "%.1f" % min(hv.values())
                md.append("| %s | %d | %.2f | %d | %.2f | %s | %.1f | %.2f%s |" % (k, r["d"], r["log2q"], r["bits"], r["log2beta"],
                          hs, min(r["binding"].values()), r["size"] / 1024, flag))
        md += ["", "Thresholds (largest N_V per rank): " + "; ".join(
            "L=%d d=%d: %s" % (x["L"], x["d"], ">=1e8" if x["N_V_star"] == 10 ** 8 else ("none" if x["N_V_star"] == 0 else "%.3g" % x["N_V_star"]))
            for x in t["thresholds"]), ""]
    md += ["## Ceilings with and without the dual attacks (MATZOV)", "",
           "; ".join("lambda %s, d %s: %.2f / %.2f" % (tuple(k.split()) + tuple(v)) for k, v in out["dual_free_ceilings_MATZOV"].items()), ""]
    with open(os.path.join(ROOT, "results", "derivation.md"), "w") as f:
        f.write("\n".join(md) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=("estimate", "select"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    out, need = selection(with_thresholds=(a.what == "select"))
    if a.what == "estimate":
        sis = [{"kind": "sis", "scheme": s, "NV": (cfg[0] if s == "T-EVOLVE" else 0), "n": (cfg[1] if s == "T-EVOLVE" else cfg[0]),
                "L": (cfg[2] if s == "T-EVOLVE" else 1), "d": d, "b": round(math.log2(q), 6), "q": q, "beta_SIS": beta, "model": m}
               for (kind, s, cfg, d, q, beta, m) in sorted(x for x in need if x[0] == "sis")]
        run(sis, os.path.join(EST, "sis_fine.jsonl"), a.jobs)
        out, need = selection(with_thresholds=False)
        lwe2 = [{"kind": "lwe", "d": d, "L": 2, "b": round(math.log2(q), 6), "q": q, "model": m}
                for (kind, d, q, m) in sorted(x for x in need if x[0] == "lwe2")]
        run(lwe2, os.path.join(EST, "lwe_L2.jsonl"), a.jobs)
        lwe1 = [{"kind": "lwe", "d": d, "L": 1, "b": round(math.log2(q), 6), "q": q, "model": m}
                for (kind, d, q, m) in sorted(x for x in need if x[0] == "lwe1")]
        run(lwe1, os.path.join(EST, "lwe_sel.jsonl"), a.jobs)
        return
    pend = [k for t in out["targets"].values() for k, r in t["cells"].items() if r.get("pending")]
    if pend:
        print("binding estimates missing for %d cells: run `derive.py estimate` first" % len(pend))
    with open(os.path.join(ROOT, "results", "derivation.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    write_md(out)


if __name__ == "__main__":
    main()
