#!/usr/bin/env python3
"""Numbers and table rows of the paper derived from results/derivation.json and param_report.

Prints (and writes to results/paper_numbers.txt): the rows of Tables 2-5, the smallest binding
margin, the effects of each input (Section 6: L, level, cost model, ballot size per modulus bit, step
7 -> 8, size per authority), the authority-side data per ballot, and the difference of beta_SIS between
T-EVOLVE (d = 7, n = 4) and EVOLVE (d = 7, N_A = 4) on the dense electorate grid.

usage: python3 scripts/paper_numbers.py
"""
import json, math, os, sys, re
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from derive import TEVOLVE, TARGETS  # noqa
from estimate_grid import prime_below  # noqa

D = json.load(open(os.path.join(REPO, "results/derivation.json")))["targets"]
NV = [10 ** 4, 10 ** 5, 10 ** 6, 6480000, 10 ** 7]
NVTEX = {10 ** 4: "$10^4$", 10 ** 5: "$10^5$", 10 ** 6: "$10^6$", 6480000: "$6.48\\cdot10^6$", 10 ** 7: "$10^7$"}
out = {}


def cell(t, nv, n=4, L=1):
    return D[t]["cells"]["T-EVOLVE (%d, %d, %d)" % (nv, n, L)]


def hid(c):
    h = c.get("hiding_L2") or c["hiding_at_q"]
    return None if (not h or None in h.values()) else min(h.values())


def kib(c):
    return c["size"] / 1024


# ---- Table 2
rows = []
for nv in NV:
    a, b = cell("128 classical", nv), cell("128 classical+quantum", nv)
    def part(c, prev_d=None):
        dd = ("\\textbf{%d}" % c["d"]) if prev_d is not None and c["d"] != prev_d else "%d" % c["d"]
        return "%s & %.2f (%d) & %.1f / %.1f & %.1f" % (dd, c["log2q"], c["bits"], hid(c), min(c["binding"].values()), kib(c))
    rows.append("%s & %s & %s\\\\" % (NVTEX[nv], part(a, cell("128 classical", NV[max(0, NV.index(nv) - 1)])["d"]),
                                       part(b, cell("128 classical+quantum", NV[max(0, NV.index(nv) - 1)])["d"])))
out["TAB2"] = "\n".join(rows)

# ---- Table 3
rows = []
for nv in NV:
    parts = []
    for t in ("128 core-SVP", "128 core-SVP classical+quantum", "192 classical", "192 classical+quantum"):
        c = cell(t, nv)
        parts.append("%d / %d / %.1f" % (c["d"], c["bits"], kib(c)))
    rows.append("%s & %s\\\\" % (NVTEX[nv], " & ".join(parts)))
out["TAB3"] = "\n".join(rows)


def fmt_nv(x):
    if x is None or x == 0:
        return "--"
    if x >= 10 ** 8:
        return "$\\ge10^8$"
    if x < 1000:
        return "%d" % int(round(float("%.2g" % x)))
    e = int(math.floor(math.log10(x)))
    return "$%.1f\\cdot10^{%d}$" % (x / 10 ** e, e) if e > 1 else "%d" % x


def thr_row(t, label):
    th = {x["d"]: x["N_V_star"] for x in D[t]["thresholds"] if x["L"] == 1}
    cells = []
    for d in range(6, 13):
        if d in th:
            cells.append(fmt_nv(th[d]))
        else:
            below = [dd for dd in th if dd > d]
            cells.append("--" if below else "")
    return label + " & " + " & ".join(cells) + "\\\\"


out["THRROWS"] = "\n".join(thr_row(t, l) for t, l in (
    ("128 classical", "128, MATZOV classical"), ("128 classical+quantum", "128, MATZOV classical + quantum"),
    ("128 core-SVP", "128, core-SVP $0.292\\beta$"), ("128 core-SVP classical+quantum", "128, core-SVP $0.292\\beta$ and $0.265\\beta$"),
    ("192 classical", "192, MATZOV classical"), ("192 classical+quantum", "192, MATZOV classical + quantum")))
out["THR_L2_RATIOS"] = {t: [(x["d"], y["N_V_star"] / x["N_V_star"]) for x in D[t]["thresholds"] for y in D[t]["thresholds"]
                            if x["L"] == 1 and y["L"] == 2 and x["d"] == y["d"] and x["N_V_star"] and 0 < x["N_V_star"] < 10 ** 8]
                        for t in D}
out["CEIL"] = {t: D[t]["ceilings"] for t in D}

# ---- binding minimum over all cells, count cells
mins, ncell = [], 0
for t, td in D.items():
    lam_b = TARGETS[t][2]
    for k, c in td["cells"].items():
        if c.get("d") is None:
            continue
        ncell += 1
        mins.append((min(c["binding"].values()) - lam_b, t, k))
    for k, c in td["cells"].items():
        h = hid(c)
        if h is not None and h < TARGETS[t][0]:
            print("HIDING BELOW TARGET", t, k, h)
        if h is None:
            print("hiding missing", t, k)
out["MINBIND"] = min(mins)
out["NCELLS"] = ncell

# ---- EVOLVE table
rows = []
for na in (2, 3, 4, 5, 7, 9, 15):
    a = D["128 classical"]["cells"]["EVOLVE (%d,)" % na]
    b = D["128 classical+quantum"]["cells"]["EVOLVE (%d,)" % na]
    from evolve_chain import evolve_bounds
    l6 = math.log2(evolve_bounds(6, na)["beta_SIS"])
    assert a["d"] == b["d"] and a["q"] == b["q"]
    rows.append("%d & %.2f & %d & %.2f (%d) & %.1f / %.1f & %.1f & %.2f\\\\" % (
        na, l6, a["d"], a["log2q"], a["bits"], a["hiding_at_q"]["MATZOV"], b["hiding_at_q"]["MATZOV_q"],
        a["binding"]["MATZOV"], kib(a)))
out["TAB7"] = "\n".join(rows)
out["EVOLVE_d"] = {t: sorted({c["d"] for k, c in D[t]["cells"].items() if k.startswith("EVOLVE")}) for t in D}
out["EVOLVE_d_detail"] = {t: [(k, c["d"]) for k, c in D[t]["cells"].items() if k.startswith("EVOLVE")] for t in ("192 classical", "192 classical+quantum", "128 core-SVP", "128 core-SVP classical+quantum")}

# ---- sensitivity numbers
def rng(v):
    return "%.2f-%.2f" % (min(v), max(v))
s = {}
s["L_beta"] = rng([cell(t, nv, 4, 2)["log2beta"] - cell(t, nv)["log2beta"] for t in D for nv in NV])
s["L_bits"] = [(t, nv, cell(t, nv, 4, 2)["bits"] - cell(t, nv)["bits"], cell(t, nv, 4, 2)["d"] - cell(t, nv)["d"]) for t in D for nv in NV
               if cell(t, nv, 4, 2)["bits"] != cell(t, nv)["bits"] or cell(t, nv, 4, 2)["d"] != cell(t, nv)["d"]]
s["L_cmp"] = len(D) * len(NV)
s["L_kib"] = sorted({(cell(t, nv)["d"], round(kib(cell(t, nv, 4, 2)) - kib(cell(t, nv)), 2)) for t in D for nv in NV
                     if cell(t, nv, 4, 2)["d"] == cell(t, nv)["d"]})
s["lvl_beta_cl"] = rng([cell("192 classical", nv)["log2beta"] - cell("128 classical", nv)["log2beta"] for nv in NV])
s["lvl_beta_q"] = rng([cell("192 classical+quantum", nv)["log2beta"] - cell("128 classical+quantum", nv)["log2beta"] for nv in NV])
s["lvl_bits"] = sorted({cell(a, nv)["bits"] - cell(b, nv)["bits"] for a, b in (("192 classical", "128 classical"), ("192 classical+quantum", "128 classical+quantum")) for nv in NV})
s["lvl_rank"] = sorted({cell(a, nv)["d"] - cell(b, nv)["d"] for a, b in (("192 classical", "128 classical"), ("192 classical+quantum", "128 classical+quantum")) for nv in NV})
s["lvl_size_cl"] = rng([100 * (kib(cell("192 classical", nv)) / kib(cell("128 classical", nv)) - 1) for nv in NV])
s["lvl_size_q"] = rng([100 * (kib(cell("192 classical+quantum", nv)) / kib(cell("128 classical+quantum", nv)) - 1) for nv in NV])
s["q_bits"] = sorted({cell("128 classical+quantum", nv)["bits"] - cell("128 classical", nv)["bits"] for nv in NV} |
                     {cell("192 classical+quantum", nv)["bits"] - cell("192 classical", nv)["bits"] for nv in NV})
s["q_size"] = [(t, nv, round(100 * (kib(cell(t + " classical+quantum", nv)) / kib(cell(t + " classical", nv)) - 1), 1)) for t in ("128", "192") for nv in NV
               if cell(t + " classical+quantum", nv)["d"] != cell(t + " classical", nv)["d"]]
s["core_rank"] = sorted({cell("128 core-SVP", nv)["d"] - cell("128 classical", nv)["d"] for nv in NV})
s["core_size"] = rng([100 * (kib(cell("128 core-SVP", nv)) / kib(cell("128 classical", nv)) - 1) for nv in NV])
s["core_bits"] = sorted({cell("128 core-SVP", nv)["bits"] - cell("128 classical", nv)["bits"] for nv in NV})
# ballot per modulus bit at d=7 and step 7->8
sz = {}
for d in (7, 8):
    for b in range(42, 48):
        sz[(d, b)] = TEVOLVE.report((10 ** 4, 4, 1), d, prime_below(b))["ballot_bytes"] / 1024
s["per_bit_d7"] = sorted({round(sz[(7, b + 1)] - sz[(7, b)], 3) for b in range(42, 47)})
s["per_bit_d8"] = sorted({round(sz[(8, b + 1)] - sz[(8, b)], 3) for b in range(42, 47)})
s["step78"] = rng([sz[(8, b)] - sz[(7, b)] for b in range(42, 48)])
s["n_kib"] = ((kib(cell("128 classical", 10 ** 6, 9, 1)) - kib(cell("128 classical", 10 ** 6))) / 5,
              (kib(cell("128 classical+quantum", 10 ** 6, 9, 1)) - kib(cell("128 classical+quantum", 10 ** 6))) / 5)
out["SENS"] = s

# ---- authority side
auth = []
for nv in NV:
    c = cell("128 classical", nv)
    r = TEVOLVE.report((nv, 4, 1), c["d"], prime_below(c["bits"]))
    r2 = TEVOLVE.report((nv, 4, 1), c["d"], prime_below(c["bits"] + 1))
    auth.append((nv, c["d"], c["bits"], r["auth_per_ballot_bytes"] / 1024, r2["auth_per_ballot_bytes"] / 1024,
                 r["auth_round1_bytes"] / 2 ** 20, r["auth_round2_bytes"] / 2 ** 30))
out["AUTH"] = auth
# beta_SIS of T-EVOLVE (d = 7, n = 4, L = 1) minus EVOLVE (d = 7, N_A = 4)
from evolve_chain import evolve_bounds  # noqa
ev = math.log2(evolve_bounds(7, 4)["beta_SIS"])
grid = json.load(open(os.path.join(REPO, "results/sweep_electorate.json")))["grid"]
def at(nv):
    for a, b in zip(grid, grid[1:]):
        if a["NV"] <= nv <= b["NV"]:
            t = (math.log10(nv) - math.log10(a["NV"])) / (math.log10(b["NV"]) - math.log10(a["NV"]))
            return a["log2_beta_SIS"] + t * (b["log2_beta_SIS"] - a["log2_beta_SIS"])
out["TEVOLVE_minus_EVOLVE_d7"] = {"EVOLVE_log2beta": ev, **{str(nv): at(nv) - ev for nv in (100, 2800, 10 ** 4, 10 ** 6, 10 ** 7)}}
# T-EVOLVE cells only: smallest binding margin
mt = min((min(c["binding"].values()) - TARGETS[t][2], t, k) for t, td in D.items() for k, c in td["cells"].items()
         if k.startswith("T-EVOLVE") and c.get("d") is not None)
out["MINBIND_TEVOLVE"] = mt
lines = []
for k, v in out.items():
    lines.append("===== " + k)
    lines.append(v if isinstance(v, str) else json.dumps(v, indent=0, default=str))
print("\n".join(lines))
with open(os.path.join(REPO, "results", "paper_numbers.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
