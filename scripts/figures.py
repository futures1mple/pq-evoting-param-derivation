#!/usr/bin/env python3
"""Figures of the paper from the results of sweep_electorate.py and derive.py.

  fig_scaling: log2 beta_SIS against the number of voters for module ranks 7 and 8, with the hiding
               ceilings of rank 7 (largest log2 q at which M-LWE still reaches 2^128, at 0.1-bit
               resolution): the rank must increase where the curve crosses the ceiling.
  fig_rank:    selected module rank against the number of voters for four targets.

The rank-8 curve is computed with param_report on first use (about 10 minutes) and stored in
results/sweep_d8.json. Output: results/figures/*.pdf and *.png.

usage: python3 scripts/figures.py
"""
import json
import math
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from estimate_grid import prime_below, BIN  # noqa: E402

RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "figures")
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3 (validated all-pairs)

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 2,
                     "legend.frameon": False, "pdf.fonttype": 42})


def sweep_d8():
    path = os.path.join(RES, "sweep_d8.json")
    if not os.path.exists(path):
        out = []
        for i in range(56):
            nv = int(round(10 ** (2 + 0.1 * i)))
            r = json.loads(subprocess.run([BIN, str(nv), "4", "1", "1", "-1", "8", str(prime_below(45)), "1", "64"],
                                          capture_output=True, text=True, check=True).stdout)
            out.append({"NV": nv, "log2_beta_SIS": r["log2beta"]})
        with open(path, "w") as f:
            json.dump(out, f)
    with open(path) as f:
        return json.load(f)


def target(sel, name):
    return sel["targets"][name]


def rank_steps(t, L=1):
    """(N_V*, d) pairs: rank d serves up to N_V*(d); the selected rank is the smallest such d."""
    thr = sorted((x["d"], x["N_V_star"]) for x in t["thresholds"] if x["L"] == L and (x["N_V_star"] or 0) > 0)
    steps, cur = [], 0
    for d, nstar in thr:
        if nstar > cur:
            steps.append((cur, nstar, d))
            cur = nstar
    return steps


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, name + "." + ext), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_scaling(sel):
    with open(os.path.join(RES, "sweep_electorate.json")) as f:
        g7 = json.load(f)["grid"]
    g8 = sweep_d8()
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.plot([r["NV"] for r in g7], [r["log2_beta_SIS"] for r in g7], color=C1, label=r"rank $d = 7$")
    ax.plot([r["NV"] for r in g8], [r["log2_beta_SIS"] for r in g8], color=C2, label=r"rank $d = 8$")
    for name, text in (("128 classical", r"$d = 7$ limit, classical"),
                       ("128 classical+quantum", r"$d = 7$ limit, classical and quantum")):
        t = target(sel, name)
        y = t["ceilings"]["7"]
        ax.axhline(y, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.text(1.3e2, y + 0.15, text + r": $q \leq 2^{%.2f}$" % y, color=INK2, fontsize=7.5)
        nstar = next(x["N_V_star"] for x in t["thresholds"] if x["L"] == 1 and x["d"] == 7)
        ax.plot([nstar], [y], marker="o", markersize=5, color=INK, zorder=5)
        ax.annotate(r"$%.1f\cdot 10^{%d}$ voters" % (nstar / 10 ** int(math.log10(nstar)), int(math.log10(nstar))), (nstar, y), textcoords="offset points", xytext=(6, -12), fontsize=7.5, color=INK)
    ax.axvline(6.48e6, color=GRID, linewidth=1)
    ax.text(6.48e6, 47.0, r"$6.48\cdot 10^6$ ", color=INK2, fontsize=7.5, va="top", ha="right")
    ax.set_xscale("log")
    ax.set_xlim(1e2, 3.2e7)
    ax.set_xlabel(r"number of voters $N_V$")
    ax.set_ylabel(r"$\log_2\beta_{\mathrm{SIS}}$ (bits)")
    ax.legend(loc="lower right")
    save(fig, "fig_scaling")


def fig_rank(sel):
    fig, ax = plt.subplots(figsize=(5.2, 2.8))
    for (name, label, c, ls) in (("128 classical", "128 bit, classical", C1, "-"),
                                 ("128 classical+quantum", "128 bit, classical and quantum", C2, (0, (4, 2))),
                                 ("128 core-SVP", "128 bit, core-SVP", INK2, (0, (1, 1.5))),
                                 ("192 classical", "192 bit (hardness), classical", C3, "-")):
        xs, ys = [], []
        for lo, hi, d in rank_steps(target(sel, name)):
            xs += [max(lo, 1e2), min(hi, 1e8)]
            ys += [d, d]
        ax.plot(xs, ys, color=c, label=label, linestyle=ls)
    ax.set_xscale("log")
    ax.set_xlim(1e2, 1e8)
    ax.set_xlabel(r"number of voters $N_V$")
    ax.set_ylabel(r"module rank $d$")
    ax.set_yticks(range(6, 13))
    ax.set_ylim(5.6, 13.9)
    ax.legend(loc="upper left", ncol=2, fontsize=7.5, handlelength=2.6, columnspacing=1.2)
    save(fig, "fig_rank")


def main():
    with open(os.path.join(RES, "derivation.json")) as f:
        sel = json.load(f)
    fig_scaling(sel)
    fig_rank(sel)
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
