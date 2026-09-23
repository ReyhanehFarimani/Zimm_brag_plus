#!/usr/bin/env python3
"""Transition figure of the TABLE J scan (runs/<geom>/jscan_db, jscan_db_gen.sh): helicity theta(J),
chi_theta = N var(theta), handedness <|m|>(J), Binder U4(J) in the top row; <Rg^2>, <E_nb>, <E_twist>
in the bottom row.  Blue curves: exact 1D chain without non-bonded interactions for the geometry's
parameter file (../sample_data.dat, transfer matrix of ../../t100/exact_fss.py), N = 100 .. 1600;
red points: MC N = 200 with the six-argument helix-helix table + registry twist term; filled = finished
runs, hollow = still running (user 2026-09-23: "plot even if they are not finished"; every production
row written so far, at least 100 of them), errors = max(SEM over seeds, combined 20-block error).  Dotted line: J*(N = 200) of the exact
chain (4.68 at theta0 = 45, 4.19 at 100).  Style: no grid, ticks inside, no title.

  cd runs/t45/jscan_db && python ../../jscan_db_plot.py      -> transitions_db.png / .pdf + numbers
"""
import contextlib
import glob
import importlib.util
import io
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.getcwd()
GEOM = os.path.basename(os.path.dirname(HERE))
N0 = 200; NS = (100, 200, 400, 800, 1600); Jg = np.arange(0.5, 10.501, 0.125)
JSTAR = {"t45": 4.68, "t100": 4.19}
# registry free energy of a same-handed junction relative to a free registry (the entropy the explicit
# registry loses when locked in the harmonic well): F = -ln[ sqrt(2 pi / kappa) erf(pi sqrt(kappa/2)) / 2 pi ]
from math import erf, sqrt, pi, log
KAPPA = {"t45": 6.6014, "t100": 1.8177}                      # R.R / L.L junction (kT/rad^2)
KAPPA_RL = {"t45": 2 * 2.465, "t100": 2 * 0.108}              # R.L junction
def f_reg(k): return -log(sqrt(2 * pi / k) * erf(pi * sqrt(k / 2)) / (2 * pi))
F_REG = {g: f_reg(k) for g, k in KAPPA.items()}
F_RL = {g: f_reg(k) for g, k in KAPPA_RL.items()}
# THEORY WITH REGISTRIES (user 2026-09-23, "the theory should be updated"): for harmonic junction
# registries the run partition function factorises, every same-handed junction contributes
# exp(-F_HH) and every R|L junction exp(-F_RL) relative to free registries (coil-helix junctions and
# chain ends are free), so the exact transfer matrix applies with E_HH = -J + F_HH, E_RL = +J + F_RL.
CN = dict(zip(NS, ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]))
RED, INK, MUTED = "#c0392f", "#0b0b0b", "#898781"

# ---------------------------------------------------------------- exact 1D reference
sys.argv = ["x", "../sample_data.dat"]
spec = importlib.util.spec_from_file_location("ex", os.path.abspath("../../t100/exact_fss.py")); ex = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(ex)
PAIRS = sorted(ex.idx, key=ex.idx.get); SECOND = np.array([b for (a, b) in PAIRS])


def exact_dist(J, N, helix=False, registries=False):
    if registries:
        ex.J0, ex.J1, ex.J2 = float(J) - F_REG[GEOM], 0.0, float(J) + F_RL[GEOM]
    else:
        ex.J0, ex.J1, ex.J2 = float(J), 0.0, float(J)
    T, v0 = ex.build(0.0)
    f = (lambda s: abs(s)) if helix else (lambda s: s); off = 0 if helix else N
    P = np.zeros((9, (N + 1) if helix else (2 * N + 1)))
    for p, (a, b) in enumerate(PAIRS):
        P[p, off + f(a) + f(b)] = v0[p]
    for _ in range(N - 2):
        Q = T.T @ P
        for q in range(9):
            Q[q] = np.roll(Q[q], f(SECOND[q]))
        P = Q / Q.sum()
    p = P.sum(0); return p / p.sum()


def exact_table(registries):
    E = {N: {k: np.zeros(len(Jg)) for k in ("th", "chi_th", "am", "U4")} for N in NS}
    cache = "exact_cache_db_reg.npz" if registries else "exact_cache_db.npz"
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime("../sample_data.dat") and len(np.load(cache)[f"{NS[0]}_th"]) == len(Jg):
        z = np.load(cache)
        for N in NS:
            for k in E[N]:
                E[N][k] = z[f"{N}_{k}"]
        return E
    for N in NS:
        m = np.arange(-N, N + 1) / N; h = np.arange(N + 1) / N
        for i, J in enumerate(Jg):
            pm = exact_dist(J, N, registries=registries); ph = exact_dist(J, N, helix=True, registries=registries)
            m2, m4 = (pm * m ** 2).sum(), (pm * m ** 4).sum(); th = (ph * h).sum()
            E[N]["am"][i] = (pm * np.abs(m)).sum(); E[N]["U4"][i] = 1 - m4 / (3 * m2 * m2); E[N]["th"][i] = th
            E[N]["chi_th"][i] = N * ((ph * h * h).sum() - th * th)
    np.savez(cache, **{f"{N}_{k}": E[N][k] for N in NS for k in E[N]})
    return E


E0 = exact_table(False)          # registry-free 1D chain (the old reference)
E = exact_table(True)            # 1D chain with the registries integrated out (the updated theory)


# ---------------------------------------------------------------- MC
def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}


def blk(x, nb=20):
    b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return b.std(ddof=1) / np.sqrt(nb)


runs = {}
for log in sorted(glob.glob("logs/J*_db_s*.log")):
    done = "summary" in open(log).read()
    b = os.path.basename(log)[:-4]; J = float(re.match(r"J([\d.]+)_db", b).group(1))
    f = os.path.join("out", b + "_obs.dat")
    if not os.path.isfile(f):
        continue
    try:
        o = read_obs(f)
    except (StopIteration, ValueError, IndexError):
        continue
    if len(o["sweep"]) < 100:
        continue
    m = (o["n_R"] - o["n_L"]) / N0
    runs.setdefault(J, []).append(dict(done=done, nsw=int(o["sweep"][-1]), th=o["helicity"].mean(), am=np.abs(m).mean(), rg2=o["Rg2"].mean(), enb=o["E_nb"].mean(),
                                       etw=o["E_twist"].mean() if "E_twist" in o else np.nan,
                                       e_th=blk(o["helicity"]), e_am=blk(np.abs(m)), e_rg2=blk(o["Rg2"]), e_enb=blk(o["E_nb"]),
                                       chi_th=N0 * o["helicity"].var(), U4=1 - (m ** 4).mean() / (3 * (m ** 2).mean() ** 2)))


def cell(J, k):
    rs = runs[J]; v = np.array([r[k] for r in rs]); sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    eb = np.sqrt(sum(r["e_" + k] ** 2 for r in rs)) / len(rs) if "e_" + k in rs[0] else 0.0
    return v.mean(), max(sem, eb), len(rs)


Js = sorted(runs); nfin = sum(r["done"] for v in runs.values() for r in v); npart = sum(not r["done"] for v in runs.values() for r in v)
ntot = len(glob.glob("inputs/J*.dat"))
def filled(J): return all(r["done"] for r in runs[J])
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
TOP = (("th", r"helicity $\theta$"), ("chi_th", r"$\chi_\theta = N\,\mathrm{var}(\theta)$"), ("am", r"handedness $\langle |m| \rangle$"), ("U4", r"Binder cumulant $U_4$"))
BOT = (("rg2", r"$\langle R_g^2 \rangle$ [$a^2$]"), ("enb", r"$\langle E_{nb} \rangle$ [$k_BT$]"), ("etw", r"$\langle E_{twist} \rangle$ [$k_BT$]"))
fig, ax = plt.subplots(2, 4, figsize=(16.0, 7.6))
for c, (k, yl) in enumerate(TOP):
    x = ax[0, c]
    for N in NS:
        x.plot(Jg, E0[N][k], color="#c9c7c1", lw=0.9, zorder=1, label="1D without registries, N = 100 … 1600" if (c == 0 and N == N0) else None)
    for N in NS:
        x.plot(Jg, E[N][k], color=CN[N], lw=1.8, zorder=2, label=f"1D with registries, N = {N}" if c == 0 else None)
    for fin, lab in ((True, "MC N = 200, table + twist (finished)"), (False, "still running (partial)")):
        pts = [(J,) + cell(J, k)[:2] for J in Js if filled(J) == fin]
        if pts:
            X, V, Er = map(np.array, zip(*pts))
            x.errorbar(X, V, yerr=Er, fmt="o", ms=5.5, mfc=RED if fin else "white", mec=RED if not fin else "white", mew=1.0,
                       ecolor=RED, elinewidth=0.9, zorder=5, label=lab if c == 0 else None)
    x.set_xlim(0.5, 10.5); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl)
    x.text(0.04, 0.94, f"({'abcd'[c]})  $\\theta_0$ = {GEOM[1:]}°", transform=x.transAxes, va="top", fontsize=9, color=INK)
ax[0, 0].legend(loc="center right", handlelength=1.5, fontsize=7)
for c, (k, yl) in enumerate(BOT):
    x = ax[1, c]
    pts = [(J,) + cell(J, k)[:2] + (filled(J),) for J in Js if np.isfinite(cell(J, k)[0])]
    if pts:
        X, V, Er, Fi = map(np.array, zip(*pts))
        x.plot(X, V, "-", lw=1.0, color=RED, zorder=4)
        x.errorbar(X, V, yerr=Er, fmt="none", ecolor=RED, elinewidth=0.9, zorder=5)
        x.scatter(X, V, s=32, facecolors=np.where(Fi, RED, "white"), edgecolors=RED, linewidths=1.0, zorder=6)
    x.set_xlim(0.5, 10.5); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl)
    x.text(0.04, 0.94, f"({'efg'[c]})", transform=x.transAxes, va="top", fontsize=9, color=INK)
ax[1, 3].axis("off")
def jhalf(tab, N, k="am"):
    i = int(np.argmax(tab[N][k] >= 0.5)); return float(Jg[i]) if tab[N][k].max() >= 0.5 else float("nan")
ax[1, 3].text(0.0, 0.95, f"{GEOM}: {nfin} / {ntot} runs finished, {npart} running (hollow)\n"
              "table helix_pair_db_" + GEOM + "_win.bin, k = 2K bonded, pairs from 1-4\n\n"
              "1D theory with registries: E_HH = -J + F_HH, E_RL = +J + F_RL,\n"
              f"F_HH = {F_REG[GEOM]:.2f}, F_RL = {F_RL[GEOM]:.2f} kT (junction registry free energies)\n"
              f"<|m|> = 1/2 at N = 200: J = {jhalf(E0, N0):.2f} without, {jhalf(E, N0):.2f} with registries",
              transform=ax[1, 3].transAxes, va="top", fontsize=9)
fig.tight_layout(); fig.savefig("transitions_db.png", dpi=170); fig.savefig("transitions_db.pdf"); plt.close(fig)

print(f"{GEOM} table scan: {nfin} / {ntot} runs finished, {npart} partial")
print(f"1D theory with registries: E_HH = -J + {F_REG[GEOM]:.3f}, E_RL = +J + {F_RL[GEOM]:.3f};  <|m|> = 1/2 at N = 200: J = {jhalf(E0, N0):.2f} (no registries) / {jhalf(E, N0):.2f} (with)")
print(f"{'J':>5} {'n':>2} {'theta MC':>9} {'theory':>7} {'<|m|> MC':>9} {'theory':>7} {'U4 MC':>7} {'theory':>7} {'Rg2':>8} {'E_nb':>7} {'E_twist':>8}   (theory = 1D with registries)")
for J in Js:
    i = int(round((J - Jg[0]) / 0.125))
    t, a, u, r, e, w = (cell(J, k) for k in ("th", "am", "U4", "rg2", "enb", "etw"))
    print(f"{J:>5g} {t[2]:>2} {t[0]:>9.4f} {E[N0]['th'][i]:>7.4f} {a[0]:>9.4f} {E[N0]['am'][i]:>7.4f} {u[0]:>7.3f} {E[N0]['U4'][i]:>7.3f} {r[0]:>8.1f} {e[0]:>7.3f} {w[0]:>8.2f}"
          + ("" if filled(J) else f"   partial: {min(r['nsw'] for r in runs[J]) // 1000}k sweeps"))
