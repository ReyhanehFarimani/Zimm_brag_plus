#!/usr/bin/env python3
"""Transition figure of the STAR J scan (runs/<geom>/star<n>, star_gen.sh), same layout as jscan_db_plot.py:
helicity theta(J), chi_theta = N var(theta), handedness <|m|>(J), Binder U4(J) in the top row; <Rg^2>, <E_nb>,
<E_twist> in the bottom row.  theta, m are of the WHOLE star (all n_arms x N residues, as in the obs file).
Reference curves (transfer matrix of ../../t100/exact_fss.py on ../sample_data.dat, no non-bonded interactions):
  blue   = n_arms INDEPENDENT arms of N residues, 1D theory with the registries integrated out (E_HH = -J + F_HH,
           E_RL = +J + F_RL as in jscan_db_plot.py): theta and chi are those of one arm, <|m|> and U4 come from the
           n_arms-fold convolution of the one-arm m distribution (the star's m is the arm average);
  dashed = ONE arm alone (what a single N-residue chain would show for <|m|>, U4);
  grey   = independent arms without registries.
Red points: MC star, filled = finished runs, hollow = still running; only the LAST half of each run's production rows
enters the averages (user 2026-09-23: "when plotting consider only the last steps"; the star relaxes slowly), at least
20 rows; errors = max(SEM over seeds, combined 20-block error).  Everything the MC does beyond the blue
curve is inter-arm (crowding on the core) physics.
  cd runs/t45/star50 && python ../../star_plot.py      -> transitions_star.png / .pdf + numbers
"""
import contextlib
import glob
import importlib.util
import io
import os
import re
import sys
from math import erf, sqrt, pi, log

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.getcwd()
GEOM = os.path.basename(os.path.dirname(HERE))
Jg = np.arange(0.5, 10.501, 0.125)
LAST = 0.5          # fraction of each run's production rows that enter the averages, counted from the end (the star relaxes slowly)


def read_kv(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
    return d


INPUTS = sorted(glob.glob("inputs/J*_s*.dat"))
if not INPUTS:
    sys.exit("no inputs/J*_s*.dat here: run from runs/<geom>/star<n>")
kv = read_kv(INPUTS[0])
NARM, NARMS, CORE = int(kv["N"]), int(kv["n_arms"]), float(kv.get("core_radius", 0.0))
BOX = float(kv.get("box", 0.0))                      # periodic box of free chains (runs/box_gen.sh) instead of a star
SYS = (f"box {BOX:g} a, {NARMS} free chains × {NARM}" if BOX > 0 else f"star {NARMS} arms × {NARM}, core {CORE:g} a")
N0 = NARM * NARMS
NSW = int(kv["n_sweeps"])

KAPPA = {"t45": 6.6014, "t100": 1.8177}                      # R.R / L.L junction (kT/rad^2)
KAPPA_RL = {"t45": 2 * 2.465, "t100": 2 * 0.108}              # R.L junction
def f_reg(k): return -log(sqrt(2 * pi / k) * erf(pi * sqrt(k / 2)) / (2 * pi))
F_REG = {g: f_reg(k) for g, k in KAPPA.items()}
F_RL = {g: f_reg(k) for g, k in KAPPA_RL.items()}
BLUE, LBLUE, GREY = "#1c5cab", "#5598e7", "#c9c7c1"
RED, INK = "#c0392f", "#0b0b0b"

# ---------------------------------------------------------------- exact 1D reference
sys.argv = ["x", "../sample_data.dat"]
spec = importlib.util.spec_from_file_location("ex", os.path.abspath("../../t100/exact_fss.py")); ex = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(ex)
PAIRS = sorted(ex.idx, key=ex.idx.get); SECOND = np.array([b for (a, b) in PAIRS])


def exact_dist(J, N, helix=False, registries=False):
    """distribution of n_R + n_L (helix) or n_R - n_L over one chain of N residues"""
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


def star_dist(p1, n):
    """n-fold convolution of the one-arm distribution of n_R - n_L: the star's n_R - n_L over [-n N, n N]"""
    p = p1.copy()
    for _ in range(n - 1):
        p = np.convolve(p, p1)
    return p / p.sum()


def exact_table(registries):
    keys = ("th", "chi_th", "am", "U4", "am1", "U41")   # am/U4: star of NARMS independent arms; am1/U41: one arm
    cache = f"exact_cache_star{'_reg' if registries else ''}.npz"
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime("../sample_data.dat"):
        z = np.load(cache)
        if len(z["th"]) == len(Jg) and int(z["NARM"]) == NARM and int(z["NARMS"]) == NARMS:
            return {k: z[k] for k in keys}
    E = {k: np.zeros(len(Jg)) for k in keys}
    m1 = np.arange(-NARM, NARM + 1) / NARM; h = np.arange(NARM + 1) / NARM; ms = np.arange(-N0, N0 + 1) / N0
    for i, J in enumerate(Jg):
        p1 = exact_dist(J, NARM, registries=registries); ph = exact_dist(J, NARM, helix=True, registries=registries)
        ps = star_dist(p1, NARMS)
        th = (ph * h).sum(); E["th"][i] = th; E["chi_th"][i] = NARM * ((ph * h * h).sum() - th * th)
        for k, p, m in (("", ps, ms), ("1", p1, m1)):
            m2, m4 = (p * m ** 2).sum(), (p * m ** 4).sum()
            E["am" + k][i] = (p * np.abs(m)).sum(); E["U4" + k][i] = 1 - m4 / (3 * m2 * m2)
    np.savez(cache, NARM=NARM, NARMS=NARMS, **E)
    return E


E0 = exact_table(False)          # independent arms, registry-free
E = exact_table(True)            # independent arms, registries integrated out (the theory)


# ---------------------------------------------------------------- MC
def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}


def blk(x, nb=20):
    b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return b.std(ddof=1) / np.sqrt(nb)


runs = {}
for log in sorted(glob.glob("logs/J*_s*.log")):
    done = "summary" in open(log).read()
    b = os.path.basename(log)[:-4]; J = float(re.match(r"J([\d.]+)_s", b).group(1))
    f = os.path.join("out", b + "_obs.dat")
    if not os.path.isfile(f):
        continue
    try:
        o = read_obs(f)
    except (StopIteration, ValueError, IndexError):
        continue
    if o["sweep"].ndim == 0 or len(o["sweep"]) < 20:
        continue
    nrow = len(o["sweep"]); keep = max(20, int(round(LAST * nrow)))     # only the last steps (user 2026-09-23)
    o = {k: v[nrow - keep:] for k, v in o.items()}
    m = (o["n_R"] - o["n_L"]) / N0
    runs.setdefault(J, []).append(dict(done=done, nsw=int(o["sweep"][-1]), nrow=keep, th=o["helicity"].mean(), am=np.abs(m).mean(), rg2=o["Rg2"].mean(), enb=o["E_nb"].mean(),
                                       etw=o["E_twist"].mean() if "E_twist" in o else np.nan,
                                       e_th=blk(o["helicity"]), e_am=blk(np.abs(m)), e_rg2=blk(o["Rg2"]), e_enb=blk(o["E_nb"]),
                                       chi_th=N0 * o["helicity"].var(), U4=1 - (m ** 4).mean() / (3 * (m ** 2).mean() ** 2)))


def cell(J, k):
    rs = runs[J]; v = np.array([r[k] for r in rs]); sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    eb = np.sqrt(sum(r["e_" + k] ** 2 for r in rs)) / len(rs) if "e_" + k in rs[0] else 0.0
    return v.mean(), max(sem, eb), len(rs)


Js = sorted(runs); nfin = sum(r["done"] for v in runs.values() for r in v); npart = sum(not r["done"] for v in runs.values() for r in v)
ntot = len(INPUTS)
def filled(J): return all(r["done"] for r in runs[J])
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
TOP = (("th", r"helicity $\theta$"), ("chi_th", r"$\chi_\theta = N\,\mathrm{var}(\theta)$"), ("am", r"handedness $\langle |m| \rangle$"), ("U4", r"Binder cumulant $U_4$"))
BOT = (("rg2", r"$\langle R_g^2 \rangle$ [$a^2$]"), ("enb", r"$\langle E_{nb} \rangle$ [$k_BT$]"), ("etw", r"$\langle E_{twist} \rangle$ [$k_BT$]"))
STAR = f"{NARMS} arms × {NARM}"
fig, ax = plt.subplots(2, 4, figsize=(16.0, 7.6))
for c, (k, yl) in enumerate(TOP):
    x = ax[0, c]
    x.plot(Jg, E0[k], color=GREY, lw=0.9, zorder=1, label="independent arms, no registries" if c == 0 else None)
    if k in ("am", "U4"):
        x.plot(Jg, E[k + "1"], color=LBLUE, lw=1.2, ls="--", zorder=2, label=f"one arm alone (N = {NARM}), 1D with registries" if c == 2 else None)
    x.plot(Jg, E[k], color=BLUE, lw=1.8, zorder=3, label="independent arms, 1D with registries" if c == 0 else None)
    for fin, lab in ((True, ("MC box" if BOX > 0 else "MC star") + ", table + twist (finished)"), (False, "still running (partial)")):
        pts = [(J,) + cell(J, k)[:2] for J in Js if filled(J) == fin]
        if pts:
            X, V, Er = map(np.array, zip(*pts))
            x.errorbar(X, V, yerr=Er, fmt="o", ms=5.5, mfc=RED if fin else "white", mec=RED if not fin else "white", mew=1.0,
                       ecolor=RED, elinewidth=0.9, zorder=5, label=lab if c == 0 else None)
    x.set_xlim(0.5, 10.5); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl)
    x.text(0.04, 0.94, f"({'abcd'[c]})  $\\theta_0$ = {GEOM[1:]}°", transform=x.transAxes, va="top", fontsize=9, color=INK)
ax[0, 0].legend(loc="lower right", handlelength=1.5, fontsize=7)
ax[0, 2].legend(loc="upper left", bbox_to_anchor=(0.02, 0.9), handlelength=1.5, fontsize=7)
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
def jhalf(v):
    i = int(np.argmax(v >= 0.5)); return float(Jg[i]) if v.max() >= 0.5 else float("nan")
ax[1, 3].text(0.0, 0.95, f"{GEOM} {SYS}: {nfin} / {ntot} runs finished, {npart} running (hollow)\n"
              f"{NSW // 1000}k production sweeps per run, averages over the LAST {LAST:.0%} of each run; table helix_pair_db_{GEOM}_win.bin\n\n"
              "1D theory with registries: E_HH = -J + F_HH, E_RL = +J + F_RL,\n"
              f"F_HH = {F_REG[GEOM]:.2f}, F_RL = {F_RL[GEOM]:.2f} kT (junction registry free energies)\n"
              f"theta = 1/2 at J = {jhalf(E['th']):.2f};  <|m|> = 1/2 at J = {jhalf(E['am1']):.2f} for one arm of {NARM},\n"
              f"the star of {NARMS} INDEPENDENT arms saturates at <|m|> = {E['am'][-1]:.3f} (random arm signs, ~ sqrt(2 / pi n_arms));\n"
              "MC points above that blue curve = chirality transfer between arms",
              transform=ax[1, 3].transAxes, va="top", fontsize=9)
fig.tight_layout(); fig.savefig("transitions_star.png", dpi=170); fig.savefig("transitions_star.pdf"); plt.close(fig)

print(f"{GEOM} {SYS}: {nfin} / {ntot} runs finished, {npart} partial; averages over the last {LAST:.0%} of each run")
print(f"1D theory with registries: theta = 1/2 at J = {jhalf(E['th']):.2f}; <|m|> = 1/2 at J = {jhalf(E['am1']):.2f} for one arm; {NARMS} independent arms saturate at <|m|> = {E['am'][-1]:.3f}")
print(f"{'J':>5} {'n':>2} {'theta MC':>9} {'theory':>7} {'<|m|> MC':>9} {'theory':>7} {'U4 MC':>7} {'theory':>7} {'Rg2':>8} {'E_nb':>7} {'E_twist':>8}   (theory = independent arms with registries)")
for J in Js:
    i = int(round((J - Jg[0]) / 0.125))
    t, a, u, r, e, w = (cell(J, k) for k in ("th", "am", "U4", "rg2", "enb", "etw"))
    print(f"{J:>5g} {t[2]:>2} {t[0]:>9.4f} {E['th'][i]:>7.4f} {a[0]:>9.4f} {E['am'][i]:>7.4f} {u[0]:>7.3f} {E['U4'][i]:>7.3f} {r[0]:>8.1f} {e[0]:>7.3f} {w[0]:>8.2f}"
          + ("" if filled(J) else f"   partial: {min(r['nsw'] for r in runs[J]) // 1000}k sweeps"))
