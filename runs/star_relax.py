#!/usr/bin/env python3
"""Relaxation check of the STAR J scan (user 2026-09-23: "but I think it is growing"): helicity and signed handedness
of the whole star versus production sweep for every run that has data, block-averaged over BLOCK sweeps, one colour
per J; dashed = the 1D theory (independent arms with registries, exact_cache_star_reg.npz from star_plot.py).
  cd runs/t45/star50 && ~/venv/bin/python ../../star_relax.py     -> relax_star.png + drift per 1000 sweeps
"""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLOCK = 250
Jg = np.arange(0.5, 10.501, 0.125)
if not os.path.exists("exact_cache_star_reg.npz"):
    sys.exit("run star_plot.py first (needs exact_cache_star_reg.npz)")
TH = np.load("exact_cache_star_reg.npz")["th"]
def read_kv(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
    return d
kv = read_kv(sorted(glob.glob("inputs/J*_s*.dat"))[0]); N0 = int(kv["N"]) * int(kv["n_arms"]); LOG = int(kv["log_every"])
runs = []
for f in sorted(glob.glob("out/J*_s*_obs.dat")):
    rows = np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith("#")])
    if rows.ndim < 2 or len(rows) < 2 * BLOCK // LOG:
        continue
    b = os.path.basename(f)[:-8]; J = float(re.match(r"J([\d.]+)_s", b).group(1)); s = int(b.split("_s")[1])
    nb = BLOCK // LOG; k = len(rows) // nb
    sw = rows[: k * nb, 0].reshape(k, nb).mean(1); th = rows[: k * nb, 3].reshape(k, nb).mean(1)
    m = ((rows[: k * nb, 1] - rows[: k * nb, 2]) / N0).reshape(k, nb).mean(1)
    runs.append((J, s, sw, th, m))
if not runs:
    sys.exit("no run has two blocks yet")
Js = sorted({r[0] for r in runs}); cmap = plt.get_cmap("viridis"); col = {J: cmap(i / max(1, len(Js) - 1)) for i, J in enumerate(Js)}
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.0))
tmax = max(r[2][-1] for r in runs)
for J, s, sw, th, m in runs:
    ax[0].plot(sw, th, "-" if s == 1 else ("--" if s == 2 else ":"), color=col[J], lw=1.2, label=f"J = {J:g}" if s == 1 else None)
    ax[1].plot(sw, m, "-" if s == 1 else ("--" if s == 2 else ":"), color=col[J], lw=1.2)
for J in Js:
    i = int(round((J - Jg[0]) / 0.125))
    ax[0].plot([0, tmax], [TH[i]] * 2, ls=(0, (2, 3)), color=col[J], lw=0.9)
ax[0].set_xlabel("production sweep"); ax[0].set_ylabel(r"helicity $\theta$ (block mean over %d sweeps)" % BLOCK)
ax[0].text(0.03, 0.05, "dotted: 1D theory, independent arms with registries\nsolid / dashed / dotted lines: seed 1 / 2 / 3", transform=ax[0].transAxes, fontsize=8, va="bottom")
ax[0].legend(loc="center right", ncol=2, handlelength=1.5)
ax[1].axhline(0, color="#898781", lw=0.8); ax[1].set_xlabel("production sweep"); ax[1].set_ylabel("signed handedness $m$ of the whole star")
ax[1].set_ylim(-0.3, 0.3)
fig.tight_layout(); fig.savefig("relax_star.png", dpi=170); plt.close(fig)
print(f"{'run':>8} {'sweeps':>6} {'theta first':>11} {'theta last':>10} {'theory':>7} {'drift / 1000 sw':>15} {'|m| mean':>8} {'m sign changes':>14}")
for J, s, sw, th, m in runs:
    i = int(round((J - Jg[0]) / 0.125)); half = len(th) // 2
    drift = (th[half:].mean() - th[:half].mean()) / ((sw[half:].mean() - sw[:half].mean()) / 1000.0)
    print(f"J{J:g}_s{s:<4} {int(sw[-1]) + BLOCK // 2:>6} {th[0]:>11.3f} {th[-1]:>10.3f} {TH[i]:>7.3f} {drift:>+15.4f} {np.abs(m).mean():>8.3f} {int((np.sign(m[1:]) != np.sign(m[:-1])).sum()):>14}")
