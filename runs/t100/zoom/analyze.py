#!/usr/bin/env python3
"""Boundary zoom: helicity vs J0 (0.25 steps) per E_helix, MC with sterics vs exact no-steric; midpoints J0*."""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["x", "../sample_data.dat"]
import importlib.util
spec = importlib.util.spec_from_file_location("ex", os.path.abspath("../exact_fss.py")); ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)
def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
cells = {}
for f in sorted(glob.glob("out/J*_E*_obs.dat")):
    m = re.search(r"J(-?[\d.]+)_E(-?\d+)_s(\d+)_obs", f); J, E = float(m.group(1)), int(m.group(2))
    o = read_obs(f)
    if len(o["helicity"]) < 40: continue
    nb = 20; b = o["helicity"][: len(o["helicity"]) // nb * nb].reshape(nb, -1).mean(1)
    cells.setdefault((J, E), []).append((o["helicity"].mean(), b.std(ddof=1) / np.sqrt(nb)))
Es = sorted({E for J, E in cells}, reverse=True)
def mid(xs, ys):
    xs, ys = np.array(xs), np.array(ys)
    k = np.argsort(xs); xs, ys = xs[k], ys[k]
    for i in range(len(xs) - 1):
        if (ys[i] - 0.5) * (ys[i + 1] - 0.5) <= 0:
            return xs[i] + (0.5 - ys[i]) * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i])
    return np.nan
INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
cols = dict(zip(Es, ("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")))
fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=160, facecolor=BG)
ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
print(f"{'E0':>5} {'J0* (MC, sterics)':>18} {'J0* (exact no-ster)':>20} {'steric shift':>13}")
for E in Es:
    pts = sorted((J, np.mean([r[0] for r in v]), np.sqrt(np.sum([r[1]**2 for r in v])) / len(v)) for (J, EE), v in cells.items() if EE == E)
    Jline = np.linspace(pts[0][0] - 0.4, pts[-1][0] + 0.4, 60); hx = []
    for Jv in Jline:
        ex.J0 = float(Jv); hx.append(ex.hel_N(float(E), 200))
    ax.plot(Jline, hx, "--", lw=1.1, color=cols[E], alpha=0.85)
    ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[p[2] for p in pts], fmt="o-", ms=4.5, lw=1.2, color=cols[E], capsize=2.5, label=f"E0 = {E}")
    Jmc = mid([p[0] for p in pts], [p[1] for p in pts]); Jex = mid(Jline, hx)
    print(f"{E:>5} {Jmc:>18.3f} {Jex:>20.3f} {Jmc - Jex:>+13.3f}")
ax.set_xlabel("J0   (E_HH = −J0)", color=INK); ax.set_ylabel("helicity", color=INK)
ax.set_title("Boundary zoom: MC with sterics (solid) vs exact no-steric N=200 (dashed)", loc="left", fontsize=9.5, color=INK)
ax.tick_params(colors=INK2, labelsize=8.5); ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("boundary_zoom.png")
