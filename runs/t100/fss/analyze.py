#!/usr/bin/env python3
"""Finite-size scaling of the t100 crossover WITH sterics: helicity and chi = N var(h) vs E_helix for N = 50, 100, 200."""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".."); sys.argv = ["x", "../sample_data.dat"]
import importlib.util
spec = importlib.util.spec_from_file_location("exact_fss", "../exact_fss.py"); ex = importlib.util.module_from_spec(spec); spec.loader.exec_module(ex)

def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def tau_int(x):
    x = x - x.mean(); n = len(x); f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); ac /= ac[0]; t = 0.5
    for M in range(1, n // 2):
        t += ac[M]
        if M > 6 * t: break
    return t
data = {}
files = [(int(re.search(r"N(\d+)_E(-?[\d.]+)_s", f).group(1)), float(re.search(r"N(\d+)_E(-?[\d.]+)_s", f).group(2)), f) for f in glob.glob("out/N*_obs.dat")]
files += [(200, float(re.search(r"E(-?[\d.]+)_s", f).group(1)), f) for f in glob.glob("../out/E*_obs.dat")]
for N, E, f in files:
    o = read_obs(f); h = o["helicity"]
    if len(h) < 40: continue
    nb = 20; b = h[: len(h) // nb * nb].reshape(nb, -1)
    hm = h.mean(); chi = N * h.var(); chib = N * b.var(1)
    data.setdefault((N, E), []).append((hm, b.mean(1).std(ddof=1) / np.sqrt(nb), chi, chib.std(ddof=1) / np.sqrt(nb), len(h) / (2 * tau_int(h))))
Ns = sorted({N for N, E in data}); Es = sorted({E for N, E in data})
print(f"{'E_helix':>8} | " + " | ".join(f"{'N='+str(N):>22}" for N in Ns) + "   (helicity +- , chi=N var(h))")
res = {}
for E in Es:
    row = []
    for N in Ns:
        if (N, E) in data:
            r = np.array(data[(N, E)]); hm = r[:, 0].mean(); he = np.sqrt((r[:, 1] ** 2).sum()) / len(r); cm = r[:, 2].mean(); ce = np.sqrt((r[:, 3] ** 2).sum()) / len(r)
            res[(N, E)] = (hm, he, cm, ce); row.append(f"{hm:6.3f}±{he:5.3f} {cm:6.2f}±{ce:4.2f}")
        else: row.append(" " * 22)
    print(f"{E:>8} | " + " | ".join(row))
print("\nexact no-steric chi peak: N=50 2.39, N=100 2.79, N=200 2.99, N=400 3.09, N=inf ~3.2 (saturates); a true transition would grow ~N")
for N in Ns:
    pts = sorted((E, res[(N, E)][2]) for (n, E) in res if n == N)
    if pts: E0, c0 = max(pts, key=lambda p: p[1]); print(f"  with sterics, N = {N}: chi_max = {c0:.2f} at E_helix = {E0}")
BLUE, ORANGE, GREEN, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
cols = {50: GREEN, 100: ORANGE, 200: BLUE, 400: "#4a3aa7"}
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.set_xlabel("E_helix", color=INK); ax.set_xlim(-11.6, -9.9)
Eline = np.linspace(-11.5, -10.0, 60)
for N in Ns:
    pts = sorted((E, res[(N, E)]) for (n, E) in res if n == N)
    axs[0].errorbar([p[0] for p in pts], [p[1][0] for p in pts], yerr=[p[1][1] for p in pts], fmt="o-", ms=4, lw=1, color=cols.get(N, INK), capsize=2, label=f"MC with sterics, N={N}")
    axs[0].plot(Eline, [ex.hel_N(E, N) for E in Eline], "--", lw=1, color=cols.get(N, INK), alpha=0.7)
    axs[1].errorbar([p[0] for p in pts], [p[1][2] for p in pts], yerr=[p[1][3] for p in pts], fmt="o-", ms=4, lw=1, color=cols.get(N, INK), capsize=2, label=f"MC with sterics, N={N}")
    axs[1].plot(Eline, [ex.hel_N(E, N, True)[1] for E in Eline], "--", lw=1, color=cols.get(N, INK), alpha=0.7)
axs[0].plot(Eline, [ex.hel_inf(E) for E in Eline], "-", lw=1.4, color=INK2, label="exact no sterics, N=∞ (dashed: finite N)")
axs[0].set_ylabel("helicity", color=INK); axs[0].set_title("Helicity vs E_helix for several N", loc="left", fontsize=10, color=INK); axs[0].legend(frameon=False, fontsize=7.5)
axs[1].set_ylabel("χ = N·var(helicity)", color=INK); axs[1].set_title("Susceptibility: saturates with N ⇒ crossover", loc="left", fontsize=10, color=INK); axs[1].legend(frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig("t100_fss.png")
