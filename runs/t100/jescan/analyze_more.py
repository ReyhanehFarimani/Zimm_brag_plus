#!/usr/bin/env python3
"""More figures for the J0 x E_helix scan: phase-diagram heatmaps with the exact no-steric boundary,
handedness statistics (<m>, <|m|>, P(m), m(t)) across the plane."""
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

N = 200; cells = {}; series = {}
for f in sorted(glob.glob("out/J*_E*_obs.dat")):
    mm = re.search(r"J(-?\d+)_E(-?\d+)_s(\d+)_obs", f); J, E = int(mm.group(1)), int(mm.group(2))
    o = read_obs(f)
    if len(o["helicity"]) < 40: continue
    m = (o["n_R"] - o["n_L"]) / N
    c = cells.setdefault((J, E), dict(hel=[], S2N=[], Re2=[], Enb=[], am=[], mn=[]))
    c["hel"].append(o["helicity"].mean()); c["S2N"].append(N * (m**2).mean()); c["Re2"].append((o["Ree"]**2).mean())
    c["Enb"].append(o["E_nb"].mean()); c["am"].append(np.abs(m).mean()); c["mn"].append(m.mean())
    series.setdefault((J, E), []).append((o["sweep"], m))
Js = sorted({J for J, E in cells}); Es = sorted({E for J, E in cells})
grid = {k: np.array([[np.mean(cells[(J, E)][k]) if (J, E) in cells else np.nan for J in Js] for E in Es]) for k in ("hel", "S2N", "Re2", "Enb", "am", "mn")}

# exact no-steric helicity = 0.5 boundary on a fine grid
Jf = np.linspace(min(Js), max(Js), 41); Ef = np.linspace(min(Es), max(Es), 41)
H = np.zeros((len(Ef), len(Jf)))
for i, Ev in enumerate(Ef):
    for j, Jv in enumerate(Jf):
        ex.J0 = float(Jv); H[i, j] = ex.hel_inf(float(Ev))

INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(1, 4, figsize=(17, 4.2), dpi=150, facecolor=BG)
maps = (("hel", "helicity", "Blues", None), ("S2N", "log₁₀ ⟨S²⟩/N", "Purples", np.log10),
        ("Re2", "⟨Rₑ²⟩ / 1000", "Oranges", lambda x: x / 1000), ("Enb", "⟨E_nb⟩ [kT]", "Greens", None))
for ax, (key, lab, cmap, tr) in zip(axs, maps):
    Z = grid[key] if tr is None else tr(grid[key])
    pc = ax.pcolormesh(Js, Es, Z, cmap=cmap, shading="nearest", edgecolors="white", linewidth=1.2)
    for i, E in enumerate(Es):
        for j, J in enumerate(Js):
            v = Z[i, j]
            ax.text(J, E, f"{v:.2f}" if abs(v) < 100 else f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                    color=("white" if (v - np.nanmin(Z)) / (np.nanmax(Z) - np.nanmin(Z) + 1e-12) > 0.6 else INK))
    ax.contour(Jf, Ef, H, levels=[0.5], colors=INK, linewidths=1.4, linestyles="--")
    ax.set_xlabel("J0  (E_HH = −J0)", color=INK); ax.set_ylabel("E_helix", color=INK)
    ax.set_title(lab + "   (-- exact no-steric boundary)", loc="left", fontsize=9, color=INK)
    fig.colorbar(pc, ax=ax, shrink=0.85).ax.tick_params(labelsize=7)
fig.suptitle("t100 phase plane: N = 200, sterics on, E_CH = +8, E_RL = +11", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig("je_maps.png")

# handedness figure: <|m|> map, P(m) histograms, m(t) traces
fig = plt.figure(figsize=(15, 8.4), dpi=140, facecolor=BG)
ax = fig.add_subplot(231)
Z = np.sqrt(N) * grid["am"]
pc = ax.pcolormesh(Js, Es, Z, cmap="Purples", shading="nearest", edgecolors="white", linewidth=1.2)
for i, E in enumerate(Es):
    for j, J in enumerate(Js): ax.text(J, E, f"{Z[i,j]:.2f}", ha="center", va="center", fontsize=6.5, color=("white" if Z[i,j] > 0.6*np.nanmax(Z) else INK))
ax.set_xlabel("J0", color=INK); ax.set_ylabel("E_helix", color=INK); fig.colorbar(pc, ax=ax, shrink=0.85).ax.tick_params(labelsize=7)
ax.set_title("√N·⟨|m|⟩  (0.67 = disordered Gaussian; ≫1 = homochiral domains)", loc="left", fontsize=8.5, color=INK)
picks = [(-6, -9), (-7, -11), (-8, -10), (-10, -11), (-12, -12)]
for k, (J, E) in enumerate(picks[:5]):
    ax = fig.add_subplot(2, 3, k + 2)
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    if (J, E) in series:
        m = np.concatenate([s[1] for s in series[(J, E)]])
        ax.hist(m, bins=51, range=(-1, 1), density=True, color="#2a78d6", alpha=0.85)
        s0, m0 = series[(J, E)][0]
        ax2 = ax.twinx(); ax2.plot(np.linspace(-1, 1, len(m0)), np.nan * m0)  # keep scales sane
        ax.set_title(f"P(m) at J0={J}, E0={E}   (hel={np.mean(cells[(J,E)]['hel']):.2f}, ⟨S²⟩/N={np.mean(cells[(J,E)]['S2N']):.1f})", loc="left", fontsize=8.5, color=INK)
        ax.set_xlabel("m = (n_R − n_L)/N", color=INK); ax.set_ylabel("P(m)", color=INK); ax.axvline(0, color=INK2, lw=0.8)
    ax.tick_params(colors=INK2, labelsize=7.5)
fig.suptitle("Handedness across the plane — histograms of the chiral order parameter m", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig("je_handedness.png")

# m(t) traces for the strongly chiral cells
fig, axs = plt.subplots(3, 1, figsize=(11, 6.5), dpi=140, sharex=True, facecolor=BG)
for ax, (J, E) in zip(axs, [(-6, -9), (-7, -11), (-8, -10)]):
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for (s, m), col in zip(series[(J, E)], ("#2a78d6", "#eb6834")):
        ax.plot(s / 1e5, m, lw=0.7, color=col)
    ax.axhline(0, color=INK2, lw=0.8); ax.set_ylabel("m", color=INK); ax.set_ylim(-1.05, 1.05)
    ax.text(0.01, 0.9, f"J0={J}, E0={E}", transform=ax.transAxes, fontsize=9, color=INK)
    ax.tick_params(colors=INK2, labelsize=8)
axs[-1].set_xlabel("sweeps [×10⁵]", color=INK)
axs[0].set_title("Chiral order parameter over the run (two seeds per panel) — slow sense dynamics in the homochiral regime", loc="left", fontsize=10, color=INK)
fig.tight_layout(); fig.savefig("je_m_traces.png")
print("wrote je_maps.png, je_handedness.png, je_m_traces.png")
