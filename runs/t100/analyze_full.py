#!/usr/bin/env python3
"""Fuller picture of the t100 E_helix scan: helicity, mean spin m = (n_R - n_L)/N and its sign statistics,
<|m|>, <S^2>/N, n_R / n_L, sense-reversal walls, energies, size, acceptance, plus histograms of m."""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def walls_from_frames(f):
    out = []
    with open(f) as fh:
        while True:
            h = fh.readline()
            if not h: break
            n = int(h); fh.readline(); sp = np.array([int(fh.readline().split()[-1]) for _ in range(n)])
            out.append((np.sum(sp[:-1] * sp[1:] < 0), np.sum((sp[:-1] == 0) != (sp[1:] == 0))))   # R/L walls, helix/coil boundaries
    return np.array(out)
def blk(x, nb=20):
    b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return x.mean(), b.std(ddof=1) / np.sqrt(nb)

N = 200; per = {}; hist = {}
for f in sorted(glob.glob("out/E*_obs.dat")):
    E = float(re.search(r"E(-?[\d.]+)_s", f).group(1)); o = read_obs(f)
    m = (o["n_R"] - o["n_L"]) / N
    d = {"helicity": o["helicity"], "m": m, "|m|": np.abs(m), "sign(m)": np.sign(m), "S2/N": N * m**2, "n_R/N": o["n_R"] / N, "n_L/N": o["n_L"] / N,
         "E_state": o["E_state"], "E_bend": o["E_bend"], "E_nb": o["E_nb"], "Re2": o["Ree"]**2, "Rg2": o["Rg2"],
         "chi_h": np.full_like(m, N * o["helicity"].var()), "acc_pos": o["acc_pos"], "acc_state": o["acc_state"], "acc_pivot": o["acc_pivot"]}
    w = walls_from_frames(f.replace("_obs.dat", "_conf.xyz"))
    d["RL walls"] = np.repeat(w[:, 0].mean(), 4); d["helix/coil boundaries"] = np.repeat(w[:, 1].mean(), 4)
    per.setdefault(E, []).append({k: blk(v) for k, v in d.items()}); hist.setdefault(E, []).append(m)
Es = sorted(per); keys = list(next(iter(per.values()))[0].keys())
mean = {k: np.array([np.mean([r[k][0] for r in per[E]]) for E in Es]) for k in keys}
err  = {k: np.array([np.sqrt(np.sum([r[k][1]**2 for r in per[E]])) / len(per[E]) for E in Es]) for k in keys}

print(f"{'E_helix':>8} {'helicity':>9} {'<m>':>8} {'+-':>6} {'<|m|>':>7} {'<sign m>':>9} {'S2/N':>7} {'n_R/N':>7} {'n_L/N':>7} {'RLwalls':>8} {'HCbound':>8} {'E_state':>8} {'E_bend':>7} {'E_nb':>6} {'Re2':>7} {'Rg2':>7} {'chi_h':>6}")
for i, E in enumerate(Es):
    print(f"{E:>8} {mean['helicity'][i]:9.4f} {mean['m'][i]:8.4f} {err['m'][i]:6.4f} {mean['|m|'][i]:7.4f} {mean['sign(m)'][i]:9.3f} {mean['S2/N'][i]:7.3f} {mean['n_R/N'][i]:7.3f} {mean['n_L/N'][i]:7.3f} "
          f"{mean['RL walls'][i]:8.2f} {mean['helix/coil boundaries'][i]:8.2f} {mean['E_state'][i]:8.1f} {mean['E_bend'][i]:7.1f} {mean['E_nb'][i]:6.2f} {mean['Re2'][i]:7.0f} {mean['Rg2'][i]:7.0f} {mean['chi_h'][i]:6.2f}")

BLUE, ORANGE, GREEN, YELLOW, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(3, 4, figsize=(16, 11), dpi=140, facecolor=BG); axs = axs.ravel()
def style(ax, title, ylab):
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8); ax.set_title(title, loc="left", fontsize=9.5, color=INK); ax.set_ylabel(ylab, color=INK); ax.set_xlabel("E_helix", color=INK)
def eb(ax, key, col, lab=None, **kw): ax.errorbar(Es, mean[key], yerr=err[key], fmt="o-", ms=4, lw=1, color=col, capsize=2, label=lab, **kw)
style(axs[0], "Helix content", "helicity"); eb(axs[0], "helicity", BLUE); axs[0].set_ylim(0, 1.02)
style(axs[1], "Mean spin  m = (n_R − n_L)/N", "⟨m⟩"); eb(axs[1], "m", ORANGE); axs[1].axhline(0, color=INK2, lw=0.8)
style(axs[2], "Sign statistics of m", ""); eb(axs[2], "|m|", GREEN, "⟨|m|⟩"); eb(axs[2], "sign(m)", YELLOW, "⟨sign m⟩"); axs[2].axhline(0, color=INK2, lw=0.8); axs[2].legend(frameon=False, fontsize=8)
style(axs[3], "Chiral order fluctuations", "⟨S²⟩/N = N⟨m²⟩"); eb(axs[3], "S2/N", GREEN)
style(axs[4], "R and L content", "fraction"); eb(axs[4], "n_R/N", BLUE, "R"); eb(axs[4], "n_L/N", ORANGE, "L"); axs[4].legend(frameon=False, fontsize=8)
style(axs[5], "Domain walls per chain", "count"); eb(axs[5], "RL walls", ORANGE, "R|L sense reversals"); eb(axs[5], "helix/coil boundaries", BLUE, "helix|coil boundaries"); axs[5].legend(frameon=False, fontsize=8)
style(axs[6], "Energies", "kT"); eb(axs[6], "E_state", BLUE, "E_state (pairs + E_helix)"); eb(axs[6], "E_bend", GREEN, "E_bend"); eb(axs[6], "E_nb", ORANGE, "E_nb (steric)"); axs[6].legend(frameon=False, fontsize=8)
style(axs[7], "Helicity susceptibility", "χ = N·var(helicity)"); eb(axs[7], "chi_h", BLUE)
style(axs[8], "Chain size", ""); eb(axs[8], "Re2", ORANGE, "⟨Rₑ²⟩"); axs[8].errorbar(Es, 6 * mean["Rg2"], yerr=6 * err["Rg2"], fmt="s-", ms=4, lw=1, mfc=BG, color=ORANGE, capsize=2, label="6⟨Rg²⟩"); axs[8].legend(frameon=False, fontsize=8)
style(axs[9], "Acceptance rates", ""); eb(axs[9], "acc_pos", BLUE, "position"); eb(axs[9], "acc_state", GREEN, "state"); eb(axs[9], "acc_pivot", ORANGE, "pivot"); axs[9].legend(frameon=False, fontsize=8)
# histograms of m at three E values
for ax, E in zip((axs[10], axs[11]), (max(Es, key=lambda e: abs(mean["helicity"][Es.index(e)] - 0.5) * -1), Es[0])):
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    m = np.concatenate(hist[E]); ax.hist(m, bins=60, color=BLUE, alpha=0.8, density=True)
    ax.set_title(f"distribution of m at E_helix = {E}  (helicity {mean['helicity'][Es.index(E)]:.2f})", loc="left", fontsize=9.5, color=INK)
    ax.set_xlabel("m = (n_R − n_L)/N", color=INK); ax.set_ylabel("P(m)", color=INK); ax.tick_params(colors=INK2, labelsize=8); ax.axvline(0, color=INK2, lw=0.8)
fig.suptitle(f"t100 scan: N = {N}, E_HH = +9.5, E_CH = +8, E_RL = +11, sterics on — 3 seeds × 2×10⁵ sweeps", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig("t100_full.png")
