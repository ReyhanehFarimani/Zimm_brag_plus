#!/usr/bin/env python3
"""2D scan J0 x E_helix of the t100 model with sterics: helicity, <S^2>/N, R|L walls, Re^2 vs J0 for each E_helix,
against the exact no-steric transfer matrix. E_HH = -J0, E_CH = +8, E_RL = +11 fixed."""
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
def blk(x, nb=20):
    b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return x.mean(), b.std(ddof=1) / np.sqrt(nb)

N = 200; data = {}
for f in sorted(glob.glob("out/J*_E*_obs.dat")):
    m = re.search(r"J(-?\d+)_E(-?\d+)_s(\d+)_obs", f); J, E = int(m.group(1)), int(m.group(2))
    o = read_obs(f)
    if len(o["helicity"]) < 40: continue
    mm = (o["n_R"] - o["n_L"]) / N
    d = dict(hel=blk(o["helicity"]), S2N=blk(N * mm**2), Re2=blk(o["Ree"]**2), Enb=blk(o["E_nb"]), chi=(N * o["helicity"].var(), 0))
    data.setdefault((J, E), []).append(d)
Js = sorted({J for J, E in data}); Es = sorted({E for J, E in data}, reverse=True)
res = {}
for k, v in data.items():
    res[k] = {q: (np.mean([r[q][0] for r in v]), np.sqrt(np.sum([r[q][1]**2 for r in v])) / len(v)) for q in v[0]}
print(f"N = {N}, sterics on; E_HH = -J0, E_CH = +8, E_RL = +11.  cells: helicity (MC with sterics / exact no-steric)")
print(f"{'J0':>4} | " + " | ".join(f"{'E0='+str(E):>17}" for E in Es))
exact = {}
for J in Js:
    row = []
    for E in Es:
        # exact no-steric: rebuild with this J0 (exact_fss module reads J0 from the file; recompute via build with E)
        ex.J0 = float(J); h = ex.hel_N(float(E), N)
        exact[(J, E)] = h
        row.append(f"{res[(J,E)]['hel'][0]:7.3f} / {h:7.3f}" if (J, E) in res else " " * 17)
    print(f"{J:>4} | " + " | ".join(row))
print("\nper-cell details: helicity +-, <S^2>/N, chi, <Re^2>, E_nb")
for E in Es:
    for J in Js:
        if (J, E) in res:
            r = res[(J, E)]
            print(f"  E0={E:>4} J0={J:>4}: hel = {r['hel'][0]:6.3f}±{r['hel'][1]:.3f}  S2/N = {r['S2N'][0]:8.3f}  chi = {r['chi'][0]:6.2f}  Re2 = {r['Re2'][0]:7.0f}  E_nb = {r['Enb'][0]:6.2f}")

BLUE, ORANGE, GREEN, VIOLET, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
cols = dict(zip(Es, (BLUE, ORANGE, GREEN, VIOLET)))
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.3), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.set_xlabel("J0   (E_HH = −J0)", color=INK)
Jline = np.linspace(min(Js) - 0.3, max(Js) + 0.3, 80)
for E in Es:
    pts = sorted((J, res[(J, E)]) for (J, EE) in res if EE == E)
    hl = []
    for Jv in Jline:
        ex.J0 = float(Jv); hl.append(ex.hel_inf(float(E)))
    axs[0].plot(Jline, hl, "--", lw=1, color=cols[E], alpha=0.8)
    axs[0].errorbar([p[0] for p in pts], [p[1]["hel"][0] for p in pts], yerr=[p[1]["hel"][1] for p in pts], fmt="o-", ms=4, lw=1, color=cols[E], capsize=2, label=f"E0 = {E}")
    axs[1].errorbar([p[0] for p in pts], [p[1]["S2N"][0] for p in pts], yerr=[p[1]["S2N"][1] for p in pts], fmt="o-", ms=4, lw=1, color=cols[E], capsize=2, label=f"E0 = {E}")
    axs[2].errorbar([p[0] for p in pts], [p[1]["Re2"][0] for p in pts], yerr=[p[1]["Re2"][1] for p in pts], fmt="o-", ms=4, lw=1, color=cols[E], capsize=2, label=f"E0 = {E}")
axs[0].set_ylabel("helicity", color=INK); axs[0].set_ylim(-0.02, 1.05)
axs[0].set_title("Helix content vs J0  (dashed: exact no-steric, N=∞)", loc="left", fontsize=9.5, color=INK)
axs[1].set_yscale("log"); axs[1].set_ylabel("⟨S²⟩/N", color=INK); axs[1].set_title("Chiral order fluctuations", loc="left", fontsize=9.5, color=INK)
axs[2].set_ylabel("⟨Rₑ²⟩", color=INK); axs[2].set_title("Chain size", loc="left", fontsize=9.5, color=INK)
for ax in axs: ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("je_scan.png")
