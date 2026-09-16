#!/usr/bin/env python3
"""t100 parameter set: observables vs the helix ground energy E_helix (with sterics), compared with the exact
no-steric transfer-matrix solution of the same model (E_HH = -J0, E_CH = +J1, E_RL = +J2 from the input)."""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here)
sys.argv = [sys.argv[0]]
ns = open("../no_steric/analyze.py").read()
exec(ns.split("inp0 = read_input")[0])                       # readers (+ deg-aware read_input)
exec(ns[ns.index("def triple_class"):ns.index("return dict(hel=hel, Es=Es, Eb=Eb, Re2=Re2, S2N=tot / N)") + len("return dict(hel=hel, Es=Es, Eb=Eb, Re2=Re2, S2N=tot / N)")])

inp0 = read_input("sample_data.dat"); N = int(inp0["N"]); J0, J1, J2 = (float(inp0[k]) for k in ("J0", "J1", "J2"))
runs = {}
for f in sorted(glob.glob("out/E*_obs.dat")):
    E = float(re.search(r"E(-?[\d.]+)_s(\d+)_obs", f).group(1))
    o = read_obs(f); n = len(o["helicity"])
    if n < 40: continue
    nb = 20; blk = lambda x: x[: n // nb * nb].reshape(nb, -1).mean(1)
    S = o["n_R"] - o["n_L"]
    obs = dict(hel=o["helicity"], Es=o["E_state"], Eb=o["E_bend"], Enb=o["E_nb"], Re2=o["Ree"] ** 2, Rg2=o["Rg2"], S2N=S ** 2 / N)
    runs.setdefault(E, []).append({k: (blk(v).mean(), blk(v).std(ddof=1) / np.sqrt(nb)) for k, v in obs.items()})
Es = sorted(runs); keys = ("hel", "Es", "Eb", "Enb", "Re2", "Rg2", "S2N")
mc  = {k: np.array([np.mean([r[k][0] for r in runs[E]]) for E in Es]) for k in keys}
err = {k: np.array([np.sqrt(np.sum([r[k][1] ** 2 for r in runs[E]])) / len(runs[E]) for E in Es]) for k in keys}
def ex(E):
    d = dict(inp0); d["E_helix"] = str(E); return exact_bend(d, J0, J2=J2)
Eline = np.linspace(min(Es) - 0.5, max(Es) + 0.5, 60); exl = [ex(E) for E in Eline]

print(f"t100: N={N}, E_HH={-J0:+g}, E_CH={J1:+g}, E_RL={J2:+g}, sterics on (coil core A={inp0['nb_A']}, rods L={inp0['rod_L']} r={inp0['rod_r']}); {len(runs[Es[0]])} seeds")
print(f"{'E_helix':>8} | {'helicity':>8} {'+-':>6} {'no-ster':>8} | {'<S2>/N':>8} {'+-':>6} {'no-ster':>8} | {'<Re2>':>7} {'+-':>5} {'no-ster':>7} | {'<Rg2>':>7} | {'E_nb':>6} {'E_bend':>7}")
for i, E in enumerate(Es):
    e = ex(E)
    print(f"{E:>8} | {mc['hel'][i]:8.4f} {err['hel'][i]:6.4f} {e['hel']:8.4f} | {mc['S2N'][i]:8.3f} {err['S2N'][i]:6.3f} {e['S2N']:8.3f} | {mc['Re2'][i]:7.0f} {err['Re2'][i]:5.0f} {e['Re2']:7.0f} | {mc['Rg2'][i]:7.1f} | {mc['Enb'][i]:6.2f} {mc['Eb'][i]:7.1f}")

BLUE, ORANGE, GREEN, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.1), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.set_xlabel("E_helix  (ground energy per helical residue)", color=INK)
for ax, key, col, ylab in ((axs[0], "hel", BLUE, "helicity"), (axs[1], "Re2", ORANGE, r"$\langle R_e^2\rangle$"), (axs[2], "S2N", GREEN, r"$\langle S^2\rangle/N$")):
    ax.plot(Eline, [e[key] for e in exl], "-", lw=1.5, color=INK2, label="no sterics (exact)")
    ax.errorbar(Es, mc[key], yerr=err[key], fmt="o", ms=5, color=col, capsize=2.5, label="with sterics (MC)")
    ax.set_ylabel(ylab, color=INK); ax.legend(frameon=False, fontsize=8)
axs[0].set_ylim(0, 1.02)
axs[0].set_title(f"t100: helix content  (E_HH={-J0:+g}, E_CH={J1:+g}, E_RL={J2:+g})", loc="left", fontsize=9.5, color=INK)
axs[1].set_title("Chain size", loc="left", fontsize=9.5, color=INK); axs[2].set_title("Chiral order fluctuations", loc="left", fontsize=9.5, color=INK)
fig.tight_layout(); fig.savefig("t100_scan.png")
