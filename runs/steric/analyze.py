#!/usr/bin/env python3
"""J-scan with sterics (Gaussian coil core + Gay-Berne helix rods) vs the same model without sterics
(exact transfer-matrix solution and the no-steric MC of runs/no_steric)."""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here)
sys.argv = [sys.argv[0]]
ns = open("../no_steric/analyze.py").read()
exec(ns.split("inp0 = read_input")[0])                                   # readers + exact()
exec(open("../no_steric/analyze.py").read()[1283:5260])                    # triple_class, hinge_moments, exact_bend
def collect(pattern, inputs):
    runs = {}
    for f in sorted(glob.glob(pattern)):
        J = float(re.search(r"J([\d.]+)_s(\d+)_obs", f).group(1))
        o = read_obs(f); n = len(o["helicity"])
        if n < 40: continue
        nb = 20; blk = lambda x: x[: n // nb * nb].reshape(nb, -1).mean(1)
        S = o["n_R"] - o["n_L"]; N = int(read_input(inputs)["N"])
        obs = dict(hel=o["helicity"], Es=o["E_state"], Eb=o["E_bend"], Enb=o.get("E_nb", np.zeros(n)), Re2=o["Ree"] ** 2, Rg2=o["Rg2"], S2N=S ** 2 / N)
        runs.setdefault(J, []).append({k: (blk(v).mean(), blk(v).std(ddof=1) / np.sqrt(nb)) for k, v in obs.items()})
    Js = sorted(runs); keys = ("hel", "Es", "Eb", "Enb", "Re2", "Rg2", "S2N")
    mc  = {k: np.array([np.mean([r[k][0] for r in runs[J]]) for J in Js]) for k in keys}
    err = {k: np.array([np.sqrt(np.sum([r[k][1] ** 2 for r in runs[J]])) / len(runs[J]) for J in Js]) for k in keys}
    return Js, mc, err, len(runs[Js[0]]) if Js else 0

inpg = read_input(sorted(glob.glob("inputs/J*_s1.dat"))[0]); N = int(inpg["N"])
Js, mc, err, nseed = collect("out/J*_obs.dat", sorted(glob.glob("inputs/J*_s1.dat"))[0])
Jn, mcn, errn, _ = collect("../no_steric/out/J*_obs.dat", sorted(glob.glob("../no_steric/inputs/J*_s1.dat"))[0])
Jline = np.linspace(min(Js), max(Js), 40); ex = [exact_bend(inpg, J) for J in Jline]; exJ = {J: exact_bend(inpg, J) for J in Js}

print(f"N = {N}, J1 = {inpg['J1']}, J0 = J2 = J, theta0 = 180 deg, sterics: coil core A={inpg['nb_A']} sigma={inpg['nb_sigma']}, helix rods L={inpg['rod_L']} D={2*float(inpg['rod_r'])}; {nseed} seeds")
print(f"{'J':>5} | {'helicity':>8} {'+-':>6} {'no-ster':>8} | {'<Re2>':>7} {'+-':>5} {'no-ster':>7} | {'<Rg2>':>7} {'+-':>5} | {'<S2>/N':>8} {'+-':>6} {'no-ster':>8} | {'E_nb':>6} {'E_bend':>7}")
for i, J in enumerate(Js):
    e = exJ[J]
    print(f"{J:>5} | {mc['hel'][i]:8.4f} {err['hel'][i]:6.4f} {e['hel']:8.4f} | {mc['Re2'][i]:7.0f} {err['Re2'][i]:5.0f} {e['Re2']:7.0f} | {mc['Rg2'][i]:7.1f} {err['Rg2'][i]:5.1f} | {mc['S2N'][i]:8.3f} {err['S2N'][i]:6.3f} {e['S2N']:8.3f} | {mc['Enb'][i]:6.2f} {mc['Eb'][i]:7.1f}")
print("(no-ster = exact transfer-matrix solution of the same model without sterics)")

BLUE, ORANGE, GREEN, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.1), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.set_xlabel("J   (J0 = J2 = J)", color=INK)
for ax, key, col, lab, ylab in ((axs[0], "hel", BLUE, "helicity", "helicity"), (axs[1], "Re2", ORANGE, "⟨Rₑ²⟩", r"$\langle R_e^2\rangle$"), (axs[2], "S2N", GREEN, "⟨S²⟩/N", r"$\langle S^2\rangle/N$")):
    ax.plot(Jline, [e[key] for e in ex], "-", lw=1.5, color=INK2, label="no sterics (exact)")
    if Jn: ax.errorbar(Jn, mcn[key], yerr=errn[key], fmt="s", ms=4, mfc=BG, color=INK2, capsize=2, label="no sterics (MC)")
    ax.errorbar(Js, mc[key], yerr=err[key], fmt="o", ms=5, color=col, capsize=2.5, label="with sterics (MC)")
    ax.set_ylabel(ylab, color=INK); ax.legend(frameon=False, fontsize=8)
axs[0].set_ylim(0, 1.02); axs[2].set_yscale("log")
axs[0].set_title(f"Helix content  (N={N}, θ0=π, coil core + GB rods)", loc="left", fontsize=9.5, color=INK)
axs[1].set_title("Chain size", loc="left", fontsize=9.5, color=INK); axs[2].set_title("Chiral order fluctuations", loc="left", fontsize=9.5, color=INK)
fig.tight_layout(); fig.savefig("jscan_steric.png")
