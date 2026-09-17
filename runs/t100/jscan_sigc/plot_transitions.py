#!/usr/bin/env python3
"""jscan_sigc transitions figure: the helix-coil and handedness crossovers at each chiral RANGE sigma_c.
  rows    = sigma_c: 0.70 a (default, the jscan_eps runs), 0.80 a, 0.93 a
  columns = helicity theta(J) | chi_theta = N var(theta) | handedness <|m|>(J) | Binder cumulant U4(J)
Blue curves   = EXACT 1D chain without non-bonded interactions, N = 100 .. 1600 (cache written by
                ../jscan_eps/plot_transitions.py).  A real transition would make the U4 curves cross at one J.
Orange marks  = MC at N = 200, ONE MARK PER eps_s (light -> dark with eps_s), side by side at each J because they
                coincide; eps_s = 0 does not depend on sigma_c and is the same set of runs in every row.
Hollow diamond = control without non-bonded interactions, at its true J (must sit on the exact N = 200 curve).
FINISHED runs only -- can be run while the scan is going.  Errors: max(SEM over seeds, combined block error) for
theta and <|m|>, SEM over seeds for chi_theta and U4.
Style: no grid, ticks inside on all four sides, no title, legend on."""
import glob, os, re, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
N0 = 200; NS = (100, 200, 400, 800, 1600); Jg = np.arange(0.5, 8.001, 0.125)
CACHE = "../jscan_eps/exact_cache.npz"
if not os.path.exists(CACHE): subprocess.run([sys.executable, "../jscan_eps/plot_transitions.py"], check=True, stdout=subprocess.DEVNULL)
z = np.load(CACHE); E = {N: {k: z[f"{N}_{k}"] for k in ("th", "chi_th", "am", "U4")} for N in NS}; th_inf = z["th_inf"]

CN = dict(zip(NS, ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]))                   # validated ordinal ramps
RAMP_E = ["#ef9c77", "#e9713f", "#d1501c", "#a53c13", "#752a0c", "#471806"]                  # eps_s, light -> dark
INK, MUTED, SURF = "#0b0b0b", "#898781", "#ffffff"

def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def blk(x, nb=20): b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return b.std(ddof=1) / np.sqrt(nb)
def load(folder, pattern, rx, fixed_sc=None):
    out = {}
    for log in sorted(glob.glob(os.path.join(folder, "logs", pattern))):
        if "summary" not in open(log).read(): continue
        g = re.search(rx, os.path.basename(log)).groups()
        key = (float(g[0]), "ctrl", None) if fixed_sc == "ctrl" else (float(g[0]), float(g[1]), float(g[2]) if fixed_sc is None else fixed_sc)
        o = read_obs(os.path.join(folder, "out", os.path.basename(log)[:-4] + "_obs.dat")); m = (o["n_R"] - o["n_L"]) / N0
        out.setdefault(key, []).append(dict(th=o["helicity"].mean(), am=np.abs(m).mean(), e_th=blk(o["helicity"]), e_am=blk(np.abs(m)), rg2=o["Rg2"].mean(), e_rg2=blk(o["Rg2"]),
                                            chi_th=N0 * o["helicity"].var(), U4=1 - (m**4).mean() / (3 * (m**2).mean()**2)))
    return out
runs = load("../jscan_eps", "J*_es*_s*.log", r"J([\d.]+)_es(\d+)_s(\d+)", fixed_sc=0.70)
runs.update(load("../jscan_eps", "J*_nonb_s*.log", r"J([\d.]+)_nonb", fixed_sc="ctrl"))
new = load(".", "J*_es*_sc*_s*.log", r"J([\d.]+)_es([\d.]+)_sc([\d.]+)_s(\d+)"); runs.update(new)
def cell(key3, k):
    rs = runs[key3]; v = np.array([r[k] for r in rs]); sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    eb = np.sqrt(sum(r["e_" + k]**2 for r in rs)) / len(rs) if "e_" + k in rs[0] else 0.0
    return v.mean(), max(sem, eb)
Js = sorted({k[0] for k in runs}); SCs = sorted({k[2] for k in runs if k[2] is not None})

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
COLS = (("th", r"helicity $\theta$"), ("chi_th", r"$\chi_\theta = N\,\mathrm{var}(\theta)$"), ("am", r"handedness $\langle |m| \rangle$"), ("U4", r"Binder cumulant $U_4$"))
fig, ax = plt.subplots(len(SCs), 4, figsize=(16.0, 4.1 * len(SCs)), squeeze=False)
for r, sc in enumerate(SCs):
    # eps_s = 0 (independent of sigma_c) + the eps_s values of this row, ranked light -> dark
    ess = [0.0] + sorted({k[1] for k in runs if k[2] == sc and k[1] != 0.0}); col = {e: RAMP_E[min(i, 5)] for i, e in enumerate(ess)}
    for c, (k, yl) in enumerate(COLS):
        x = ax[r, c]
        for N in NS: x.plot(Jg, E[N][k], color=CN[N], lw=1.8, label=f"exact 1D, N = {N}" if (r == 0 and c == 0) else None)
        if k == "th": x.plot(Jg, th_inf, color=INK, lw=1.1, ls="--", label=r"exact 1D, N $\to\infty$" if r == 0 else None)
        pts = [(J,) + cell((J, "ctrl", None), k) for J in Js if (J, "ctrl", None) in runs]
        if pts:
            X, V, Er = map(np.array, zip(*pts))
            x.errorbar(X, V, yerr=Er, fmt="D", ms=9, mfc="none", mec=INK, mew=1.2, ecolor=INK, elinewidth=0.9, zorder=4, label="control: no non-bonded" if k == "am" else None)
        for n, e in enumerate(ess):
            pts = [(J + (n - (len(ess) - 1) / 2) * 0.09,) + cell((J, e, 0.70 if e == 0.0 else sc), k) for J in Js if (J, e, 0.70 if e == 0.0 else sc) in runs]
            if not pts: continue
            X, V, Er = map(np.array, zip(*pts))
            x.errorbar(X, V, yerr=Er, fmt="o", ms=5, mfc=col[e], mec=SURF, mew=0.6, ecolor=col[e], elinewidth=0.9, zorder=5, label=rf"$\epsilon_s$ = {e:g}" if k == "am" else None)
        x.set_xlim(0.5, 8); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl)
        # tag where each panel is empty: helicity and chi_theta on the right, <|m|> and U4 top left
        tx, ty, ha = {"th": (0.96, 0.66, "right"), "chi_th": (0.96, 0.94, "right"), "am": (0.04, 0.94, "left"), "U4": (0.04, 0.84, "left")}[k]
        x.text(tx, ty, rf"({'abcdefghijkl'[4 * r + c]})  $\sigma_c$ = {sc:.2f} a" + ("  (default)" if sc == 0.70 else ""), transform=x.transAxes, va="top", ha=ha, fontsize=9, color=INK)
        if k == "U4": x.axhline(2 / 3, color=MUTED, lw=0.8, ls=":"); x.text(0.7, 0.655, "2/3 = fully ordered", color=MUTED, fontsize=8, va="top")
        if k == "am": x.legend(loc="upper left", bbox_to_anchor=(0.0, 0.89), handletextpad=0.3, labelspacing=0.35, title="MC, N = 200", title_fontsize=8)
ax[0, 0].legend(loc="lower right")
fig.tight_layout(); fig.savefig("transitions_range.png", dpi=150); fig.savefig("transitions_range.pdf")

# ---------------------------------------------------------------- second figure: NO exact curves -- one line per eps_s across J
COLS5 = COLS + (("rg2", r"chain size $\langle R_g^2 \rangle$ [$a^2$]"),)
fig, ax = plt.subplots(len(SCs), 5, figsize=(19.5, 4.1 * len(SCs)), squeeze=False)
for r, sc in enumerate(SCs):
    ess = [0.0] + sorted({k[1] for k in runs if k[2] == sc and k[1] != 0.0}); col = {e: RAMP_E[min(i, 5)] for i, e in enumerate(ess)}
    for c, (k, yl) in enumerate(COLS5):
        x = ax[r, c]
        for e in ess:
            pts = [(J,) + cell((J, e, 0.70 if e == 0.0 else sc), k) for J in Js if (J, e, 0.70 if e == 0.0 else sc) in runs]
            if not pts: continue
            X, V, Er = map(np.array, zip(*pts))
            x.errorbar(X, V, yerr=Er, fmt="o-", ms=4.5, lw=1.4, color=col[e], mfc=col[e], mec=SURF, mew=0.5, elinewidth=0.9, zorder=3 + ess.index(e), label=rf"$\epsilon_s$ = {e:g}")
        x.set_xlim(0.5, 6.5); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl)
        tx, ty, ha = {"th": (0.96, 0.10, "right"), "chi_th": (0.96, 0.94, "right"), "am": (0.04, 0.94, "left"), "U4": (0.04, 0.94, "left"), "rg2": (0.04, 0.94, "left")}[k]
        x.text(tx, ty, rf"({'abcdefghijklmno'[5 * r + c]})  $\theta_0$ = 100,  $\sigma_c$ = {sc:.2f} a", transform=x.transAxes, va="top" if ty > 0.5 else "bottom", ha=ha, fontsize=9, color=INK)
        if k == "am": x.legend(loc="upper left", bbox_to_anchor=(0.0, 0.89), handletextpad=0.4, labelspacing=0.35, title="MC, N = 200", title_fontsize=8)
fig.tight_layout(); fig.savefig("eps_lines_range.png", dpi=140); fig.savefig("eps_lines_range.pdf"); plt.close(fig)

nfin = sum(len(v) for v in new.values())
print(f"transitions_range: {nfin} / {len(glob.glob('inputs/J*.dat'))} range-scan runs finished;  rows: sigma_c = {SCs}")
print("  does eps_s or sigma_c move the crossovers?  value at J = 2 (helix-coil midpoint) and J = 4 (handedness crossover, N = 200)")
for sc in SCs:
    for J, k, nm in ((2, "th", "theta(J=2)"), (4, "am", "<|m|>(J=4)"), (4, "U4", "U4(J=4)")):
        ess = [0.0] + sorted({kk[1] for kk in runs if kk[2] == sc and kk[1] != 0.0 and kk[0] == J})
        print(f"    sigma_c {sc:.2f}  {nm:11s} " + "  ".join(f"eps {e:g}: {cell((J, e, 0.70 if e == 0.0 else sc), k)[0]:.4f}+-{cell((J, e, 0.70 if e == 0.0 else sc), k)[1]:.4f}"
                                                                for e in ess if (J, e, 0.70 if e == 0.0 else sc) in runs)
              + f"   | exact 1D: {E[N0][k][int(round((J - Jg[0]) / 0.125))]:.4f}")
