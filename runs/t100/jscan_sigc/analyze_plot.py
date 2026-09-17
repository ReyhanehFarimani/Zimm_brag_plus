#!/usr/bin/env python3
"""jscan_sigc: does a longer chiral RANGE sigma_c let the chiral amplitude eps_s change anything?
Relative change [%] of helicity theta, handedness <|m|> and chain size <Rg^2> against eps_s, one line per J,
one column per sigma_c.  The reference is the eps_s = 0 cell of the same J (independent of sigma_c) and the
sigma_c = 0.70 a column is the jscan_eps scan, both read from ../jscan_eps.  FINISHED runs only, so this can be
run while the scan is going.  Errors: max(SEM over seeds, combined block error).
Vertical line = the eps_s at which a radially attractive configuration first appears for that sigma_c.
Style: no grid, ticks inside on all four sides, no title, legend on."""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
N0 = 200
LIMIT = {0.70: 8.2, 0.80: 9.2, 0.93: 10.4}                       # first attraction (clamped potential, theta0 = 100)
from matplotlib.colors import LinearSegmentedColormap
_CJ = LinearSegmentedColormap.from_list("J", ["#52c79e", "#22b07f", "#149166", "#0e7150", "#09523a", "#043325"])   # validated ramp, J = 1 .. 6
def cj(J): return _CJ((float(J) - 1.0) / 5.0)                      # J may be non-integer (4.5, 5.5)
INK, MUTED, SURF = "#0b0b0b", "#898781", "#ffffff"
KEYS = (("th", r"$\Delta\theta/\theta$  [%]", "helicity"), ("am", r"$\Delta\langle|m|\rangle/\langle|m|\rangle$  [%]", "handedness <|m|>"),
        ("rg2", r"$\Delta\langle R_g^2\rangle/\langle R_g^2\rangle$  [%]", "chain size <Rg2>"))

def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def blk(x, nb=20): b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return b.std(ddof=1) / np.sqrt(nb)
def load(folder, pattern, rx, sc_default=None):
    out = {}
    for log in sorted(glob.glob(os.path.join(folder, "logs", pattern))):
        if "summary" not in open(log).read(): continue
        g = re.search(rx, os.path.basename(log)).groups(); J, es = float(g[0]), float(g[1]); sc = float(g[2]) if sc_default is None else sc_default
        o = read_obs(os.path.join(folder, "out", os.path.basename(log)[:-4] + "_obs.dat")); m = np.abs(o["n_R"] - o["n_L"]) / N0
        out.setdefault((J, es, sc), []).append(dict(th=o["helicity"].mean(), am=m.mean(), rg2=o["Rg2"].mean(),
                                                    e_th=blk(o["helicity"]), e_am=blk(m), e_rg2=blk(o["Rg2"])))
    return out
runs = load("../jscan_eps", "J*_es*_s*.log", r"J([\d.]+)_es(\d+)_s(\d+)", sc_default=0.70)
new = load(".", "J*_es*_sc*_s*.log", r"J([\d.]+)_es([\d.]+)_sc([\d.]+)_s(\d+)"); runs.update(new)
def cell(key3, k):
    rs = runs[key3]; v = np.array([r[k] for r in rs]); sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    return v.mean(), max(sem, np.sqrt(sum(r["e_" + k]**2 for r in rs)) / len(rs)), len(rs)
Js = sorted({k[0] for k in runs}); SCs = sorted({k[2] for k in runs})
nfin = sum(len(v) for v in new.values()); ntot = len(glob.glob("inputs/J*.dat"))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(3, len(SCs), figsize=(4.2 * len(SCs), 10.5), squeeze=False); REL = {}
for r, (k, yl, name) in enumerate(KEYS):
    ymax = 0.0
    for c, sc in enumerate(SCs):
        x = ax[r, c]
        for n, J in enumerate(Js):
            if (J, 0.0, 0.70) not in runs: continue
            v0, e0, _ = cell((J, 0.0, 0.70), k); ess = sorted(e for (j, e, s) in runs if j == J and s == sc and e > 0)
            if not ess: continue
            V, Er, Ns = map(np.array, zip(*[cell((J, e, sc), k) for e in ess])); rel, rerr = 100 * (V - v0) / v0, 100 * np.hypot(Er, e0) / v0
            REL[(k, sc, J)] = (np.array(ess), rel, rerr, Ns); ymax = max(ymax, np.abs(rel).max() + rerr.max())
            x.errorbar(np.array(ess) + (n - 2.5) * 0.06, rel, yerr=rerr, fmt="o-", ms=4.5, lw=1.2, color=cj(J), mfc=cj(J), mec=SURF, mew=0.5, elinewidth=0.9, label=f"J = {J:g}")
        x.axhline(0, color=MUTED, lw=0.8); x.axvline(LIMIT[sc], color=MUTED, lw=0.8, ls=":")
        x.set_xlim(0, 13.5); x.set_xlabel(r"chiral amplitude $\epsilon_s$"); x.set_ylabel(yl + r"   relative to $\epsilon_s = 0$" if c == 0 else "")
        x.text(0.04, 0.95, rf"({'abcdefghi'[3 * r + c]}) {name.split(' ')[0]},  $\sigma_c$ = {sc:.2f} a" + ("  (default)" if sc == 0.70 else ""), transform=x.transAxes, va="top", fontsize=9, color=INK)
    for c in range(len(SCs)):
        ax[r, c].set_ylim(-1.5 * ymax, 1.5 * ymax)
        ax[r, c].text(LIMIT[SCs[c]] - 0.15, -1.42 * ymax, "first attraction", rotation=90, va="bottom", ha="right", fontsize=7.5, color=MUTED)
ax[0, 0].legend(loc="lower left", ncol=3, columnspacing=1.0, handlelength=1.5)
fig.tight_layout(); fig.savefig("range_scan.png", dpi=170); fig.savefig("range_scan.pdf")

print(f"range scan: {nfin} / {ntot} runs finished;  sigma_c columns: {SCs}")
for k, _, name in KEYS:
    print(f"\n== {name}: relative change vs eps_s = 0 in %  [z]   (n = seeds finished)")
    for sc in SCs:
        print(f"  sigma_c = {sc:.2f} a")
        zs = []
        for J in Js:
            if (k, sc, J) not in REL: continue
            ess, rel, rerr, Ns = REL[(k, sc, J)]; zs += list(rel / rerr)
            print(f"    J = {J:g}: " + "   ".join(f"eps {e:g}: {a:+6.3f}+-{b:5.3f} [{a / b:+4.1f}] n={n_}" for e, a, b, n_ in zip(ess, rel, rerr, Ns)))
        if zs: zs = np.array(zs); print(f"    -> RMS z = {np.sqrt((zs**2).mean()):.2f} (1 = noise), mean z = {zs.mean():+.2f}, largest |z| = {np.abs(zs).max():.1f}, cells = {len(zs)}")
    top = {sc: max(e for (kk, s, J), v in REL.items() if kk == k and s == sc for e in v[0]) for sc in SCs if any(kk == k and s == sc for (kk, s, J) in REL)}
    print("  at the largest eps_s of each sigma_c, all J combined sum(z)/sqrt(n):  " + "   ".join(
        f"sigma_c {sc:.2f} (eps {top[sc]:g}): {np.sum([v[1][list(v[0]).index(top[sc])] / v[2][list(v[0]).index(top[sc])] for (kk, s, J), v in REL.items() if kk == k and s == sc and top[sc] in v[0]]) / np.sqrt(sum(1 for (kk, s, J), v in REL.items() if kk == k and s == sc and top[sc] in v[0])):+.2f}" for sc in top))
