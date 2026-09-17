#!/usr/bin/env python3
"""t45, ZOOM on the transition region J = 2.5 .. 5 (step 0.25): does a same-handed attraction change anything there?
Potential exactly as fitted, ONE eps_s for every pair type, chiral range 0.87 a:
    eps_s = 0     reference (no chiral term)                                   [../jscan_sigc, J*_es0_sc0_s*]
    eps_s = 7.6   same-handed well about -0.5 kT, R.L purely repulsive
    eps_s = 10.5  same-handed well about -2.7 kT, R.L purely repulsive
    eps_s = 13    STRONGEST: same-handed well about -4.6 kT, R.L well about -0.7 kT
Rows 1-2: helicity, its susceptibility, chain size / handedness <|m|>, Binder cumulant, non-bonded energy, with the exact
1D curve (no non-bonded interactions, N = 200) where it exists.  Row 3: the DIFFERENCE to eps_s = 0 at the same J with
its error -- the sensitive view.  FINISHED runs only; errors = max(SEM over seeds, combined block error).
Style: no grid, ticks inside on all four sides, no title, legend on (outside the panels)."""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
N0 = 200; JLO, JHI = 2.5, 5.0
EPS = [0.0, 7.6, 10.5, 13.0]
COL = {0.0: "#eb6834", 7.6: "#8f82e0", 10.5: "#4c3cae", 13.0: "#2f2378"}        # orange = reference; validated violet ramp = attraction
MRK = {0.0: "o", 7.6: "^", 10.5: "s", 13.0: "D"}
LAB = {0.0: r"$\epsilon_s$ = 0 (reference)", 7.6: r"$\epsilon_s$ = 7.6  (well $-0.5\,k_BT$)", 10.5: r"$\epsilon_s$ = 10.5  (well $-2.7\,k_BT$)",
       13.0: r"$\epsilon_s$ = 13  (well $-4.6\,k_BT$; R–L $-0.7$)"}
INK, MUTED, SURF, EXACT = "#0b0b0b", "#898781", "#ffffff", "#2a78d6"

def read_obs(f):
    L = open(f).read().splitlines(); h = next(l for l in L if l.startswith("#")).split()[1:]
    a = np.array([[float(x) for x in l.split()] for l in L if not l.startswith("#") and len(l.split()) == len(h)]); return {k: a[:, i] for i, k in enumerate(h)}
def blk(x, nb=20): b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return b.std(ddof=1) / np.sqrt(nb)
runs = {}
for d, pat in (("../jscan_sigc", "J*_es0_sc0_s*.log"), (".", "J*_sc0.87_s*.log")):
    for log in glob.glob(f"{d}/logs/{pat}"):
        if "summary" not in open(log).read(): continue
        J, es = map(float, re.search(r"J([\d.]+)_es([\d.]+)_", os.path.basename(log)).groups())
        if es not in EPS or not (JLO - 1e-9 <= J <= JHI + 1e-9): continue
        o = read_obs(f"{d}/out/" + os.path.basename(log)[:-4] + "_obs.dat"); m = (o["n_R"] - o["n_L"]) / N0; am = np.abs(m)
        runs.setdefault((J, es), []).append(dict(th=o["helicity"].mean(), e_th=blk(o["helicity"]), am=am.mean(), e_am=blk(am), rg2=o["Rg2"].mean(), e_rg2=blk(o["Rg2"]),
            enb=o["E_nb"].mean(), e_enb=blk(o["E_nb"]), chi=N0 * o["helicity"].var(), U4=1 - (m**4).mean() / (3 * (m**2).mean()**2)))
def cell(key, k):
    rs = runs[key]; v = np.array([r[k] for r in rs]); n = len(v); sem = v.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0
    eb = np.sqrt(sum(r["e_" + k]**2 for r in rs)) / n if "e_" + k in rs[0] else 0.0
    return v.mean(), max(sem, eb), n
Js = sorted({k[0] for k in runs})
z = np.load("../jscan_sigc/exact_cache.npz"); Jg = np.arange(0.5, 10.001, 0.125); sel = (Jg >= JLO - 1e-9) & (Jg <= JHI + 1e-9)
EX = {"th": z["200_th"], "chi": z["200_chi_th"], "am": z["200_am"], "U4": z["200_U4"]}

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(3, 3, figsize=(15.5, 11.5))
PAN = [("th", r"helicity $\theta$"), ("chi", r"$\chi_\theta = N\,\mathrm{var}(\theta)$"), ("rg2", r"chain size $\langle R_g^2\rangle$ [$a^2$]"),
       ("am", r"handedness $\langle|m|\rangle$"), ("U4", r"Binder cumulant $U_4$"), ("enb", r"non-bonded energy $\langle E_{nb}\rangle$ [$k_BT$]")]
for i, (k, yl) in enumerate(PAN):
    x = ax[i // 3, i % 3]
    if k in EX: x.plot(Jg[sel], EX[k][sel], color=EXACT, lw=1.6, label="exact 1D, no non-bonded, N = 200")
    for n, e in enumerate(EPS):
        pts = [(J + (n - 1.5) * 0.03,) + cell((J, e), k)[:2] for J in Js if (J, e) in runs]
        if not pts: continue
        X, V, Er = map(np.array, zip(*pts)); x.errorbar(X, V, yerr=Er, fmt=MRK[e] + "-", ms=5, lw=1.2, color=COL[e], mfc=COL[e], mec=SURF, mew=0.5, elinewidth=0.9, label=LAB[e], zorder=3 + n)
    x.set_xlim(JLO - 0.15, JHI + 0.15); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl)
    x.text(0.04, 0.94, f"({'abcdef'[i]})", transform=x.transAxes, va="top", fontsize=9, color=INK)
DIF = [("th", r"$\theta(\epsilon_s) - \theta(0)$"), ("am", r"$\langle|m|\rangle(\epsilon_s) - \langle|m|\rangle(0)$"), ("rg2", r"$\Delta\langle R_g^2\rangle / \langle R_g^2\rangle$  [%]")]
ROWS = {}
for i, (k, yl) in enumerate(DIF):
    x = ax[2, i]; x.axhline(0, color=MUTED, lw=0.8)
    for n, e in enumerate(EPS[1:]):
        pts = []
        for J in Js:
            if (J, e) in runs and (J, 0.0) in runs:
                v, er, nn = cell((J, e), k); v0, e0, n0 = cell((J, 0.0), k); sc = 100 / v0 if k == "rg2" else 1.0
                pts.append((J + (n - 1) * 0.03, (v - v0) * sc, np.hypot(er, e0) * sc)); ROWS.setdefault((k, e), []).append((J, (v - v0) * sc, np.hypot(er, e0) * sc, nn))
        if not pts: continue
        X, V, Er = map(np.array, zip(*pts)); x.errorbar(X, V, yerr=Er, fmt=MRK[e] + "-", ms=5, lw=1.2, color=COL[e], mfc=COL[e], mec=SURF, mew=0.5, elinewidth=0.9, zorder=3 + n)
    x.set_xlim(JLO - 0.15, JHI + 0.15); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl + r"   (difference to $\epsilon_s = 0$)")
    x.text(0.04, 0.94, f"({'ghi'[i]})", transform=x.transAxes, va="top", fontsize=9, color=INK)
    yl_ = 1.3 * max(abs(np.array(x.get_ylim()))); x.set_ylim(-yl_, yl_)
h, l = ax[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="center left", bbox_to_anchor=(0.835, 0.5), title=r"$\theta_0$ = 45, range 0.87 a, N = 200", title_fontsize=8)
fig.tight_layout(rect=(0, 0, 0.83, 1)); fig.savefig("zoom_transition.png", dpi=150); fig.savefig("zoom_transition.pdf")

print(f"zoom J = {JLO} .. {JHI}: finished runs per eps_s: " + ", ".join(f"{e:g}: {sum(len(v) for kk, v in runs.items() if kk[1] == e)}" for e in EPS))
for k, nm in (("th", "helicity"), ("am", "<|m|>"), ("rg2", "Rg2 [%]")):
    for e in EPS[1:]:
        if (k, e) not in ROWS: continue
        r = ROWS[(k, e)]; zz = np.array([d / er for (_, d, er, _) in r])
        print(f"  {nm:9s} eps_s = {e:<5g}: " + "  ".join(f"J{J:g}: {d:+.4f}[{d / er:+.1f}]n{nn}" for (J, d, er, nn) in r) + f"   | RMS z {np.sqrt((zz**2).mean()):.2f}, sum z/sqrt(n) {zz.sum() / np.sqrt(len(zz)):+.2f}")
