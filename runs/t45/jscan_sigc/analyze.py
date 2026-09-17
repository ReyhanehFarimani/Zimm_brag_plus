#!/usr/bin/env python3
"""theta0 = 45 scan: transitions and the effect of the chiral amplitude eps_s at three chiral ranges sigma_c.
  transitions_t45.png : rows = sigma_c (0.51, 0.62, 0.87 a); columns = helicity theta(J), chi_theta = N var(theta),
                        handedness <|m|>(J), Binder cumulant U4(J).  Blue curves = EXACT 1D chain without non-bonded
                        interactions for the t45 parameters (../sample_data.dat), N = 100 .. 1600; orange marks =
                        MC N = 200, one per eps_s (light -> dark), side by side; hollow diamond = control without
                        non-bonded interactions (runs J*_nonb_s*, if present) at its true J.
  eps_t45.png         : relative change [%] of theta, <|m|>, <Rg^2> against eps_s, one line per J, one column per
                        sigma_c; reference = the eps_s = 0 runs (J*_es0_sc0_s*), which do not depend on sigma_c.
FINISHED runs only.  Errors: max(SEM over seeds, combined block error); SEM over seeds for chi_theta and U4.
Style: no grid, ticks inside on all four sides, no title, legend on."""
import contextlib, glob, importlib.util, io, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
N0 = 200; NS = (100, 200, 400, 800, 1600); Jg = np.arange(0.5, 10.001, 0.125)
LIMIT = {0.51: 3.1, 0.62: 4.5, 0.87: 6.6}                         # first attraction (clamped potential, theta0 = 45)
CN = dict(zip(NS, ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]))                   # validated ordinal ramps
RAMP_E = ["#ef9c77", "#e9713f", "#d1501c", "#a53c13", "#752a0c", "#471806"]
CJ = dict(zip([1, 2, 3, 4, 5, 6], ["#52c79e", "#22b07f", "#149166", "#0e7150", "#09523a", "#043325"]))
INK, MUTED, SURF = "#0b0b0b", "#898781", "#ffffff"

# ---------------------------------------------------------------- exact 1D reference for the t45 parameters
sys.argv = ["x", "../sample_data.dat"]
spec = importlib.util.spec_from_file_location("ex", os.path.abspath("../../t100/exact_fss.py")); ex = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(ex)
PAIRS = sorted(ex.idx, key=ex.idx.get); SECOND = np.array([b for (a, b) in PAIRS])
def exact_dist(J, N, helix=False):
    ex.J0, ex.J1, ex.J2 = float(J), 0.0, float(J); T, v0 = ex.build(0.0)
    f = (lambda s: abs(s)) if helix else (lambda s: s); off = 0 if helix else N
    P = np.zeros((9, (N + 1) if helix else (2 * N + 1)))
    for p, (a, b) in enumerate(PAIRS): P[p, off + f(a) + f(b)] = v0[p]
    for _ in range(N - 2):
        Q = T.T @ P
        for q in range(9): Q[q] = np.roll(Q[q], f(SECOND[q]))
        P = Q / Q.sum()
    p = P.sum(0); return p / p.sum()
E = {N: {k: np.zeros(len(Jg)) for k in ("th", "chi_th", "am", "U4")} for N in NS}; CACHE = "exact_cache.npz"
if os.path.exists(CACHE) and os.path.getmtime(CACHE) > os.path.getmtime("../sample_data.dat"):
    z = np.load(CACHE)
    for N in NS:
        for k in E[N]: E[N][k] = z[f"{N}_{k}"]
else:
    for N in NS:
        m = np.arange(-N, N + 1) / N; h = np.arange(N + 1) / N
        for i, J in enumerate(Jg):
            pm = exact_dist(J, N); ph = exact_dist(J, N, helix=True); m2, m4 = (pm * m**2).sum(), (pm * m**4).sum(); th = (ph * h).sum()
            E[N]["am"][i] = (pm * np.abs(m)).sum(); E[N]["U4"][i] = 1 - m4 / (3 * m2 * m2); E[N]["th"][i] = th; E[N]["chi_th"][i] = N * ((ph * h * h).sum() - th * th)
    np.savez(CACHE, **{f"{N}_{k}": E[N][k] for N in NS for k in E[N]})

# ---------------------------------------------------------------- MC
def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def blk(x, nb=20): b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return b.std(ddof=1) / np.sqrt(nb)
runs = {}
for log in sorted(glob.glob("logs/J*.log")):
    if "summary" not in open(log).read(): continue
    b = os.path.basename(log); g = re.search(r"J(\d+)_es([\d.]+)_sc([\d.]+)_s(\d+)", b)
    key = (int(g.group(1)), float(g.group(2)), float(g.group(3))) if g else (int(re.search(r"J(\d+)_nonb", b).group(1)), "ctrl", None)
    if g and key[1] == 0.0: key = (key[0], 0.0, None)                                # eps_s = 0: no range
    o = read_obs("out/" + b[:-4] + "_obs.dat"); m = (o["n_R"] - o["n_L"]) / N0
    runs.setdefault(key, []).append(dict(th=o["helicity"].mean(), am=np.abs(m).mean(), rg2=o["Rg2"].mean(), e_th=blk(o["helicity"]), e_am=blk(np.abs(m)),
        e_rg2=blk(o["Rg2"]), chi_th=N0 * o["helicity"].var(), U4=1 - (m**4).mean() / (3 * (m**2).mean()**2), enb=o["E_nb"].mean()))
def cell(key3, k):
    rs = runs[key3]; v = np.array([r[k] for r in rs]); sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    eb = np.sqrt(sum(r["e_" + k]**2 for r in rs)) / len(rs) if "e_" + k in rs[0] else 0.0
    return v.mean(), max(sem, eb), len(rs)
Js = sorted({k[0] for k in runs}); SCs = sorted({k[2] for k in runs if k[2] is not None}) or [0.51]
kk = lambda J, e, sc: (J, e, None if e == 0.0 else sc)
nfin = sum(len(v) for v in runs.values()); ntot = len(glob.glob("inputs/J*.dat"))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
# ---------------------------------------------------------------- figure 1: transitions
COLS = (("th", r"helicity $\theta$"), ("chi_th", r"$\chi_\theta = N\,\mathrm{var}(\theta)$"), ("am", r"handedness $\langle |m| \rangle$"), ("U4", r"Binder cumulant $U_4$"))
fig, ax = plt.subplots(len(SCs), 4, figsize=(16.0, 4.1 * len(SCs)), squeeze=False)
for r, sc in enumerate(SCs):
    ess = [0.0] + sorted({k[1] for k in runs if k[2] == sc}); col = {e: RAMP_E[min(i, 5)] for i, e in enumerate(ess)}
    for c, (k, yl) in enumerate(COLS):
        x = ax[r, c]
        for N in NS: x.plot(Jg, E[N][k], color=CN[N], lw=1.8, label=f"exact 1D, N = {N}" if (r == 0 and c == 0) else None)
        pts = [(J,) + cell((J, "ctrl", None), k)[:2] for J in Js if (J, "ctrl", None) in runs]
        if pts:
            X, V, Er = map(np.array, zip(*pts)); x.errorbar(X, V, yerr=Er, fmt="D", ms=9, mfc="none", mec=INK, mew=1.2, ecolor=INK, elinewidth=0.9, zorder=4, label="control: no non-bonded" if k == "am" else None)
        for n, e in enumerate(ess):
            pts = [(J + (n - (len(ess) - 1) / 2) * 0.09,) + cell(kk(J, e, sc), k)[:2] for J in Js if kk(J, e, sc) in runs]
            if not pts: continue
            X, V, Er = map(np.array, zip(*pts)); x.errorbar(X, V, yerr=Er, fmt="o", ms=5, mfc=col[e], mec=SURF, mew=0.6, ecolor=col[e], elinewidth=0.9, zorder=5, label=rf"$\epsilon_s$ = {e:g}" if k == "am" else None)
        x.set_xlim(0.5, 10); x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl)
        tx, ty, ha = {"th": (0.96, 0.66, "right"), "chi_th": (0.96, 0.94, "right"), "am": (0.04, 0.94, "left"), "U4": (0.04, 0.84, "left")}[k]
        x.text(tx, ty, rf"({'abcdefghijkl'[4 * r + c]})  $\theta_0$ = 45,  $\sigma_c$ = {sc:.2f} a", transform=x.transAxes, va="top", ha=ha, fontsize=9, color=INK)
        if k == "U4": x.axhline(2 / 3, color=MUTED, lw=0.8, ls=":"); x.text(0.7, 0.655, "2/3 = fully ordered", color=MUTED, fontsize=8, va="top")
        if k == "am": x.legend(loc="upper left", bbox_to_anchor=(0.0, 0.89), handletextpad=0.3, labelspacing=0.35, title="MC, N = 200", title_fontsize=8)
ax[0, 0].legend(loc="lower right"); fig.tight_layout(); fig.savefig("transitions_t45.png", dpi=150); fig.savefig("transitions_t45.pdf"); plt.close(fig)

# ---------------------------------------------------------------- figure 2: relative change vs eps_s
KEYS = (("th", r"$\Delta\theta/\theta$  [%]", "helicity"), ("am", r"$\Delta\langle|m|\rangle/\langle|m|\rangle$  [%]", "handedness <|m|>"), ("rg2", r"$\Delta\langle R_g^2\rangle/\langle R_g^2\rangle$  [%]", "chain size <Rg2>"))
fig, ax = plt.subplots(3, len(SCs), figsize=(4.2 * len(SCs), 10.5), squeeze=False); REL = {}
for r, (k, yl, name) in enumerate(KEYS):
    ymax = 1e-9
    for c, sc in enumerate(SCs):
        x = ax[r, c]
        for n, J in enumerate(Js):
            if (J, 0.0, None) not in runs: continue
            v0, e0, _ = cell((J, 0.0, None), k); ess = sorted(e for (j, e, s) in runs if j == J and s == sc)
            if not ess: continue
            V, Er, Ns = map(np.array, zip(*[cell((J, e, sc), k) for e in ess])); rel, rerr = 100 * (V - v0) / v0, 100 * np.hypot(Er, e0) / v0
            REL[(k, sc, J)] = (np.array(ess), rel, rerr, Ns); ymax = max(ymax, np.abs(rel).max() + rerr.max())
            x.errorbar(np.array(ess) + (n - 2.5) * 0.04, rel, yerr=rerr, fmt="o-", ms=4.5, lw=1.2, color=CJ[J], mfc=CJ[J], mec=SURF, mew=0.5, elinewidth=0.9, label=f"J = {J}")
        x.axhline(0, color=MUTED, lw=0.8); x.axvline(LIMIT.get(sc, np.nan), color=MUTED, lw=0.8, ls=":"); x.axvline(1.0, color=MUTED, lw=0.8, ls="--")
        x.set_xlim(0, 8.5); x.set_xlabel(r"chiral amplitude $\epsilon_s$"); x.set_ylabel(yl + r"   relative to $\epsilon_s = 0$" if c == 0 else "")
        x.text(0.04, 0.95, rf"({'abcdefghi'[3 * r + c]}) {name.split(' ')[0]},  $\sigma_c$ = {sc:.2f} a", transform=x.transAxes, va="top", fontsize=9, color=INK)
    for c in range(len(SCs)):
        ax[r, c].set_ylim(-1.5 * ymax, 1.5 * ymax)
        ax[r, c].text(LIMIT.get(SCs[c], 0) - 0.1, -1.42 * ymax, "first attraction", rotation=90, va="bottom", ha="right", fontsize=7.5, color=MUTED)
        ax[r, c].text(1.0 - 0.1, -1.42 * ymax, "as measured", rotation=90, va="bottom", ha="right", fontsize=7.5, color=MUTED)
ax[0, 0].legend(loc="lower right", ncol=2, columnspacing=1.0, handlelength=1.5)
fig.tight_layout(); fig.savefig("eps_t45.png", dpi=170); fig.savefig("eps_t45.pdf"); plt.close(fig)

# ---------------------------------------------------------------- numbers
print(f"t45 scan: {nfin} / {ntot} runs finished;  sigma_c = {SCs};  eps_s = 0 reference cells: {sum(1 for k in runs if k[1] == 0.0)}")
for k, _, name in KEYS:
    print(f"\n== {name}: relative change vs eps_s = 0 in %  [z]  n = seeds")
    for sc in SCs:
        zs = []
        for J in Js:
            if (k, sc, J) not in REL: continue
            ess, rel, rerr, Ns = REL[(k, sc, J)]; zs += list(rel / rerr)
            print(f"  sigma_c {sc:.2f}  J = {J}: " + "   ".join(f"eps {e:g}: {a:+6.3f}+-{b:5.3f} [{a / b:+4.1f}] n={n_}" for e, a, b, n_ in zip(ess, rel, rerr, Ns)))
        if zs: zs = np.array(zs); print(f"  sigma_c {sc:.2f}  -> RMS z = {np.sqrt((zs**2).mean()):.2f} (1 = noise), mean z = {zs.mean():+.2f}, largest |z| = {np.abs(zs).max():.1f}, cells = {len(zs)}")
print("\n== MC vs exact 1D (N = 200) at eps_s = 0:  J | theta MC / exact | <|m|> MC / exact | <E_nb> [kT]")
for J in Js:
    if (J, 0.0, None) not in runs: continue
    i = int(round((J - Jg[0]) / 0.125)); t, a = cell((J, 0.0, None), "th"), cell((J, 0.0, None), "am")
    print(f"   J = {J}:  {t[0]:.4f}+-{t[1]:.4f} / {E[N0]['th'][i]:.4f}    {a[0]:.4f}+-{a[1]:.4f} / {E[N0]['am'][i]:.4f}    {cell((J, 0.0, None), 'enb')[0]:.3f}")
