#!/usr/bin/env python3
"""jscan_eps figure: is any candidate transition a real phase transition, and what does eps_s do to it?
  row 1  helix-coil  : helicity theta(J), susceptibility chi_theta = N var(theta), sharpness vs N
  row 2  handedness  : <|m|>(J), Binder cumulant U4(J), location of the crossover vs N
  row 3  chain size  : <Rg^2>(J), its relative fluctuation, and P(m) at J = 3, 4, 5
  row 4  eps_s       : relative change [%] of theta, <|m|>, <Rg^2> against eps_s, one line per J
Blue curves  = EXACT 1D chain without non-bonded interactions (polynomial transfer matrix, weights from
               ../exact_fss.py), N = 100 .. 1600, light -> dark.
Orange marks = MC at N = 200 (sterics + fitted helix-helix potential), ONE MARK PER eps_s, light -> dark with
               eps_s, drawn side by side at each J because they coincide.  Finished runs only.
Errors: max(SEM over seeds, combined block error) for theta, <|m|>, <Rg^2>; SEM over seeds for the fluctuation
quantities.  Style: no grid, ticks inside on all four sides, no title, legend on."""
import contextlib, glob, importlib.util, io, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["x", "../sample_data.dat"]
spec = importlib.util.spec_from_file_location("ex", os.path.abspath("../exact_fss.py")); ex = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(ex)
PAIRS = sorted(ex.idx, key=ex.idx.get); SECOND = np.array([b for (a, b) in PAIRS])

def set_J(J): ex.J0, ex.J1, ex.J2 = float(J), 0.0, float(J)                     # E_HH = -J, E_CH = 0, E_RL = +J
def exact_dist(J, N, helix=False):
    """exact distribution of M = n_R - n_L (support -N..N) or of n_h = n_R + n_L (support 0..N)"""
    set_J(J); T, v0 = ex.build(0.0); f = (lambda s: abs(s)) if helix else (lambda s: s); off = 0 if helix else N
    P = np.zeros((9, (N + 1) if helix else (2 * N + 1)))
    for p, (a, b) in enumerate(PAIRS): P[p, off + f(a) + f(b)] = v0[p]
    for _ in range(N - 2):
        Q = T.T @ P
        for q in range(9): Q[q] = np.roll(Q[q], f(SECOND[q]))
        P = Q / Q.sum()
    p = P.sum(0); return p / p.sum()
def hel_inf(J): set_J(J); return ex.hel_inf(0.0)

# validated ramps (dataviz validate_palette.js --ordinal): one hue each, light -> dark
NS = (100, 200, 400, 800, 1600);  CN = dict(zip(NS, ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]))          # N
CE = dict(zip([0, 1, 3, 4, 7, 8], ["#ef9c77", "#e9713f", "#d1501c", "#a53c13", "#752a0c", "#471806"]))                 # eps_s
CJ = dict(zip([1, 2, 3, 4, 5, 6], ["#52c79e", "#22b07f", "#149166", "#0e7150", "#09523a", "#043325"]))                 # J
CP = {3: "#eb6834", 4: "#1baf7a", 5: "#4a3aa7"}                                                                       # P(m)
INK, MUTED, C1, C2, SURF = "#0b0b0b", "#898781", "#2a78d6", "#eb6834", "#ffffff"

Jg = np.arange(0.5, 8.001, 0.125)
E = {N: {k: np.zeros(len(Jg)) for k in ("th", "chi_th", "am", "U4")} for N in NS}
for N in NS:
    m = np.arange(-N, N + 1) / N; h = np.arange(N + 1) / N
    for i, J in enumerate(Jg):
        pm = exact_dist(J, N); ph = exact_dist(J, N, helix=True)
        m2, m4 = (pm * m**2).sum(), (pm * m**4).sum(); th = (ph * h).sum()
        E[N]["am"][i] = (pm * np.abs(m)).sum(); E[N]["U4"][i] = 1 - m4 / (3 * m2 * m2)
        E[N]["th"][i] = th; E[N]["chi_th"][i] = N * ((ph * h * h).sum() - th * th)
th_inf = np.array([hel_inf(J) for J in Jg])

# ---------------------------------------------------------------- MC, N = 200, per (J, eps_s)
N0 = 200
def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def blk(x, nb=20): b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return b.std(ddof=1) / np.sqrt(nb)
runs = {}
for log in sorted(glob.glob("logs/J*_es*_s*.log")):
    if "summary" not in open(log).read(): continue
    J, es, s = map(int, re.search(r"J(\d+)_es(\d+)_s(\d+)", log).groups())
    o = read_obs("out/" + os.path.basename(log)[:-4] + "_obs.dat"); m = (o["n_R"] - o["n_L"]) / N0
    runs.setdefault((J, es), []).append(dict(m=m, th=o["helicity"].mean(), am=np.abs(m).mean(), rg2=o["Rg2"].mean(),
        e_th=blk(o["helicity"]), e_am=blk(np.abs(m)), e_rg2=blk(o["Rg2"]), chi_th=N0 * o["helicity"].var(),
        U4=1 - (m**4).mean() / (3 * (m**2).mean()**2), rg2_rel=o["Rg2"].std() / o["Rg2"].mean()))
Js = sorted({k[0] for k in runs}); ESs = sorted({k[1] for k in runs})
def cell(J, es, key):
    rs = runs[(J, es)]; v = np.array([r[key] for r in rs]); sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    eb = np.sqrt(sum(r["e_" + key]**2 for r in rs)) / len(rs) if "e_" + key in rs[0] else 0.0
    return v.mean(), max(sem, eb)

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(4, 3, figsize=(12.5, 14.0)); (a, b, c), (d, e, f), (g, h, i_), (j, k, l) = ax
def tag(x, s, loc=(0.04, 0.94)): x.text(*loc, s, transform=x.transAxes, va="top", ha="left", fontsize=9, color=INK)
def mc_by_eps(x, key, legend=False):
    """one mark per eps_s at every J, side by side (they coincide)"""
    for n, es in enumerate(ESs):
        pts = [(J + (n - (len(ESs) - 1) / 2) * 0.09,) + cell(J, es, key) for J in Js if (J, es) in runs]
        if not pts: continue
        X, V, Er = map(np.array, zip(*pts))
        x.errorbar(X, V, yerr=Er, fmt="o", ms=5, mfc=CE[es], mec=SURF, mew=0.6, ecolor=CE[es], elinewidth=0.9, zorder=5,
                   label=rf"MC N = 200, $\epsilon_s$ = {es}" if legend else None)
def jaxis(x, yl): x.set_xlabel(r"coupling $J$ [$k_BT$]"); x.set_ylabel(yl); x.set_xlim(0.5, 8)

for x, key, yl, t in ((a, "th", r"helicity $\theta$", "(a) helix–coil"), (b, "chi_th", r"$\chi_\theta = N\,\mathrm{var}(\theta)$", "(b)"),
                      (d, "am", r"handedness $\langle |m| \rangle$", "(d) handedness"), (e, "U4", r"Binder cumulant $U_4$", "(e)")):
    for N in NS: x.plot(Jg, E[N][key], color=CN[N], lw=1.8, label=f"exact 1D, N = {N}" if key == "th" else None)
    if key == "th": x.plot(Jg, th_inf, color=INK, lw=1.1, ls="--", label=r"exact 1D, N $\to\infty$")
    mc_by_eps(x, key, legend=(key == "am")); jaxis(x, yl); tag(x, t, loc=(0.04, 0.84) if key == "U4" else (0.04, 0.94))
a.legend(loc="lower right"); d.legend(loc="lower right")
e.axhline(2 / 3, color=MUTED, lw=0.8, ls=":"); e.text(0.7, 0.655, "2/3 = fully ordered", color=MUTED, fontsize=8, va="top")
e.text(0.97, 0.05, "curves of different N never cross", transform=e.transAxes, ha="right", fontsize=8, color=INK)

Jf = np.linspace(Jg[0], Jg[-1], 6001)
sl_th = [np.gradient(E[N]["th"], Jg).max() for N in NS]; sl_am = [np.gradient(E[N]["am"], Jg).max() for N in NS]
js_th = [Jf[np.argmin(np.abs(np.interp(Jf, Jg, E[N]["th"]) - 0.5))] for N in NS]; js_am = [Jf[np.argmin(np.abs(np.interp(Jf, Jg, E[N]["am"]) - 0.5))] for N in NS]
for x, y1, y2 in ((c, sl_th, sl_am), (f, js_th, js_am)):
    x.plot(NS, y1, "o-", color=C1, lw=1.8, ms=6, label="helix–coil"); x.plot(NS, y2, "s--", color=C2, lw=1.8, ms=6, label="handedness")
    x.set_xscale("log"); x.set_xlabel("chain length $N$ (exact 1D)"); x.legend(loc="lower right")
c.set_ylim(0, 1.0); c.set_ylabel(r"steepest slope, max $d\,\cdot/dJ$  [1/$k_BT$]"); tag(c, "(c) sharpness does not grow with N")
kf = np.polyfit(np.log(NS), js_am, 1); f.set_ylim(0, 7.5); f.set_ylabel(r"crossover location $J^*$ (value = 0.5)  [$k_BT$]")
tag(f, "(f) handedness crossover drifts with N"); tag(f, rf"handedness: $J^* = {kf[1]:.2f} + {kf[0]:.2f}\,\ln N$", loc=(0.04, 0.86))

mc_by_eps(g, "rg2"); g.set_ylim(0, 10500); jaxis(g, r"$\langle R_g^2 \rangle$ [code units$^2$]"); tag(g, "(g) chain size (MC only)")
mc_by_eps(h, "rg2_rel"); h.set_ylim(0, 0.7); jaxis(h, r"std$(R_g^2)\,/\,\langle R_g^2 \rangle$"); tag(h, "(h) size fluctuation (MC only)")

bins = np.linspace(-1, 1, 42); ctr = 0.5 * (bins[1:] + bins[:-1]); mgrid = np.arange(-N0, N0 + 1) / N0
nint = np.histogram(mgrid, bins=bins)[0]                       # integer M values per bin (9 or 10): density per value, no zig-zag
for J in (3, 4, 5):
    if J not in Js: continue
    pe = np.histogram(mgrid, bins=bins, weights=exact_dist(J, N0))[0] / nint * N0
    mm = np.concatenate([r["m"] for (jj, ee), rs in runs.items() if jj == J for r in rs]); pmc = np.histogram(mm, bins=bins)[0] / len(mm) / nint * N0
    pmc[pmc == 0] = np.nan; i_.plot(ctr, pe, color=CP[J], lw=1.8); i_.plot(ctr, pmc, "o", ms=3.5, mfc=SURF, mec=CP[J], mew=1.0)
    i_.text(0.0, pe[len(ctr) // 2] * 1.3, f"J = {J}", color=INK, fontsize=8, ha="center", va="bottom")
i_.plot([], [], color=INK, lw=1.8, label="exact 1D, N = 200"); i_.plot([], [], "o", ms=3.5, mfc=SURF, mec=INK, label=r"MC, all $\epsilon_s$ pooled")
i_.set_yscale("log"); i_.set_ylim(3e-3, 80); i_.legend(loc="upper center", ncol=2, columnspacing=1.0, handlelength=1.5)
i_.set_xlabel(r"$m = (n_R - n_L)/N$"); i_.set_ylabel(r"$P(m)$"); i_.set_xlim(-1, 1); tag(i_, "(i)", loc=(0.04, 0.82))

REL = {}
for x, key, yl, t in ((j, "th", r"$\Delta\theta/\theta$  [%]", r"(j) helicity vs $\epsilon_s$"), (k, "am", r"$\Delta\langle|m|\rangle/\langle|m|\rangle$  [%]", r"(k) handedness vs $\epsilon_s$"),
                      (l, "rg2", r"$\Delta\langle R_g^2\rangle/\langle R_g^2\rangle$  [%]", r"(l) chain size vs $\epsilon_s$")):
    for n, J in enumerate(Js):
        if (J, 0) not in runs: continue
        v0, e0 = cell(J, 0, key); pts = [(es + (n - 2.5) * 0.05,) + cell(J, es, key) for es in ESs if es != 0 and (J, es) in runs]
        if not pts: continue
        X, V, Er = map(np.array, zip(*pts)); rel, rerr = 100 * (V - v0) / v0, 100 * np.hypot(Er, e0) / v0
        REL[(key, J)] = (np.array([es for es in ESs if es != 0 and (J, es) in runs]), rel, rerr)
        x.errorbar(X, rel, yerr=rerr, fmt="o-", ms=4.5, lw=1.2, color=CJ[J], mfc=CJ[J], mec=SURF, mew=0.5, elinewidth=0.9, label=f"J = {J}")
    x.axhline(0, color=MUTED, lw=0.8); x.set_xlim(0, max(ESs) + 1); x.set_xlabel(r"chiral amplitude $\epsilon_s$"); x.set_ylabel(yl + r"   relative to $\epsilon_s = 0$"); tag(x, t)
    yl_ = 1.45 * max(abs(np.array(x.get_ylim()))); x.set_ylim(-yl_, yl_)      # headroom for the tag and the legend
j.legend(loc="lower left", ncol=3, columnspacing=1.0, handlelength=1.5)
fig.tight_layout(); fig.savefig("transitions.png", dpi=170); fig.savefig("transitions.pdf")

# ---------------------------------------------------------------- numbers behind the figure
print(f"MC runs used: {sum(len(v) for v in runs.values())};  eps_s values: {ESs}")
print("\n     N   max dtheta/dJ   J*(theta=0.5)   max chi_theta   |  max d<|m|>/dJ   J*(<|m|>=0.5)")
for n, N in enumerate(NS): print(f"  {N:5d}   {sl_th[n]:13.3f}   {js_th[n]:13.3f}   {E[N]['chi_th'].max():13.3f}   |  {sl_am[n]:14.3f}   {js_am[n]:13.3f}")
print(f"  N->inf: max dtheta/dJ = {np.gradient(th_inf, Jg).max():.3f}  (finite -> crossover)")
for key, name in (("th", "helicity"), ("am", "<|m|>"), ("rg2", "<Rg2>")):
    print(f"\n  relative change of {name} vs eps_s = 0, in %  (value +- error)")
    print("     J  " + "".join(f"{'eps_s=' + str(es):>18s}" for es in ESs if es != 0))
    comb = {es: [] for es in ESs if es != 0}
    for J in Js:
        if (key, J) not in REL: continue
        es_, rel, rerr = REL[(key, J)]; dd = dict(zip(es_, zip(rel, rerr)))
        for es in dd: comb[es].append(dd[es][0] / dd[es][1])
        print(f"  {J:4d}  " + "".join(f"{dd[es][0]:+9.3f} +- {dd[es][1]:5.3f}" if es in dd else f"{'--':>18s}" for es in ESs if es != 0))
    print("   all J combined, sum(z)/sqrt(n)  (|.| > 3 = a significant common shift):  "
          + "   ".join(f"eps_s={es}: {np.sum(z) / np.sqrt(len(z)):+.2f}" for es, z in comb.items() if z))
