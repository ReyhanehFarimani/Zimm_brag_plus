#!/usr/bin/env python3
"""Fig. 2 of the paper: <Re^2>_id / <Re^2>_0 vs M, MC vs eq 7 (exact rigid-rod theory) for J = 0.5 and 2.
Also <N_d> vs eq 9. Our exact finite-kappa reference is evaluated for M <= 300 (O(M^2) transfer-matrix sums)."""
import glob, os, re, sys
import numpy as np
from math import comb, lgamma, log, exp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv = [sys.argv[0]]
exec(open("../paper_ideal/analyze.py").read().split("# ---------------- collect runs")[0])   # reuse readers, exact_model, tau_int

def eq7(M, J):
    """<Re^2>_id / l0^2 (eq 7), evaluated in log space to avoid overflow at large M."""
    c = np.exp(2 * J); lc = 2 * J
    terms = []
    for n in range(1, M + 1):
        lt = (M - n) * lc + lgamma(M) - lgamma(n) - lgamma(M - n + 1) + log(2 * M - (n - 1)) - log(n + 1)
        terms.append(lt)
    m = max(terms)
    return M * exp(m - (M - 1) * np.log1p(c)) * sum(exp(t - m) for t in terms)

def eq8(M, J): return M * (1 + 2 * np.exp(2 * J))
def eq9(M, J): return (M - 1) * np.exp(-J) / (np.exp(J) + np.exp(-J))

runs = {}
for f in sorted(glob.glob("out/M*_obs.dat")):
    m = re.search(r"M(\d+)_J([\d.]+)_s(\d+)_obs", f); M, J = int(m.group(1)), float(m.group(2))
    o = read_obs(f); n = len(o["Ree"])
    if n < 100: continue
    Re2 = o["Ree"] ** 2
    dt = o["sweep"][1] - o["sweep"][0]; tau = tau_int(Re2) * dt
    nb = 20; blk = Re2[: n // nb * nb].reshape(nb, -1).mean(1); Re2_err = blk.std(ddof=1) / np.sqrt(nb)   # blocking error
    frames = f.replace("_obs.dat", "_conf.xyz")
    Nd = RB = RB_err = np.nan
    if os.path.exists(frames):
        S = read_frames(frames)[:, :M]
        if len(S) >= 20:
            Nd = np.mean(np.sum(S[:, :-1] * S[:, 1:] < 0, axis=1))
            # Rao-Blackwell estimator: given the spins, rigid rods of m_k bonds are freely jointed, <Re^2|spins> = sum_k m_k^2
            rb = []
            for row in S:
                walls = np.flatnonzero(row[:-1] * row[1:] < 0)
                lens = np.diff(np.concatenate(([0], walls + 1, [M])))
                rb.append(np.sum(lens ** 2))
            rb = np.array(rb, float); b2 = rb[: len(rb) // 10 * 10].reshape(10, -1).mean(1)
            RB, RB_err = rb.mean(), b2.std(ddof=1) / np.sqrt(10)
    runs.setdefault((J, M), []).append(dict(Re2=Re2.mean(), Re2_err=Re2_err, tau=tau, nsw=o["sweep"][-1], Nd=Nd,
                                            RB=RB, RB_err=RB_err, acc_h=o["acc_hinge"][-1]))

print("ratio = <Re^2>/<Re^2>_0 ; 'direct' = measured end-to-end distance, 'RB' = <Re^2 | spins> = sum of squared rod lengths;")
print("errors: blocking within each run, combined over seeds (seed-to-seed scatter shown as 'scat')")
print(f"{'J':>4} {'M':>5} {'sd':>3} {'direct':>8} {'+-':>7} {'scat':>7} {'RB':>8} {'+-':>7} {'eq7':>8} {'dev_dir':>7} {'dev_RB':>7} {'<Nd>':>8} {'eq9':>8} {'acc_h':>6}")
res = {}
for (J, M), rs in sorted(runs.items()):
    n = len(rs); norm = eq8(M, J)
    y = np.array([r["Re2"] for r in rs]) / norm; ye = np.sqrt(np.sum([(r["Re2_err"] / norm) ** 2 for r in rs])) / n
    scat = y.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    rb = np.array([r["RB"] for r in rs]) / norm; rbe = np.sqrt(np.nansum([(r["RB_err"] / norm) ** 2 for r in rs])) / n
    th = eq7(M, J) / norm; nd = np.nanmean([r["Nd"] for r in rs])
    res[(J, M)] = (y.mean(), ye, th, np.nanmean(rb), rbe)
    print(f"{J:>4} {M:>5} {n:>3} {y.mean():>8.4f} {ye:>7.4f} {scat:>7.4f} {np.nanmean(rb):>8.4f} {rbe:>7.4f} {th:>8.4f} {(y.mean()-th)/ye:>7.2f} {(np.nanmean(rb)-th)/rbe:>7.2f} {nd:>8.2f} {eq9(M,J):>8.2f} {rs[0]['acc_h']:>6.2f}")

# ---------- figure (paper Fig. 2 layout) ----------
BLUE, ORANGE, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=160, facecolor=BG); ax.set_facecolor(BG)
Ms = np.arange(50, 1600, 10)
for J, col, lab in ((0.5, BLUE, "J = 0.5"), (2.0, ORANGE, "J = 2")):
    ax.plot(Ms, [eq7(int(M), J) / eq8(int(M), J) for M in Ms], "-", lw=1.5, color=col, alpha=0.9, label=f"theory (eq 7), {lab}")
    pts = sorted((M, v) for (Jj, M), v in res.items() if Jj == J)
    if pts:
        ax.errorbar([p[0] - 8 for p in pts], [p[1][0] for p in pts], yerr=[p[1][1] for p in pts], fmt="o", ms=5, mfc=BG,
                    color=col, capsize=2.5, label=f"MC direct, {lab}")
        ax.errorbar([p[0] + 8 for p in pts], [p[1][3] for p in pts], yerr=[p[1][4] for p in pts], fmt="s", ms=4.5,
                    color=col, capsize=2.5, label=f"MC ⟨Rₑ²|spins⟩, {lab}")
ax.set_xlabel("M  (number of bonds)", color=INK); ax.set_ylabel(r"$\langle R_e^2\rangle_{id}\,/\,\langle R_e^2\rangle_0$", color=INK)
ax.set_title("Ideal chain with reversible hinges — paper Fig. 2", loc="left", fontsize=10.5, color=INK)
ax.set_xlim(0, 1600); ax.set_ylim(0.7, 1.01); ax.grid(True, color=GRID, lw=0.8)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
ax.tick_params(colors=INK2); ax.legend(frameon=False, fontsize=8.5, loc="lower right")
fig.tight_layout(); fig.savefig("paper_fig2.png")
