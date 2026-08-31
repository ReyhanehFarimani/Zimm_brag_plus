#!/usr/bin/env python3
"""Average Rg^2 per run, then mean +/- SEM over seeds per N; write rg_scaling.csv and rg_scaling.png."""
import argparse, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
ap = argparse.ArgumentParser()
ap.add_argument("--data", default="out", help="directory with *_obs.dat files")
ap.add_argument("--tag", default="", help="suffix for output files")
args = ap.parse_args()
DATA, TAG = args.data, args.tag

# columns: sweep n_R n_L helicity E_bond Ree Rg2 acc_state acc_pos
def tau_int(x):
    """Integrated autocorrelation time in units of samples (Sokal automatic window)."""
    x = np.asarray(x, float); x = x - x.mean(); n = len(x)
    if n < 16 or x.std() == 0: return float("nan")
    f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); ac /= ac[0]
    t = 0.5
    for M in range(1, n // 2):
        t += ac[M]
        if M > 6 * t: break
    return t

runs, series, taus = {}, {}, {}
for f in sorted(glob.glob(f"{DATA}/N*_s*_obs.dat")):
    m = re.search(r"N(\d+)_s(\d+)_obs\.dat$", f)
    N, s = int(m.group(1)), int(m.group(2))
    # parse by hand so a file still being written (torn last line) doesn't break the analysis;
    # column positions come from the header line so old and new file layouts both work
    lines = open(f).read().splitlines()
    hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    ncol, i_rg = len(hdr), hdr.index("Rg2")
    rows = [(float(t[0]), float(t[i_rg])) for t in (l.split() for l in lines) if len(t) == ncol and t[0] != "#"]
    if not rows:
        continue
    t, rg2 = np.array(rows).T
    runs.setdefault(N, []).append(float(rg2.mean()))
    series.setdefault(N, []).append((s, t, rg2))
    dt = t[1] - t[0] if len(t) > 1 else 1.0
    taus.setdefault(N, []).append((tau_int(rg2) * dt, len(rg2) * dt))

Ns = np.array(sorted(runs))
mean = np.array([np.mean(runs[N]) for N in Ns])
sem  = np.array([np.std(runs[N], ddof=1) / np.sqrt(len(runs[N])) if len(runs[N]) > 1 else 0.0 for N in Ns])
nrep = np.array([len(runs[N]) for N in Ns])

# Freely-jointed reference (bonds only, no angular terms => independent bond vectors):
#   <Rg2> = b^2 (N^2 - 1) / (6 N),  b^2 = <r^2> of one harmonic bond in 3D (with the r^2 Jacobian)
R0, KB, KT = 1.0, 100.0, 1.0
r = np.linspace(max(0.0, R0 - 12 / np.sqrt(KB / KT)), R0 + 12 / np.sqrt(KB / KT), 20001)
w = r**2 * np.exp(-0.5 * KB / KT * (r - R0) ** 2)
b2 = np.trapezoid(r**2 * w, r) / np.trapezoid(w, r)
theory = lambda N: b2 * (N**2 - 1) / (6 * N)
Nline = np.geomspace(Ns[0] / 1.3, Ns[-1] * 1.3, 50)
ref = theory(Nline)

with open(f"rg_scaling{TAG}.csv", "w") as fh:
    fh.write("N,n_seeds,Rg2_mean,Rg2_sem,Rg2_theory,Rg2_per_seed\n")
    for N, m_, e_, n_ in zip(Ns, mean, sem, nrep):
        fh.write(f"{N},{n_},{m_:.4f},{e_:.4f},{theory(N):.4f},\"{' '.join(f'{x:.3f}' for x in runs[N])}\"\n")

# power-law fit  Rg2 = A N^(2nu)
p = np.polyfit(np.log(Ns), np.log(mean), 1)
two_nu, logA = p


# ---- figure ----
BLUE, INK, INK2, GRID = "#2a78d6", "#0b0b0b", "#52514e", "#e6e5e1"
fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=160, facecolor="#fcfcfb")
ax.set_facecolor("#fcfcfb")
ax.plot(Nline, ref, ls="--", lw=1.4, color=INK2, zorder=1)
ax.errorbar(Ns, mean, yerr=sem, fmt="o", ms=6, color=BLUE, ecolor=BLUE,
            elinewidth=1.4, capsize=3, zorder=3, label=None)
for N, v in runs.items():                                   # individual seeds, faint
    ax.plot([N] * len(v), v, "o", ms=3, mfc="none", mec=BLUE, alpha=0.45, zorder=2)
ax.set_xscale("log", base=2); ax.set_yscale("log")
ax.set_xticks(Ns); ax.set_xticklabels([str(N) for N in Ns])
ax.set_xlabel("N (residues)", color=INK)
ax.set_ylabel(r"$\langle R_g^2\rangle$  [bond length$^2$]", color=INK)
ax.set_title("Radius of gyration vs chain length", loc="left", color=INK, fontsize=11)
ax.text(Nline[-1], ref[-1] * 1.12, r"freely-jointed  $b^2(N^2-1)/6N$", ha="right", va="bottom",
        fontsize=8.5, color=INK2)
seeds_lbl = str(nrep.min()) if nrep.min() == nrep.max() else f"{nrep.min()}–{nrep.max()}"
ax.text(0.02, 0.97, f"fit: $R_g^2 \\propto N^{{{two_nu:.2f}}}$   ({seeds_lbl} seeds per N, error bar = SEM over seeds)",
        transform=ax.transAxes, va="top", fontsize=8.5, color=INK2)
ax.grid(True, which="major", color=GRID, lw=0.8); ax.grid(False, which="minor")
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
ax.tick_params(colors=INK2)
fig.tight_layout(); fig.savefig(f"rg_scaling{TAG}.png")

# ---- time series: Rg2(t) / theory for every seed, one panel per N ----
fig, axes = plt.subplots(len(Ns), 1, figsize=(6.4, 1.6 * len(Ns) + 0.6), dpi=160,
                         sharex=True, facecolor="#fcfcfb")
for ax, N in zip(np.atleast_1d(axes), Ns):
    ax.set_facecolor("#fcfcfb")
    for s_, t, rg2 in sorted(series[N]):
        ax.plot(t / 1e6, rg2 / theory(N), lw=0.9, color=BLUE, alpha=0.75)
    ax.axhline(1.0, ls="--", lw=1.0, color=INK2)
    ax.text(0.01, 0.92, f"N = {N}", transform=ax.transAxes, va="top", fontsize=9, color=INK)
    ax.set_ylim(0, 3); ax.set_yticks([0, 1, 2, 3])
    ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
np.atleast_1d(axes)[0].set_title("Rg² over the run, per seed  (divided by freely-jointed value)",
                                 loc="left", fontsize=10, color=INK)
np.atleast_1d(axes)[-1].set_xlabel("sweeps  [×10⁶]", color=INK)
fig.tight_layout(); fig.savefig(f"rg_timeseries{TAG}.png")

print(f"b^2 = <r^2> = {b2:.4f}")
print(f"{'N':>6} {'seeds':>5} {'<Rg2>':>10} {'SEM':>8} {'theory':>8} {'ratio':>7} {'tau_int':>10} {'run/tau':>8} {'n_indep':>8}")
for N, m_, e_, n_ in zip(Ns, mean, sem, nrep):
    tau = np.nanmean([a for a, _ in taus[N]]); L = np.mean([b for _, b in taus[N]])
    print(f"{N:>6} {n_:>5} {m_:>10.2f} {e_:>8.2f} {theory(N):>8.2f} {m_/theory(N):>7.3f} {tau:>10.3g} {L/tau:>8.1f} {n_*L/(2*tau):>8.0f}")
print(f"power-law exponent 2nu = {two_nu:.3f}")
