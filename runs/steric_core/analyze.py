#!/usr/bin/env python3
"""Gyration-tensor eigenvalues of the soft-core coil (x beads) vs the WCA chain (7x beads).
lambda_1 >= lambda_2 >= lambda_3 per frame from the dumped configurations, averaged; errors from 10 time blocks
combined over seeds. If one soft bead ~ 7 WCA monomers, lambda_i^WCA(7x) / lambda_i^gauss(x) is one common
constant (length-scale factor^2) for all i and x, and the shape ratios lambda_2/lambda_1, lambda_3/lambda_1 agree."""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def read_frames(f):
    frames = []
    with open(f) as fh:
        while True:
            head = fh.readline()
            if not head: break
            n = int(head); fh.readline()
            frames.append([[float(v) for v in fh.readline().split()[1:4]] for _ in range(n)])
    return np.array(frames)

def read_obs(f):
    lines = open(f).read().splitlines()
    hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}

def tau_int(x):
    x = np.asarray(x, float); x = x - x.mean(); n = len(x)
    if n < 16 or x.std() == 0: return np.nan
    f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); ac /= ac[0]
    t = 0.5
    for M in range(1, n // 2):
        t += ac[M]
        if M > 6 * t: break
    return t

runs = {}
for f in sorted(glob.glob("out/*_conf.xyz")):
    m = re.match(r"out/(gauss|wca)_x(\d+)_s(\d+)_conf", f); model, x = m.group(1), int(m.group(2))
    P = read_frames(f)
    if len(P) < 20: continue
    P = P[len(P) // 5:]                                           # drop the first 20 % as extra burn-in
    lam = []
    for cfg in P:
        c = cfg - cfg.mean(0); T = c.T @ c / len(c)
        lam.append(np.sort(np.linalg.eigvalsh(T))[::-1])
    lam = np.array(lam); nb = 10; blk = lam[: len(lam) // nb * nb].reshape(nb, -1, 3).mean(1)
    o = read_obs(f.replace("_conf.xyz", "_obs.dat")); dt = o["sweep"][1] - o["sweep"][0]
    runs.setdefault((model, x), []).append(dict(lam=lam.mean(0), err=blk.std(0, ddof=1) / np.sqrt(nb), nfr=len(lam),
                                              tau=tau_int(o["Rg2"]) * dt, nsw=o["sweep"][-1], acc_piv=o["acc_pivot"][-1]))

res = {}
for key, rs in sorted(runs.items()):
    n = len(rs)
    res[key] = (np.mean([r["lam"] for r in rs], 0), np.sqrt(np.sum([r["err"] ** 2 for r in rs], 0)) / n, rs)
xs = sorted({x for (_, x) in res})
print(f"{'model':>6} {'x':>4} {'N':>4} {'seeds':>5} {'lam1':>9} {'+-':>6} {'lam2':>8} {'+-':>6} {'lam3':>8} {'+-':>6} {'Rg2':>8} {'l2/l1':>6} {'l3/l1':>6} {'run/tau':>8} {'acc_piv':>7}")
for model in ("gauss", "wca"):
    for x in xs:
        if (model, x) not in res: continue
        lam, err, rs = res[(model, x)]; N = x if model == "gauss" else 7 * x
        tau = np.nanmean([r["tau"] for r in rs])
        print(f"{model:>6} {x:>4} {N:>4} {len(rs):>5} {lam[0]:>9.2f} {err[0]:>6.2f} {lam[1]:>8.2f} {err[1]:>6.2f} {lam[2]:>8.2f} {err[2]:>6.2f} {lam.sum():>8.2f} {lam[1]/lam[0]:>6.3f} {lam[2]/lam[0]:>6.3f} {rs[0]['nsw']/tau:>8.1f} {rs[0]['acc_piv']:>7.2f}")
print("\nratio lambda_i(WCA, 7x) / lambda_i(gauss, x)   [should be one x-independent constant for all i]")
print(f"{'x':>4} {'lam1':>8} {'+-':>6} {'lam2':>8} {'+-':>6} {'lam3':>8} {'+-':>6}")
ratios = {}
for x in xs:
    if ("gauss", x) in res and ("wca", x) in res:
        lg, eg, _ = res[("gauss", x)]; lw, ew, _ = res[("wca", x)]
        r = lw / lg; e = r * np.sqrt((ew / lw) ** 2 + (eg / lg) ** 2); ratios[x] = (r, e)
        print(f"{x:>4} {r[0]:>8.3f} {e[0]:>6.3f} {r[1]:>8.3f} {e[1]:>6.3f} {r[2]:>8.3f} {e[2]:>6.3f}")

# power-law exponents 2nu from lambda_1
for model in ("gauss", "wca"):
    pts = [(x, res[(model, x)][0][0]) for x in xs if (model, x) in res]
    if len(pts) >= 3:
        p = np.polyfit(np.log([q[0] for q in pts]), np.log([q[1] for q in pts]), 1)
        print(f"{model}: lambda_1 ~ x^{p[0]:.3f}  (self-avoiding: 2nu = 1.176; ideal: 1)")

# ---------------- figure ----------------
BLUE, ORANGE, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, which="major", color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.set_xscale("log", base=2); ax.set_xticks(xs); ax.set_xticklabels([str(x) for x in xs])
    ax.set_xlabel("x   (soft beads; WCA chain has 7x beads)", color=INK)
ax = axs[0]; ax.set_yscale("log")
mk = ["o", "s", "^"]
for model, col, lab in (("gauss", BLUE, "soft core, x beads"), ("wca", ORANGE, "WCA, 7x beads")):
    for i in range(3):
        pts = [(x, res[(model, x)][0][i], res[(model, x)][1][i]) for x in xs if (model, x) in res]
        if pts:
            ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[p[2] for p in pts], fmt=mk[i] + "-", ms=5, lw=1,
                        color=col, mfc=(col if model == "wca" else BG), capsize=2, label=f"{lab}: λ{i+1}" if True else None)
ax.set_ylabel("gyration-tensor eigenvalues", color=INK); ax.set_title("λ₁ ≥ λ₂ ≥ λ₃ vs chain length", loc="left", fontsize=10, color=INK)
ax.legend(frameon=False, fontsize=7.5, ncol=2)
ax = axs[1]
for i in range(3):
    pts = [(x, ratios[x][0][i], ratios[x][1][i]) for x in xs if x in ratios]
    if pts: ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[p[2] for p in pts], fmt=mk[i] + "-", ms=5, lw=1, color=INK2, capsize=2, label=f"λ{i+1}")
ax.set_ylabel("λᵢ(WCA, 7x) / λᵢ(soft, x)", color=INK); ax.set_title("Length-scale factor between the models", loc="left", fontsize=10, color=INK)
ax.legend(frameon=False, fontsize=8.5)
ax = axs[2]
for model, col, lab in (("gauss", BLUE, "soft core"), ("wca", ORANGE, "WCA")):
    for i, ls in ((1, "-"), (2, "--")):
        pts = [(x, res[(model, x)][0][i] / res[(model, x)][0][0]) for x in xs if (model, x) in res]
        if pts: ax.plot([p[0] for p in pts], [p[1] for p in pts], mk[i] + ls, ms=5, lw=1, color=col, mfc=(col if model == "wca" else BG), label=f"{lab}: λ{i+1}/λ₁")
ax.set_ylabel("shape ratios", color=INK); ax.set_title("Shape anisotropy (dimensionless)", loc="left", fontsize=10, color=INK)
ax.set_ylim(0, 0.35); ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("steric_core.png")
