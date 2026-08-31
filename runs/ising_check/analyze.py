#!/usr/bin/env python3
"""Compare the MC state statistics with the exact transfer-matrix solution of the 3-state chain."""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))

ST = np.array([0, 1, -1])            # C, R, L  (spin values)

def pair_E(a, b, J0, J1, J2):
    if a == 0 and b == 0: return 0.0
    if a == 0 or b == 0:  return J1
    return -J0 if a == b else J2

def tmatrix(beta, J0, J1, J2):
    return np.array([[np.exp(-beta * pair_E(a, b, J0, J1, J2)) for b in ST] for a in ST])

def lnZ(N, beta, J0, J1, J2):
    """log partition function of the open chain, with per-step renormalisation (no overflow at large N)."""
    T = tmatrix(beta, J0, J1, J2)
    v, lg = np.ones(3), 0.0
    for _ in range(1, N):
        v = v @ T; m = v.max(); v /= m; lg += np.log(m)
    return lg + np.log(v.sum())

def exact(N, kT, J0, J1, J2, rmax=None):
    """Return helicity, <E>/(N-1), and <s_i s_{i+r}> (averaged over i) for the open chain."""
    beta = 1.0 / kT
    T = tmatrix(beta, J0, J1, J2)
    T = T / np.abs(np.linalg.eigvals(T)).max()       # rescale: ratios of weights are unchanged, no overflow
    L = [np.ones(3)]
    for _ in range(1, N): L.append(L[-1] @ T)
    R = [np.ones(3)]
    for _ in range(1, N): R.append(T @ R[-1])
    Z = L[N - 1] @ np.ones(3)
    hel = sum((L[i] * R[N - 1 - i])[1:].sum() for i in range(N)) / Z / N
    h = 1e-5
    E = -(lnZ(N, beta + h, J0, J1, J2) - lnZ(N, beta - h, J0, J1, J2)) / (2 * h)
    S = np.diag(ST)
    corr, Tr = [], np.eye(3)
    for r in range(0, (N if rmax is None else rmax + 1)):
        c = np.mean([L[i] @ S @ Tr @ S @ R[N - 1 - i - r] for i in range(N - r)]) / Z
        corr.append(c); Tr = Tr @ T
    return hel, E / (N - 1), np.array(corr)

def read_obs(f):
    lines = open(f).read().splitlines()
    hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}

def read_input(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]
            d[k] = v
    return d

# ---------- scan A: kT ----------
A = []
for f in sorted(glob.glob("out/A_kT*_obs.dat")):
    tag = re.search(r"(A_kT[\d.]+)_obs", f).group(1)
    inp = read_input(f"inputs/{tag}.dat"); N = int(inp["N"]); kT = float(inp["kT"])
    J0, J1, J2 = (float(inp[k]) for k in ("J0", "J1", "J2"))
    o = read_obs(f)
    n = len(o["helicity"]); nb = 20
    blk = lambda x: (x[: n // nb * nb].reshape(nb, -1).mean(1))
    hel_b, e_b = blk(o["helicity"]), blk(o["E_state"] / (N - 1))
    he, Ee, _ = exact(N, kT, J0, J1, J2)
    A.append((kT, hel_b.mean(), hel_b.std(ddof=1) / np.sqrt(nb), e_b.mean(), e_b.std(ddof=1) / np.sqrt(nb), he, Ee))
A = np.array(sorted(A))
kTs = np.linspace(A[0, 0] * 0.9, A[-1, 0] * 1.05, 120)
ex = np.array([exact(N, k, J0, J1, J2, rmax=0)[:2] for k in kTs])

# ---------- scan B: correlation function ----------
Bfiles = sorted(glob.glob("out/B_*_conf.xyz"))
inpB = read_input("inputs/" + re.search(r"(B_\w+)_conf", Bfiles[0]).group(1) + ".dat")
NB = int(inpB["N"]); kTB = float(inpB["kT"]); J0B, J1B, J2B = (float(inpB[k]) for k in ("J0", "J1", "J2"))
rmax = 20
def read_frames(f):
    frames = []
    with open(f) as fh:
        while True:
            head = fh.readline()
            if not head: break
            n = int(head); fh.readline()
            frames.append([{"R": 1, "L": -1}.get(fh.readline()[0], 0) for _ in range(n)])
    return np.array(frames, float)
per_seed = []                                     # correlation function of each independent run
for f in Bfiles:
    S = read_frames(f)
    per_seed.append([np.mean(S[:, :NB - r] * S[:, r:]) for r in range(rmax + 1)])
per_seed = np.array(per_seed); nseed, nf = len(per_seed), len(S)
corr_mc = per_seed.mean(0)
corr_err = per_seed.std(0, ddof=1) / np.sqrt(nseed)   # error over independent seeds
_, _, corr_ex = exact(NB, kTB, J0B, J1B, J2B, rmax)
Jeff = 0.5 * (J0B + J2B)
ising2 = np.tanh(Jeff / kTB) ** np.arange(rmax + 1)
coil_frac = 1 - np.mean([read_obs(f.replace("_conf.xyz", "_obs.dat"))["helicity"].mean() for f in Bfiles])

# ---------- report ----------
print(f"scan A: N={N}  J0={J0} J1={J1} J2={J2}")
print(f"{'kT':>5} {'hel_MC':>9} {'+-':>7} {'hel_exact':>10} {'dev/sig':>8} {'E_MC':>9} {'+-':>7} {'E_exact':>9} {'dev/sig':>8}")
for kT, hm, he_, em, ee_, hx, Ex in A:
    print(f"{kT:>5.2f} {hm:>9.5f} {he_:>7.5f} {hx:>10.5f} {(hm-hx)/he_:>8.2f} {em:>9.5f} {ee_:>7.5f} {Ex:>9.5f} {(em-Ex)/ee_:>8.2f}")
print(f"\nscan B: N={NB} kT={kTB} J0={J0B} J1={J1B} J2={J2B}  seeds={nseed} x {nf} frames  coil fraction={coil_frac:.4f}")
print(f"{'r':>3} {'<s s>_MC':>10} {'+-':>7} {'exact3':>8} {'ising2 tanh^r':>14}")
for r in range(0, rmax + 1, 2):
    print(f"{r:>3} {corr_mc[r]:>10.4f} {corr_err[r]:>7.4f} {corr_ex[r]:>8.4f} {ising2[r]:>14.4f}")

# ---------- figure ----------
BLUE, ORANGE, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.0), dpi=160, facecolor=BG)
def style(ax):
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5)
ax = axs[0]; style(ax)
ax.plot(kTs, ex[:, 0], "-", lw=1.6, color=INK2, label="exact transfer matrix")
ax.errorbar(A[:, 0], A[:, 1], yerr=A[:, 2], fmt="o", ms=5, color=BLUE, capsize=2.5, label="MC")
ax.set_xlabel("kT", color=INK); ax.set_ylabel("helicity  (fraction R or L)", color=INK)
ax.set_title(f"Helicity vs kT   (J0={J0}, J1={J1}, J2={J2}, N={N})", loc="left", fontsize=9.5, color=INK)
ax.legend(frameon=False, fontsize=8.5)
ax = axs[1]; style(ax)
ax.plot(kTs, ex[:, 1], "-", lw=1.6, color=INK2, label="exact transfer matrix")
ax.errorbar(A[:, 0], A[:, 3], yerr=A[:, 4], fmt="o", ms=5, color=BLUE, capsize=2.5, label="MC")
ax.set_xlabel("kT", color=INK); ax.set_ylabel(r"$\langle E_{state}\rangle$ per bond", color=INK)
ax.set_title("State energy vs kT", loc="left", fontsize=9.5, color=INK)
ax.legend(frameon=False, fontsize=8.5)
ax = axs[2]; style(ax)
r = np.arange(rmax + 1)
ax.plot(r, corr_ex, "-", lw=1.6, color=INK2, label="exact 3-state")
ax.plot(r, ising2, "--", lw=1.4, color=ORANGE, label=rf"2-state Ising  tanh$^r$($J_{{eff}}/kT$), $J_{{eff}}$=(J0+J2)/2")
ax.errorbar(r, corr_mc, yerr=corr_err, fmt="o", ms=4.5, color=BLUE, capsize=2, label=f"MC ({nseed} seeds)")
ax.set_yscale("log"); ax.set_ylim(min(corr_ex[rmax], ising2[rmax]) * 0.5, 1.3)
ax.set_xlabel("r  (residue separation)", color=INK); ax.set_ylabel(r"$\langle s_i\, s_{i+r}\rangle$", color=INK)
ax.set_title(f"Spin correlation  (N={NB}, kT={kTB})", loc="left", fontsize=9.5, color=INK)
ax.text(0.98, 0.95, f"J0={J0B}, J1={J1B}, J2={J2B}\ncoil fraction {100*coil_frac:.2f}%", transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("ising_check.png")
