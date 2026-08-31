#!/usr/bin/env python3
"""MC vs exact results for the freely-jointed chain with reversible hinges (ideal chain, no excluded volume).

Reference 1 (exact for the simulated model): 3-state transfer matrix with the hinge partition function folded
into the pair weights, finite kappa, coil included, dangling last bead included.
Reference 2 (paper, Chauhan/Mueller/Daoulas 2025): 2-state NN Ising with J = J0 + 1/2 ln(w_r/w_f):
  eq 7  <Re^2>_id,  eq 8 <Re^2>_0 = M(1+2c),  eq 9 <N_d>,  eq 13 C(s) = exp(-s/xi),  eq 16 <S^2>/M.
"""
import argparse, glob, os, re, sys
import numpy as np
from math import comb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--data", default="out"); ap.add_argument("--inputs", default="inputs")
ap.add_argument("--burn", type=float, default=0.0, help="fraction of each run to discard")
ap.add_argument("--tag", default="")
args = ap.parse_args()

ST = np.array([0, 1, -1])   # C, R, L

def read_input(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v
    return d

def read_obs(f):
    lines = open(f).read().splitlines()
    hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}

def read_frames(f):
    frames = []
    with open(f) as fh:
        while True:
            head = fh.readline()
            if not head: break
            n = int(head); fh.readline()
            frames.append([{"R": 1, "L": -1}.get(fh.readline()[0], 0) for _ in range(n)])
    return np.array(frames, float)

def hinge_integrals(kappa, theta0=np.pi):
    """w_r/w_f (partition function of a stiff hinge relative to a free one) and <cos(bond-bond angle)> of the
    stiff hinge. theta is the valence angle at the bead (straight = pi), so cos(bond-bond angle) = -cos(theta)."""
    t = np.linspace(0, np.pi, 400001); w = np.sin(t) * np.exp(-0.5 * kappa * (t - theta0)**2)
    Z = np.trapezoid(w, t)
    return 0.5 * Z, -np.trapezoid(np.cos(t) * w, t) / Z

def bond_moments(k, r0=1.0):
    r = np.linspace(max(0, r0 - 12 / np.sqrt(k)), r0 + 12 / np.sqrt(k), 20001)
    w = r**2 * np.exp(-0.5 * k * (r - r0)**2)
    return np.trapezoid(r**2 * w, r) / np.trapezoid(w, r), (np.trapezoid(r * w, r) / np.trapezoid(w, r))**2

def tau_int(x):
    x = np.asarray(x, float); x = x - x.mean(); n = len(x)
    if n < 16 or x.std() == 0: return np.nan
    f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); ac /= ac[0]
    t = 0.5
    for M in range(1, n // 2):
        t += ac[M]
        if M > 6 * t: break
    return t

def exact_model(N, J0, J1, J2, kappa, kb, kT=1.0, smax=20, two_state=False, theta0=np.pi):
    """Exact TM for the simulated model. Spins live on beads 0..N-2 (bond spins), bead N-1 dangles.
    two_state=True removes the coil state (n_states = 2 runs)."""
    beta = 1 / kT
    hr, cr = hinge_integrals(kappa, theta0)
    def E(a, b):
        if a == 0 and b == 0: return 0.0
        if a == 0 or b == 0:  return J1
        return -J0 if a == b else J2
    T0 = np.array([[np.exp(-beta * E(a, b)) for b in ST] for a in ST])       # plain pair weight (last bond)
    if two_state: T0[0, :] = 0.0; T0[:, 0] = 0.0
    Th = T0.copy()
    for a in (1, 2): Th[a, a] *= hr                                            # HH pairs carry a stiff hinge
    # products; coupling n (between beads n, n+1) uses Th for n = 0..N-3 and T0 for n = N-2
    Ts = [Th] * (N - 2) + [T0]
    L = [np.ones(3)]
    for n in range(N - 1): L.append(L[-1] @ Ts[n])
    R = [np.ones(3)]
    for n in range(N - 2, -1, -1): R.append(Ts[n] @ R[-1])
    R = R[::-1]                                # R[n] = product of Ts[n:] applied to ones
    Z = L[N - 1] @ np.ones(3)
    M = N - 1
    prob = lambda i: L[i] * R[i] / Z          # marginal of bead i
    hel = np.mean([prob(i)[1:].sum() for i in range(M)])
    coil = 1 - hel
    S = np.diag(ST)
    # C(s) over bond spins 0..M-1
    C = []
    for s in range(smax + 1):
        vals = []
        for i in range(M - s):
            v = L[i] @ S
            for n in range(i, i + s): v = v @ Ts[n]
            vals.append(v @ S @ R[i + s] / Z)
        C.append(np.mean(vals))
    C = np.array(C)
    # <N_d>: RL walls among bond spins (couplings 0..M-2)
    Nd = 0.0
    for n in range(M - 1):
        Trl = np.zeros((3, 3)); Trl[1, 2] = Ts[n][1, 2]; Trl[2, 1] = Ts[n][2, 1]
        Nd += L[n] @ Trl @ R[n + 1] / Z
    # <S^2>/M with S = sum of bond spins: sum_i sum_j <s_i s_j>
    tot = 0.0
    for i in range(M):
        v = L[i] @ S; tot += v @ S @ R[i] / Z
        for j in range(i + 1, M):
            v = v @ Ts[j - 1]; tot += 2 * (v @ S @ R[j] / Z)
    S2M = tot / M
    # <Re^2>: sum_i <r^2> + 2 sum_{i<j} <r>^2 cr^(j-i) P(beads i..j all R or all L)
    b2, rm2 = bond_moments(kb)
    Re2 = M * b2
    for i in range(M):
        vR = np.zeros(3); vR[1] = L[i][1]; vL = np.zeros(3); vL[2] = L[i][2]
        for j in range(i + 1, M):
            vR = vR @ Ts[j - 1] * np.array([0, 1, 0]); vL = vL @ Ts[j - 1] * np.array([0, 0, 1])
            P = (vR @ R[j] + vL @ R[j]) / Z
            Re2 += 2 * rm2 * cr ** (j - i) * P
    return dict(hel=hel, coil=coil, C=C, Nd=Nd, S2M=S2M, Re2=Re2, hr=hr, cr=cr, b2=b2, J=J0 + 0.5 * np.log(hr))

def paper(M, J, smax=20):
    c = np.exp(2 * J)
    # eq 7 as printed: <Re^2>_id = M l0^2 /(1+c)^(M-1) * sum_{N=1}^{M} c^(M-N) C(M-1,N-1) (2M-(N-1))/(N+1)
    Re2_id = M / (1 + c) ** (M - 1) * sum(c ** (M - n) * comb(M - 1, n - 1) * (2 * M - (n - 1)) / (n + 1) for n in range(1, M + 1))
    Re2_0 = M * (1 + 2 * c)
    Nd = (M - 1) * np.exp(-J) / (np.exp(J) + np.exp(-J))
    xi = 1 / abs(np.log(np.tanh(J)))
    C = np.exp(-np.arange(smax + 1) / xi)
    t = np.tanh(J)
    S2M = (1 + t) / (1 - t) - 2 * t / M * (1 - t ** M) / (1 - t) ** 2
    return dict(Re2_id=Re2_id, Re2_0=Re2_0, Nd=Nd, C=C, S2M=S2M, xi=xi)

# ---------------- collect runs, grouped by (J, kappa) ----------------
groups = {}
for f in sorted(glob.glob(f"{args.data}/*_obs.dat")):
    tag = os.path.basename(f)[:-8]
    m = re.match(r"J([\d.]+)_k(\d+)_s(\d+)", tag)
    groups.setdefault((float(m.group(1)), int(m.group(2))), []).append(tag)

smax = 20
report, curves = [], {}
for (Jt, kappa), tags in sorted(groups.items()):
    inp = read_input(f"{args.inputs}/{tags[0]}.dat"); N = int(inp["N"]); M = N - 1
    J0, J1, J2, kb = (float(inp[k]) for k in ("J0", "J1", "J2", "k_bond_CC"))
    ex = exact_model(N, J0, J1, J2, kappa, kb, smax=smax, two_state=inp.get("n_states", "3") == "2",
                     theta0=float(inp.get("theta0_HHH", np.pi))); pp = paper(M, ex["J"], smax)
    per = []
    for tag in tags:
        o = read_obs(f"{args.data}/{tag}_obs.dat"); S = read_frames(f"{args.data}/{tag}_conf.xyz")[:, :M]   # bond spins only
        b = int(args.burn * len(S)); S = S[b:]; bo = int(args.burn * len(o["sweep"])); o = {k: v[bo:] for k, v in o.items()}
        Re2 = o["Ree"] ** 2
        C = np.array([np.mean(S[:, :M - s] * S[:, s:]) for s in range(smax + 1)])
        Nd = np.mean(np.sum(S[:, :-1] * S[:, 1:] < 0, axis=1))
        Ssum = S.sum(1); S2M = np.mean(Ssum ** 2) / M
        per.append(dict(Re2=Re2.mean(), tau_Re2=tau_int(Re2) * (o["sweep"][1] - o["sweep"][0]), C=C, Nd=Nd, S2M=S2M,
                        hel=o["helicity"].mean(), tau_S=tau_int(Ssum) * 1000, nsw=o["sweep"][-1], acc=o["acc_pos"][-1]))
    n = len(per); mean = lambda k: np.mean([p[k] for p in per], axis=0); sem = lambda k: np.std([p[k] for p in per], axis=0, ddof=1) / np.sqrt(n)
    curves[(Jt, kappa)] = (mean("C"), sem("C"), ex["C"], pp["C"])
    report.append((Jt, kappa, n, ex, pp, mean, sem, per))

for Jt, kappa, n, ex, pp, mean, sem, per in report:
    print(f"\n=== target J = {Jt}, kappa = {kappa}: J0 = {read_input(f'{args.inputs}/J{Jt}_k{kappa}_s1.dat')['J0']}, "
          f"w_r/w_f = {ex['hr']:.3e}, <cos theta>_HH = {ex['cr']:.4f}, J_eff = {ex['J']:.4f}   ({n} seeds, "
          f"{per[0]['nsw']:.0e} sweeps, acc_pos = {per[0]['acc']:.2f})")
    print(f"  {'observable':<12} {'MC':>12} {'+-':>9} {'exact model':>12} {'paper (J)':>12}")
    print(f"  {'helicity':<12} {mean('hel'):>12.5f} {sem('hel'):>9.5f} {ex['hel']:>12.5f} {1.0:>12.5f}")
    print(f"  {'<Re^2>':<12} {mean('Re2'):>12.2f} {sem('Re2'):>9.2f} {ex['Re2']:>12.2f} {pp['Re2_id']:>12.2f}   (paper eq 8, M(1+2c): {pp['Re2_0']:.2f})")
    print(f"  {'<N_d>':<12} {mean('Nd'):>12.3f} {sem('Nd'):>9.3f} {ex['Nd']:>12.3f} {pp['Nd']:>12.3f}")
    print(f"  {'<S^2>/M':<12} {mean('S2M'):>12.3f} {sem('S2M'):>9.3f} {ex['S2M']:>12.3f} {pp['S2M']:>12.3f}")
    Cm, Ce = mean("C"), sem("C")
    print(f"  C(s):  s      MC        +-     exact   paper exp(-s/xi), xi={pp['xi']:.3f}")
    for s in (1, 2, 4, 6, 8, 10):
        print(f"        {s:>2} {Cm[s]:>9.4f} {Ce[s]:>8.4f} {ex['C'][s]:>9.4f} {pp['C'][s]:>9.4f}")
    print(f"  tau_int(Re^2) ~ {np.mean([p['tau_Re2'] for p in per]):.3g} sweeps, tau_int(S) ~ {np.mean([p['tau_S'] for p in per]):.3g} sweeps  (run = {per[0]['nsw']:.0e})")

# ---------------- figure: C(s) and Re^2 summary ----------------
BLUE, ORANGE, GREEN, YELLOW, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
cols = [BLUE, ORANGE, GREEN, YELLOW]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5)
ax = axs[0]; smax_plot = 12; s = np.arange(smax_plot + 1)
for k, ((Jt, kappa), (Cm, Ce, Cx, Cp)) in enumerate(sorted(curves.items())):
    ax.errorbar(s + 0.08 * k, Cm[:smax_plot + 1], yerr=Ce[:smax_plot + 1], fmt="o", ms=4, color=cols[k % 4], capsize=2, label=f"MC  J={Jt}, κ={kappa}")
    ax.plot(s, Cx[:smax_plot + 1], "-", lw=1.2, color=cols[k % 4], alpha=0.8)
    if k == 0: ax.plot(s, Cp[:smax_plot + 1], "--", lw=1.4, color=INK2, label="paper eq 13, exp(−s/ξ)")
ax.set_yscale("log"); ax.set_ylim(1e-4, 1.3); ax.set_xlabel("s  (bond separation)", color=INK); ax.set_ylabel("C(s) = ⟨σᵢ σᵢ₊ₛ⟩", color=INK)
ax.set_title("Spin–spin correlation, ideal chain with hinges", loc="left", fontsize=10, color=INK)
ax.text(0.02, 0.04, "solid: exact 3-state model (finite κ, coil)", transform=ax.transAxes, fontsize=8, color=INK2)
ax.legend(frameon=False, fontsize=8)
ax = axs[1]
labels, y, e, yx, yp = [], [], [], [], []
for Jt, kappa, n, ex, pp, mean, sem, per in report:
    labels.append(f"J={Jt}\nκ={kappa}"); y.append(mean("Re2")); e.append(sem("Re2")); yx.append(ex["Re2"]); yp.append(pp["Re2_id"])
x = np.arange(len(labels)); w = 0.26
ax.bar(x - w, y, w, yerr=e, color=BLUE, capsize=3, label="MC")
ax.bar(x, yx, w, color=INK2, label="exact model (finite κ, coil)")
ax.bar(x + w, yp, w, color=ORANGE, label="paper eq 7 (rigid rods)")
ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylabel("⟨R_e²⟩  [l₀²]", color=INK)
ax.set_title("Mean squared end-to-end distance (M = 100)", loc="left", fontsize=10, color=INK)
ax.set_ylim(0, max(max(y), max(yp)) * 1.25)
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"paper_ideal{args.tag}.png")
