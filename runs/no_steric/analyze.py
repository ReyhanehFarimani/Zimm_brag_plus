#!/usr/bin/env python3
"""J-scan of the no-steric chain (soft state-dependent bonds + state-dependent bending, coil allowed).
MC observables vs J with blocking errors. No closed-form reference for the bent chain; the exact no-bend
solution (runs/no_steric_no_bend/analyze.py) is drawn dashed for comparison where it applies."""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here)
sys.argv = [sys.argv[0]]
src = open("../no_steric_no_bend/analyze.py").read().split("# ---------------- collect")[0]
exec(src)                                    # read_input, read_obs, exact(inp, J)  (no-bend reference)
_read_input = read_input
def read_input(f):
    d = _read_input(f)
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l and "deg" in l:
            d[l.split("=")[0].strip() + "_deg"] = True
    return d

inp0 = read_input("data_sample.dat"); N = int(inp0["N"])
gen = sorted(glob.glob("inputs/J*_s1.dat"))
inpg = read_input(gen[0]) if gen else inp0                 # a generated input (theta0 etc. as actually run)

# ---------------- exact solution WITH bending: transfer matrix on consecutive spin pairs ----------------
TRIPLE = {}
def triple_class(a, b, c):
    h = lambda x: x != 0
    nh = h(a) + h(b) + h(c)
    if nh == 0: return "CCC"
    if nh == 1: return "CHC" if h(b) else "CCH"
    if nh == 2:
        if not h(b): return "HCH" if a == c else "RCL"
        return "CHH" if (a if h(a) else c) == b else "CRL"
    if a == b == c: return "HHH"
    return "RLR" if a == c else "RRL"

def hinge_moments(kappa, theta0, kT=1.0):
    """W = 2pi int sin(t) e^{-U/kT} dt,  <U>,  and <cos(bond-bond angle)> = -<cos(valence angle)>"""
    t = np.linspace(0, np.pi, 200001); U = 0.5 * kappa * (t - theta0) ** 2; w = np.sin(t) * np.exp(-U / kT)
    Z = np.trapezoid(w, t)
    return 2 * np.pi * Z, np.trapezoid(U * w, t) / Z, -np.trapezoid(np.cos(t) * w, t) / Z

def theta_of(inp, key):
    v = inp[key]
    return float(v) * (np.pi / 180 if "deg" in inp.get(key + "_unit", "") else 1.0)

def exact_bend(inp, J, J2=None):
    """J = J0 (E_HH = -J); E_RL = +J2 (defaults to J, i.e. the J0 = J2 scans); E_CH = +J1 from the input."""
    N = int(inp["N"]); kT = float(inp["kT"]); J1 = float(inp["J1"]); Eh = float(inp.get("E_helix", 0.0))
    E = {"CC": 0.0, "CH": J1, "HH": -J, "RL": (J if J2 is None else J2)}
    sw = lambda a: np.exp(-Eh / kT) if a != 0 else 1.0                      # on-site weight
    Zb, r2, rm = {}, {}, {}
    for c in ("CC", "CH", "HH", "RL"):
        r0, k = float(inp[f"bond_len_{c}"]), float(inp[f"k_bond_{c}"])
        sig = np.sqrt(kT / k); r = np.linspace(0.0, r0 + 14 * sig, 60001); w = r**2 * np.exp(-0.5 * k / kT * (r - r0)**2)
        Zb[c] = np.trapezoid(w, r); r2[c] = np.trapezoid(r**2 * w, r) / Zb[c]; rm[c] = np.trapezoid(r * w, r) / Zb[c]
    hw = {}
    for t_ in ("CCC", "CCH", "CHC", "CHH", "CRL", "HCH", "RCL", "HHH", "RRL", "RLR"):
        th0 = float(inp[f"theta0_{t_}"]); th0 = np.deg2rad(th0) if inp.get(f"theta0_{t_}_deg") else th0
        hw[t_] = hinge_moments(float(inp[f"kappa_{t_}"]), th0, kT)
    idx = {(a, b): 3 * i + j for i, a in enumerate(ST) for j, b in enumerate(ST)}
    B = {(a, b): np.exp(-E[pair_name(a, b)] / kT) * Zb[pair_name(a, b)] for (a, b) in idx}
    T, Tc, Ub = np.zeros((9, 9)), np.zeros((9, 9)), np.zeros((9, 9))
    for (a, b), p in idx.items():
        for (b2, c), q in idx.items():
            if b2 != b: continue
            W, Umean, cosab = hw[triple_class(a, b, c)]
            T[p, q] = W * B[(b, c)] * sw(c); Tc[p, q] = T[p, q] * cosab; Ub[p, q] = T[p, q] * Umean
    scale = np.abs(np.linalg.eigvals(T)).max(); T /= scale; Tc /= scale; Ub /= scale
    v0 = np.zeros(9)
    for (a, b), p in idx.items(): v0[p] = B[(a, b)] * sw(a) * sw(b)
    L = [v0]
    for _ in range(N - 2): L.append(L[-1] @ T)
    R = [np.ones(9)]
    for _ in range(N - 2): R.append(T @ R[-1])
    R = R[::-1]                                   # R[k] = T^{N-2-k} 1
    Z = L[N - 2] @ np.ones(9)
    first = np.array([a for (a, b) in sorted(idx, key=idx.get)]); second = np.array([b for (a, b) in sorted(idx, key=idx.get)])
    # marginals of pairs k = 0..N-2
    P = [L[k] * R[k] / Z for k in range(N - 1)]
    hel = (sum((P[k] * (first != 0)).sum() for k in range(N - 1)) + (P[N - 2] * (second != 0)).sum()) / N
    Es = sum(sum(P[k][idx[(a, b)]] * E[pair_name(a, b)] for (a, b) in idx) for k in range(N - 1)) + Eh * hel * N
    Eb = sum(L[j - 1] @ Ub @ R[j] / Z for j in range(1, N - 1))
    # <Re^2> = sum_k <r_k^2> + 2 sum_{i<j} <r_i><r_j> prod cos
    D2 = np.array([r2[pair_name(a, b)] for (a, b) in sorted(idx, key=idx.get)])
    Dm = np.array([rm[pair_name(a, b)] for (a, b) in sorted(idx, key=idx.get)])
    Re2 = sum((P[k] * D2).sum() for k in range(N - 1))
    for i in range(N - 1):
        w = L[i] * Dm
        for j in range(i + 1, N - 1):
            w = w @ Tc; Re2 += 2 * (w * Dm) @ R[j] / Z
    # <S^2>/N, S = sum of all N spins (bead k = first of pair k; bead N-1 = second of pair N-2)
    tot = 0.0
    for k in range(N - 1):
        tot += (P[k] * first**2).sum()
        w = L[k] * first
        for l in range(k + 1, N - 1):
            w = w @ T; tot += 2 * (w * first) @ R[l] / Z
        tot += 2 * (L[k] * first) @ np.linalg.matrix_power(T, N - 2 - k) @ (second * R[N - 2]) / Z if k < N - 2 else 2 * (L[k] * first * second) @ R[k] / Z
    tot += (P[N - 2] * second**2).sum()
    return dict(hel=hel, Es=Es, Eb=Eb, Re2=Re2, S2N=tot / N)
runs = {}
for f in sorted(glob.glob("out/J*_obs.dat")):
    J = float(re.search(r"J([\d.]+)_s(\d+)_obs", f).group(1))
    o = read_obs(f); n = len(o["helicity"])
    if n < 40: continue
    nb = 20; blk = lambda x: x[: n // nb * nb].reshape(nb, -1).mean(1)
    S = o["n_R"] - o["n_L"]
    obs = dict(hel=o["helicity"], Es=o["E_state"], Eb=o["E_bend"], Re2=o["Ree"] ** 2, Rg2=o["Rg2"], S2N=S ** 2 / N)
    runs.setdefault(J, []).append({k: (blk(v).mean(), blk(v).std(ddof=1) / np.sqrt(nb)) for k, v in obs.items()})

Js = sorted(runs); keys = ("hel", "Es", "Eb", "Re2", "Rg2", "S2N")
mc  = {k: np.array([np.mean([r[k][0] for r in runs[J]]) for J in Js]) for k in keys}
err = {k: np.array([np.sqrt(np.sum([r[k][1] ** 2 for r in runs[J]])) / len(runs[J]) for J in Js]) for k in keys}
Jline = np.linspace(min(Js), max(Js), 50)
nobend = [exact(inp0, J) for J in Jline]           # same bonds / J's, kappa = 0
withb  = [exact_bend(inpg, J) for J in Jline]      # exact solution of the simulated (bent) model
exb    = {J: exact_bend(inpg, J) for J in Js}

print(f"N = {N}, J1 = {inpg['J1']}, J0 = J2 = J, theta0_HHH = {inpg['theta0_HHH']}{' deg' if inpg.get('theta0_HHH_deg') else ''}, {len(runs[Js[0]])} seeds")
print("exact = transfer matrix on spin pairs with bond and hinge partition functions folded in (the simulated model)")
print(f"{'J':>5} {'helicity':>9} {'+-':>7} {'exact':>8} | {'E_state':>8} {'exact':>8} {'E_bend':>7} {'exact':>7} | {'<Re2>':>7} {'+-':>5} {'exact':>7} | {'<S2>/N':>8} {'+-':>7} {'exact':>8}")
for i, J in enumerate(Js):
    e = exb[J]
    print(f"{J:>5} {mc['hel'][i]:>9.4f} {err['hel'][i]:>7.4f} {e['hel']:>8.4f} | {mc['Es'][i]:>8.1f} {e['Es']:>8.1f} {mc['Eb'][i]:>7.1f} {e['Eb']:>7.1f} |"
          f" {mc['Re2'][i]:>7.0f} {err['Re2'][i]:>5.0f} {e['Re2']:>7.0f} | {mc['S2N'][i]:>8.3f} {err['S2N'][i]:>7.3f} {e['S2N']:>8.3f}")

BLUE, ORANGE, GREEN, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.1), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.set_xlabel("J   (J0 = J2 = J)", color=INK)
ax = axs[0]
ax.plot(Jline, [e["hel"] for e in nobend], "--", lw=1.3, color=INK2, label="exact, same bonds, no bending")
ax.plot(Jline, [e["hel"] for e in withb], "-", lw=1.5, color=INK2, label="exact, with bending")
ax.errorbar(Js, mc["hel"], yerr=err["hel"], fmt="o", ms=5, color=BLUE, capsize=2.5, label="MC")
ax.set_ylabel("helicity  (fraction R or L)", color=INK); ax.set_ylim(0, 1.02)
ax.set_title(f"Helix content   (N={N}, J1={inpg['J1']}, θ0=π)", loc="left", fontsize=9.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
ax = axs[1]
ax.plot(Jline, [e["Re2"] for e in nobend], "--", lw=1.3, color=INK2, label="no bending (exact)")
ax.plot(Jline, [e["Re2"] for e in withb], "-", lw=1.5, color=INK2, label="with bending (exact ⟨Rₑ²⟩)")
ax.errorbar(Js, mc["Re2"], yerr=err["Re2"], fmt="o", ms=5, color=ORANGE, capsize=2.5, label="MC ⟨Rₑ²⟩")
ax.errorbar(Js, 6 * mc["Rg2"], yerr=6 * err["Rg2"], fmt="s", ms=4.5, mfc=BG, color=ORANGE, capsize=2.5, label="MC 6⟨Rg²⟩")
ax.set_ylabel(r"$\langle R_e^2\rangle$,  $6\langle R_g^2\rangle$", color=INK)
ax.set_title("Chain size", loc="left", fontsize=9.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
ax = axs[2]
ax.plot(Jline, [e["S2N"] for e in nobend], "--", lw=1.3, color=INK2, label="no bending (exact)")
ax.plot(Jline, [e["S2N"] for e in withb], "-", lw=1.5, color=INK2, label="with bending (exact)")
ax.errorbar(Js, mc["S2N"], yerr=err["S2N"], fmt="o", ms=5, color=GREEN, capsize=2.5, label="MC")
ax.set_yscale("log"); ax.set_ylabel(r"$\langle S^2\rangle/N$,  $S=n_R-n_L$", color=INK)
ax.set_title("Chiral order fluctuations", loc="left", fontsize=9.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("jscan.png")
