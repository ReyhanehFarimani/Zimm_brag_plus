#!/usr/bin/env python3
"""J-scan of the no-steric / no-bend chain vs the exact solution.

With no bending, bond directions are independent, but the SOFT state-dependent bonds couple to the spins:
integrating the bond length of pair class c out gives a weight Z_c = int_0^inf r^2 exp(-k_c (r-r0_c)^2 / 2kT) dr.
Spin sector: exact 3-state transfer matrix with T[a][b] = exp(-E_pair(a,b)/kT) * Z_class(a,b).
Geometry:    <Re^2> = sum_bonds <r^2>_class, weighted by the exact per-bond class probabilities.
"""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))

ST = [0, 1, -1]                    # C, R, L

def read_input(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
    return d

def read_obs(f):
    lines = open(f).read().splitlines()
    hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}

def pair_name(a, b):
    if a == 0 and b == 0: return "CC"
    if a == 0 or b == 0:  return "CH"
    return "HH" if a == b else "RL"

def bond_int(r0, k, kT=1.0):
    """Z = int r^2 e^{-k(r-r0)^2/2kT} dr  and  <r^2> under that weight."""
    sig = np.sqrt(kT / k)
    r = np.linspace(0.0, r0 + 14 * sig, 60001)
    w = r**2 * np.exp(-0.5 * k / kT * (r - r0)**2)
    Z = np.trapezoid(w, r)
    return Z, np.trapezoid(r**2 * w, r) / Z

def exact(inp, J):
    N = int(inp["N"]); kT = float(inp["kT"]); J1 = float(inp["J1"]); Eh = float(inp.get("E_helix", 0.0))
    E = {"CC": 0.0, "CH": J1, "HH": -J, "RL": J}
    site = np.array([np.exp(-Eh / kT) if a != 0 else 1.0 for a in ST])      # on-site weight of each state
    Zb, r2 = {}, {}
    for c in ("CC", "CH", "HH", "RL"):
        Zb[c], r2[c] = bond_int(float(inp[f"bond_len_{c}"]), float(inp[f"k_bond_{c}"]), kT)
    T = np.array([[np.exp(-E[pair_name(a, b)] / kT) * Zb[pair_name(a, b)] for b in ST] for a in ST]) * site[None, :]
    T /= np.abs(np.linalg.eigvals(T)).max()
    L = [site.copy()]                                                        # first bead carries its own site weight
    for _ in range(1, N): L.append(L[-1] @ T)
    R = [np.ones(3)]
    for _ in range(1, N): R.append(T @ R[-1])
    Z = L[N - 1] @ np.ones(3)
    hel = np.mean([(L[i] * R[N - 1 - i])[1:].sum() / Z for i in range(N)])
    # per-bond pair-class probabilities -> <E_state>, <Re^2>, and <S^2>/N (S = sum of spins)
    Es = Eh * hel * N; Re2 = 0.0
    for n in range(N - 1):
        for ia, a in enumerate(ST):
            for ib, b in enumerate(ST):
                p = L[n][ia] * T[ia, ib] * R[N - 2 - n][ib] / Z
                c = pair_name(a, b)
                Es += p * E[c]; Re2 += p * r2[c]
    Sdiag = np.diag(ST)
    tot = 0.0
    for i in range(N):
        v = L[i] @ Sdiag; tot += v @ Sdiag @ R[N - 1 - i] / Z
        for j in range(i + 1, N):
            v = v @ T; tot += 2 * (v @ Sdiag @ R[N - 1 - j] / Z)
    return dict(hel=hel, Es=Es, Re2=Re2, S2N=tot / N)

# ---------------- collect ----------------
inp0 = read_input("input_1.dat"); N = int(inp0["N"])
runs = {}
for f in sorted(glob.glob("out/J*_obs.dat")):
    J = float(re.search(r"J([\d.]+)_s(\d+)_obs", f).group(1))
    o = read_obs(f); n = len(o["helicity"])
    nb = 20; blk = lambda x: x[: n // nb * nb].reshape(nb, -1).mean(1)
    S = o["n_R"] - o["n_L"]
    runs.setdefault(J, []).append(dict(
        hel=(blk(o["helicity"]).mean(), blk(o["helicity"]).std(ddof=1) / np.sqrt(nb)),
        Es=(blk(o["E_state"]).mean(), blk(o["E_state"]).std(ddof=1) / np.sqrt(nb)),
        Re2=(blk(o["Ree"]**2).mean(), blk(o["Ree"]**2).std(ddof=1) / np.sqrt(nb)),
        S2N=(blk(S**2 / N).mean(), blk(S**2 / N).std(ddof=1) / np.sqrt(nb))))

Js = sorted(runs)
mc, err = {}, {}
for key in ("hel", "Es", "Re2", "S2N"):
    mc[key]  = np.array([np.mean([r[key][0] for r in runs[J]]) for J in Js])
    err[key] = np.array([np.sqrt(np.sum([r[key][1] ** 2 for r in runs[J]])) / len(runs[J]) for J in Js])
Jline = np.linspace(min(Js), max(Js), 60)
exl = [exact(inp0, J) for J in Jline]
exp_ = {k: np.array([e[k] for e in exl]) for k in ("hel", "Es", "Re2", "S2N")}
ext  = {k: np.array([exact(inp0, J)[k] for J in Js]) for k in ("hel", "Es", "Re2", "S2N")}

print(f"N = {N}, J1 = {inp0['J1']}, J0 = J2 = J, no bending, soft state-dependent bonds ({len(runs[Js[0]])} seeds)")
print(f"{'J':>5} {'helicity':>9} {'+-':>7} {'exact':>8} | {'E_state':>9} {'+-':>6} {'exact':>9} | {'<Re2>':>8} {'+-':>6} {'exact':>8} | {'<S2>/N':>8} {'+-':>7} {'exact':>8}")
for i, J in enumerate(Js):
    print(f"{J:>5} {mc['hel'][i]:>9.4f} {err['hel'][i]:>7.4f} {ext['hel'][i]:>8.4f} |"
          f" {mc['Es'][i]:>9.1f} {err['Es'][i]:>6.1f} {ext['Es'][i]:>9.1f} |"
          f" {mc['Re2'][i]:>8.1f} {err['Re2'][i]:>6.1f} {ext['Re2'][i]:>8.1f} |"
          f" {mc['S2N'][i]:>8.3f} {err['S2N'][i]:>7.3f} {ext['S2N'][i]:>8.3f}")

# ---------------- figure ----------------
BLUE, ORANGE, GREEN, INK, INK2, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.1), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5)
    ax.set_xlabel("J   (J0 = J2 = J)", color=INK)
ax = axs[0]
ax.plot(Jline, exp_["hel"], "-", lw=1.5, color=INK2, label="exact transfer matrix")
ax.errorbar(Js, mc["hel"], yerr=err["hel"], fmt="o", ms=5, color=BLUE, capsize=2.5, label="MC")
ax.set_ylabel("helicity  (fraction R or L)", color=INK); ax.set_ylim(0, 1.02)
ax.set_title(f"Helix content   (N={N}, J1={inp0['J1']}, kT=1)", loc="left", fontsize=9.5, color=INK)
ax.legend(frameon=False, fontsize=8.5)
ax = axs[1]
ax.plot(Jline, exp_["Re2"], "-", lw=1.5, color=INK2, label="exact")
ax.errorbar(Js, mc["Re2"], yerr=err["Re2"], fmt="o", ms=5, color=ORANGE, capsize=2.5, label="MC")
ax.set_ylabel(r"$\langle R_e^2\rangle$", color=INK)
ax.set_title("End-to-end distance  (soft state-dependent bonds)", loc="left", fontsize=9.5, color=INK)
ax.legend(frameon=False, fontsize=8.5)
ax = axs[2]
ax.plot(Jline, exp_["S2N"], "-", lw=1.5, color=INK2, label="exact")
ax.errorbar(Js, mc["S2N"], yerr=err["S2N"], fmt="o", ms=5, color=GREEN, capsize=2.5, label="MC")
ax.set_ylabel(r"$\langle S^2\rangle/N$,  $S=n_R-n_L$", color=INK)
ax.set_title("Chiral order fluctuations", loc="left", fontsize=9.5, color=INK)
ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("jscan.png")
