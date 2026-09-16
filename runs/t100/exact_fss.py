#!/usr/bin/env python3
"""Exact (no-steric) finite-size behaviour of the t100 model vs E_helix: helicity for several N, the N->inf limit,
the correlation length xi, and the helicity susceptibility chi = N (<h^2> - <h>^2)."""
import numpy as np, sys
ST = [0, 1, -1]
def read_input(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
            if "deg" in v: d[k + "_deg"] = True
    return d
def pair_name(a, b):
    if a == 0 and b == 0: return "CC"
    if a == 0 or b == 0:  return "CH"
    return "HH" if a == b else "RL"
def triple_class(a, b, c):
    h = lambda x: x != 0; nh = h(a) + h(b) + h(c)
    if nh == 0: return "CCC"
    if nh == 1: return "CHC" if h(b) else "CCH"
    if nh == 2:
        if not h(b): return "HCH" if a == c else "RCL"
        return "CHH" if (a if h(a) else c) == b else "CRL"
    if a == b == c: return "HHH"
    return "RLR" if a == c else "RRL"
def hinge_W(kappa, theta0):
    t = np.linspace(0, np.pi, 200001); return 2 * np.pi * np.trapezoid(np.sin(t) * np.exp(-0.5 * kappa * (t - theta0)**2), t)
inp = read_input(sys.argv[1] if len(sys.argv) > 1 else "sample_data.dat")
J0, J1, J2 = (float(inp[k]) for k in ("J0", "J1", "J2"))
Zb = {}
for c in ("CC", "CH", "HH", "RL"):
    r0, k = float(inp[f"bond_len_{c}"]), float(inp[f"k_bond_{c}"]); r = np.linspace(0, r0 + 14 / np.sqrt(k), 60001)
    Zb[c] = np.trapezoid(r**2 * np.exp(-0.5 * k * (r - r0)**2), r)
hw = {}
for t_ in ("CCC", "CCH", "CHC", "CHH", "CRL", "HCH", "RCL", "HHH", "RRL", "RLR"):
    th0 = float(inp[f"theta0_{t_}"]); th0 = np.deg2rad(th0) if inp.get(f"theta0_{t_}_deg") else th0
    hw[t_] = hinge_W(float(inp[f"kappa_{t_}"]), th0)
idx = {(a, b): 3 * i + j for i, a in enumerate(ST) for j, b in enumerate(ST)}
first = np.array([a for (a, b) in sorted(idx, key=idx.get)])
def build(Eh):
    E = {"CC": 0.0, "CH": J1, "HH": -J0, "RL": J2}; sw = lambda a: np.exp(-Eh) if a != 0 else 1.0
    T = np.zeros((9, 9))
    for (a, b), p in idx.items():
        for (b2, c), q in idx.items():
            if b2 == b: T[p, q] = np.exp(-E[pair_name(b, c)]) * Zb[pair_name(b, c)] * sw(c) * hw[triple_class(a, b, c)]
    v0 = np.array([np.exp(-E[pair_name(a, b)]) * Zb[pair_name(a, b)] * sw(a) * sw(b) for (a, b) in sorted(idx, key=idx.get)])
    return T, v0
def hel_N(Eh, N, moments=False):
    T, v0 = build(Eh); T = T / np.abs(np.linalg.eigvals(T)).max()
    L = [v0]
    for _ in range(N - 2): L.append(L[-1] @ T)
    R = [np.ones(9)]
    for _ in range(N - 2): R.append(T @ R[-1])
    R = R[::-1]; Z = L[N - 2] @ np.ones(9)
    H = (first != 0).astype(float)
    h1 = sum((L[k] * R[k] / Z * H).sum() for k in range(N - 1))
    if not moments: return h1 / N
    # <n_h^2> = sum_k sum_l <h_k h_l>  (over beads 0..N-2)
    tot = 0.0
    for k in range(N - 1):
        w = L[k] * H; tot += (w * R[k]).sum() / Z
        for l in range(k + 1, N - 1):
            w = w @ T; tot += 2 * (w * H * R[l]).sum() / Z
    return h1 / N, (tot - h1 * h1) / N            # helicity, chi = N var(h)
def hel_inf(Eh):
    T, _ = build(Eh); w, v = np.linalg.eig(T); i = np.argmax(np.abs(w)); wl, vl = np.linalg.eig(T.T); j = np.argmax(np.abs(wl))
    p = np.abs(v[:, i]) * np.abs(vl[:, j]); p /= p.sum(); return p[first != 0].sum()
def xi(Eh):
    T, _ = build(Eh); lam = np.sort(np.abs(np.linalg.eigvals(T)))[::-1]; return 1 / np.log(lam[0] / lam[1])
if __name__ == "__main__":
    print(f"exact, no sterics: E_HH={-J0:+g} E_CH={J1:+g} E_RL={J2:+g}")
    print(f"{'E_helix':>8} {'N=50':>7} {'N=100':>7} {'N=200':>7} {'N=400':>7} {'N=inf':>7} {'xi':>8} | chi=N var(h): {'N=50':>6} {'N=100':>6} {'N=200':>6} {'N=400':>6}")
    for Eh in (-11.5, -11.2, -11.0, -10.9, -10.8, -10.7, -10.6, -10.5, -10.4, -10.2, -10.0):
        hs = [hel_N(Eh, N, True) for N in (50, 100, 200, 400)]
        print(f"{Eh:>8.1f} " + " ".join(f"{h:7.4f}" for h, c in hs) + f" {hel_inf(Eh):7.4f} {xi(Eh):8.2f} |               " + " ".join(f"{c:6.2f}" for h, c in hs))
    Es = np.linspace(-11.5, -10.0, 151); h = np.array([hel_inf(E) for E in Es]); dh = -np.gradient(h, Es); k = np.argmax(dh)
    print(f"\nN -> inf: midpoint E* = {Es[np.argmin(abs(h - 0.5))]:.3f}, max slope dh/dE = {dh[k]:.2f} per kT  (finite => crossover, not a phase transition)")
    for N in (50, 100, 200, 400, 800):
        hN = np.array([hel_N(E, N) for E in Es]); d = -np.gradient(hN, Es); print(f"  N = {N:4d}: midpoint {Es[np.argmin(abs(hN - 0.5))]:.3f}, max |slope| {d.max():.2f} per kT, width {1/d.max():.3f} kT")
