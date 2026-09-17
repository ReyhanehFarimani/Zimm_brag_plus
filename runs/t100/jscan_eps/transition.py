#!/usr/bin/env python3
"""jscan_eps: is there a phase transition in the handedness m = (n_R - n_L)/N, and does eps_s affect it?

EXACT reference: the same chain WITHOUT non-bonded interactions is a 1D nearest-neighbour model; its full
distribution P(M) is computed here by a polynomial transfer matrix (weights from ../exact_fss.py: state energies,
bond partition functions, hinge weights).  From P(M): <|m|>, <m^2>, Binder cumulant U4 = 1 - <m^4>/(3<m^2>^2)
(0 = disordered Gaussian, 2/3 = two sharp peaks at +-m), susceptibility chi = N(<m^2> - <|m|>^2), and the
R<->L-odd correlation length xi_m.  A true transition needs the U4(J) curves of different N to CROSS at a fixed
J_c and chi_max to grow as a power of N; a 1D crossover has no crossing and J*(N) drifts ~ ln N.

MC (N = 200, sterics + fitted helix-helix potential, every eps_s): same quantities, finished runs only, error =
max(seed scatter, combined block error).  z = (MC - exact)/err measures the total effect of ALL non-bonded
interactions; the spread between eps_s columns measures the chiral sector."""
import glob, importlib.util, os, re, sys
import numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["x", "../sample_data.dat"]
spec = importlib.util.spec_from_file_location("ex", os.path.abspath("../exact_fss.py")); ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)
PAIRS = sorted(ex.idx, key=ex.idx.get); SECOND = np.array([b for (a, b) in PAIRS]); FIRST = np.array([a for (a, b) in PAIRS])

def set_J(J): ex.J0, ex.J1, ex.J2 = float(J), 0.0, float(J)          # E_HH = -J, E_CH = 0, E_RL = +J, E_helix = 0

def exact_PM(J, N):
    """exact P(M), M = sum of spins = n_R - n_L, M = -N..N"""
    set_J(J); T, v0 = ex.build(0.0)
    P = np.zeros((9, 2 * N + 1))
    for p, (a, b) in enumerate(PAIRS): P[p, N + a + b] = v0[p]
    for _ in range(N - 2):
        Q = T.T @ P                                                  # Q[q] = sum_p T[p, q] P[p]
        for q in range(9): Q[q] = np.roll(Q[q], SECOND[q])           # the new bead adds its spin
        P = Q / Q.sum()
    pm = P.sum(0); return pm / pm.sum()

def exact_hel(J, N): set_J(J); return ex.hel_N(0.0, N)

def xi_m(J):
    """correlation length (in residues) of the R<->L-odd sector"""
    set_J(J); T, _ = ex.build(0.0)
    sw = {0: 0, 1: -1, -1: 1}; S = np.zeros((9, 9))
    for p, (a, b) in enumerate(PAIRS): S[ex.idx[(sw[a], sw[b])], p] = 1.0
    A = 0.5 * (np.eye(9) - S); lam0 = np.abs(np.linalg.eigvals(T)).max(); lam1 = np.abs(np.linalg.eigvals(A @ T @ A)).max()
    return 1.0 / np.log(lam0 / lam1)

def moments(pm, N):
    m = np.arange(-N, N + 1) / N; m2 = (pm * m**2).sum(); m4 = (pm * m**4).sum(); am = (pm * np.abs(m)).sum()
    return am, m2, 1 - m4 / (3 * m2 * m2), N * (m2 - am * am)

# ---------------------------------------------------------------- exact: N dependence
Jg = np.arange(1.0, 8.01, 0.25); NS = (100, 200, 400, 800, 1600)
ex_tab = {N: np.array([moments(exact_PM(J, N), N) for J in Jg]) for N in NS}
print("EXACT 1D chain (no non-bonded interactions), E_HH = -J, E_CH = 0, E_RL = +J, E_helix = 0")
print("\n  <|m|> vs J and N;   xi_m = R/L correlation length in residues")
print("     J     xi_m " + "".join(f"   N={N:<5d}" for N in NS))
for i, J in enumerate(Jg):
    if J % 0.5 == 0: print(f"  {J:4.2f} {xi_m(J):8.1f} " + "".join(f"  {ex_tab[N][i, 0]:8.4f}" for N in NS))
print("\n  Binder cumulant U4 vs J and N  (a transition would make these columns CROSS at one J)")
print("     J " + "".join(f"   N={N:<5d}" for N in NS))
for i, J in enumerate(Jg):
    if J % 0.5 == 0: print(f"  {J:4.2f} " + "".join(f"  {ex_tab[N][i, 2]:8.4f}" for N in NS))
print("\n  crossover location and sharpness vs N")
print("      N   J*(<|m|> = 0.5)   J(chi max)    chi_max    chi_max/N   max d<|m|>/dJ")
Jf = np.linspace(Jg[0], Jg[-1], 2801); prev = None
for N in NS:
    am = np.interp(Jf, Jg, ex_tab[N][:, 0]); chi = ex_tab[N][:, 3]; k = chi.argmax()
    js = Jf[np.argmin(np.abs(am - 0.5))]
    print(f"  {N:5d}   {js:15.3f}   {Jg[k]:10.2f}   {chi[k]:8.2f}   {chi[k] / N:9.4f}   {np.gradient(ex_tab[N][:, 0], Jg).max():13.3f}"
          + (f"     (J* shift for doubling N: {js - prev:+.3f})" if prev is not None else ""))
    prev = js

# ---------------------------------------------------------------- MC vs exact at N = 200
N = 200
def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def jack(m, nb=20):
    """block jackknife of (<|m|>, <m^2>, U4, chi)"""
    f = lambda x: np.array([np.abs(x).mean(), (x**2).mean(), 1 - (x**4).mean() / (3 * (x**2).mean()**2), N * ((x**2).mean() - np.abs(x).mean()**2)])
    B = m[: len(m) // nb * nb].reshape(nb, -1); full = f(m); jk = np.array([f(np.delete(B, i, 0).ravel()) for i in range(nb)])
    return full, np.sqrt((nb - 1) / nb * ((jk - jk.mean(0))**2).sum(0))
mc, hist = {}, {}
for log in sorted(glob.glob("logs/J*_es*_s*.log")):
    if "summary" not in open(log).read(): continue
    J, es, s = map(int, re.search(r"J(\d+)_es(\d+)_s(\d+)", log).groups())
    o = read_obs("out/" + os.path.basename(log)[:-4] + "_obs.dat"); m = (o["n_R"] - o["n_L"]) / N
    mc.setdefault((J, es), []).append(jack(m)); hist.setdefault(J, []).append(m)
Js = sorted({k[0] for k in mc}); ESs = sorted({k[1] for k in mc}); names = ["<|m|>", "<m^2>", "U4", "chi"]
ex200 = {J: moments(exact_PM(J, N), N) for J in Js}
print(f"\nMC (N = 200) vs EXACT 1D:  value [z = (MC - exact)/err];  {sum(len(v) for v in mc.values())} finished runs")
allz = {n: [] for n in names}; spread = {n: [] for n in names}
for q, name in enumerate(names):
    print(f"\n  {name}\n     J      exact " + "".join(f"{'eps_s=' + str(e):>19s}" for e in ESs))
    for J in Js:
        row = f"  {J:4d} {ex200[J][q]:10.4f} "; vals, errs = [], []
        for es in ESs:
            if (J, es) not in mc: row += f"{'--':>19s}"; continue
            v = np.array([r[0][q] for r in mc[(J, es)]]); e = np.array([r[1][q] for r in mc[(J, es)]])
            err = max(v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0, np.sqrt((e**2).sum()) / len(v))
            z = (v.mean() - ex200[J][q]) / err; allz[name].append(z); vals.append(v.mean()); errs.append(err)
            row += f"{v.mean():10.4f} [{z:+5.1f}]"
        vals, errs = np.array(vals), np.array(errs); w = 1 / errs**2; mu = (w * vals).sum() / w.sum()
        spread[name].append(((vals - mu)**2 * w).sum() / (len(vals) - 1))          # chi2/dof of "no eps_s dependence"
        print(row)
print("\n  summary over all cells:")
for n in names:
    z = np.array(allz[n]); print(f"    {n:6s} MC vs exact 1D: RMS z = {np.sqrt((z**2).mean()):5.2f}, mean z = {z.mean():+5.2f}   |   "
                                 f"eps_s dependence at fixed J: chi2/dof = {np.mean(spread[n]):4.2f} (1 = none)")
print("\n  shape of P(m) in the MC (all eps_s pooled): fraction of samples with |m| < 0.1 / 0.1-0.9 / > 0.9")
for J in Js:
    m = np.abs(np.concatenate(hist[J])); pe = exact_PM(J, N); me = np.abs(np.arange(-N, N + 1) / N)
    print(f"    J = {J}:  MC {np.mean(m < 0.1):.3f} / {np.mean((m >= 0.1) & (m <= 0.9)):.3f} / {np.mean(m > 0.9):.3f}"
          f"     exact {pe[me < 0.1].sum():.3f} / {pe[(me >= 0.1) & (me <= 0.9)].sum():.3f} / {pe[me > 0.9].sum():.3f}")
np.savez("transition_exact.npz", Jg=Jg, NS=np.array(NS), **{f"N{N}": ex_tab[N] for N in NS})
