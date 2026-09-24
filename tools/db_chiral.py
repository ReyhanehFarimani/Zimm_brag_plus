#!/usr/bin/env python3
"""RR vs RL of the tabulated helix-helix potential (user 2026-09-23: "extract the RR vs RL difference from the table"):
Boltzmann angular averages per r (isotropic rod orientations, free registries), the parallel side-by-side contact
(betaA = betaB = 90, psi = 0) with registries free / optimal, association constants and the crossing-angle profile.
  ~/venv/bin/python tools/db_chiral.py [helix_pair_db_t45_win.bin]      (numpy only; memmaps the 853 MB table, ~1 min)
Result (t45 and t100 win tables): in the reliable region r >= 1.9 a the same- and opposite-handed contacts are equal to
< 0.02 kT with free registries; the only handedness dependence sits inside r < 1.6 a in the parallel side-by-side
geometry with optimal registries, where both are repulsive (t45: RR softer by 0.3 kT at 1.5 a, 2.5 kT at 1.2 a;
t100: RL softer by 0.5 kT at 1.2 a).  The crossing-angle profile is even in psi: no chiral sector to speak of."""
import numpy as np, sys
path = sys.argv[1] if len(sys.argv) > 1 else "/home/reyhaneh/Documents/Zimm_brag_plus/helix_pair_db_t45_win.bin"
raw = open(path, "rb").read(4096); off = raw.index(b"END\n") + 4
hdr = raw[:off].decode().splitlines()
shape = [int(x) for x in [l for l in hdr if l.startswith("shape")][0].split()[1:]]
ax = {}
for l in hdr:
    if l.startswith("axis"):
        t = l.split(); ax[t[1]] = np.arange(int(t[3])) * float(t[7]) + float(t[5])
U = np.memmap(path, dtype="<f4", mode="r", offset=off, shape=tuple(shape))
r, bA, bB, psi, aA, aB = (ax[k] for k in ("r", "betaA", "betaB", "psi", "alphaA", "alphaB"))
TYPES = ["RR", "LL", "RL", "LR"]
wb = np.sin(np.radians(bA)); wb[0] = wb[-1] = 0.0            # isotropic rod orientations: sin beta weights
W = (wb[:, None] * wb[None, :]); W /= W.sum()
i90 = int(np.argmin(abs(bA - 90))); i0 = int(np.argmin(abs(psi - 0))); i180 = int(np.argmin(abs(psi - 180)))
F = {}; Fpar = {}; Upar = {}; Fanti = {}; Uanti = {}; Fpsi = {}
for t, name in enumerate(TYPES):
    u = np.asarray(U[0, t], dtype=np.float32)                 # (r, bA, bB, psi, aA, aB)
    b = np.exp(-u.astype(np.float64))
    ba = b.mean(axis=(3, 4, 5))                               # average over psi, alphaA, alphaB (uniform)
    F[name] = -np.log(np.einsum("rab,ab->r", ba, W))          # isotropic angular average, registries free
    side = u[:, i90, i90]                                     # (r, psi, aA, aB) side by side
    Fpar[name] = -np.log(np.exp(-side[:, i0].astype(np.float64)).mean(axis=(1, 2)))     # parallel, registries free
    Upar[name] = side[:, i0].min(axis=(1, 2))                                              # parallel, best registries
    Fanti[name] = -np.log(np.exp(-side[:, i180].astype(np.float64)).mean(axis=(1, 2)))
    Uanti[name] = side[:, i180].min(axis=(1, 2))
    Fpsi[name] = -np.log(np.exp(-side.astype(np.float64)).mean(axis=(2, 3)))             # (r, psi), registries free
    del u, b
print("table:", path.split("/")[-1]); print("\n".join(l for l in hdr if l.startswith(("level", "validity"))))
print("\nmirror check (should vanish): max|F_RR-F_LL| = %.3f, max|F_RL-F_LR| = %.3f (isotropic averages)" % (abs(F["RR"] - F["LL"]).max(), abs(F["RL"] - F["LR"]).max()))
print("\n(1) isotropic orientation + free registries:  F(r) = -ln<exp(-U)>  [kT]")
print(f"{'r[a]':>5} {'F_RR':>7} {'F_RL':>7} {'dF=RL-RR':>9}   | (2) side by side (beta 90/90), PARALLEL psi=0:  {'F_RR':>6} {'F_RL':>6} {'dF':>6} | best registries {'U_RR':>6} {'U_RL':>6} {'dU':>6} | ANTIparallel free-reg {'F_RR':>6} {'F_RL':>6} {'dF':>6}")
for i, x in enumerate(r):
    print(f"{x:>5.1f} {F['RR'][i]:>7.2f} {F['RL'][i]:>7.2f} {F['RL'][i]-F['RR'][i]:>9.2f}   |{'':>52}{Fpar['RR'][i]:>6.2f} {Fpar['RL'][i]:>6.2f} {Fpar['RL'][i]-Fpar['RR'][i]:>6.2f} |{'':>16}{Upar['RR'][i]:>6.2f} {Upar['RL'][i]:>6.2f} {Upar['RL'][i]-Upar['RR'][i]:>6.2f} |{'':>22}{Fanti['RR'][i]:>6.2f} {Fanti['RL'][i]:>6.2f} {Fanti['RL'][i]-Fanti['RR'][i]:>6.2f}")
# association constants: K = int 4 pi r^2 <exp(-U)>_iso dr over the well, relative to ideal
for lo in (1.2, 1.9):
    m = r >= lo - 1e-6
    K0 = np.trapezoid(4 * np.pi * r[m] ** 2, r[m])
    K = {n: np.trapezoid(4 * np.pi * r[m] ** 2 * np.exp(-F[n][m]), r[m]) for n in ("RR", "RL")}
    print(f"\nassociation over r in [{lo}, 3.0] a (isotropic, free registries): K_RR/K_ideal = {K['RR']/K0:.3f}, K_RL/K_ideal = {K['RL']/K0:.3f}, "
          f"-ln(K_RR/K_RL) = {-np.log(K['RR']/K['RL']):+.2f} kT  (B2 excess: RR {K['RR']-K0:+.2f}, RL {K['RL']-K0:+.2f} a^3)")
print("\n(3) crossing-angle profile side by side, free registries, F(psi) - F(psi=0) [kT]  (psi > 0: (u1 x u2).r > 0)")
for x in (2.0, 2.2, 2.5):
    i = int(np.argmin(abs(r - x)))
    print(f" r = {x}:  psi  " + " ".join(f"{p:>6.0f}" for p in psi))
    for n in ("RR", "RL"):
        print(f"          {n}   " + " ".join(f"{v:>6.2f}" for v in Fpsi[n][i] - Fpsi[n][i, i0]) + f"   (F(0) = {Fpsi[n][i, i0]:.2f})")
