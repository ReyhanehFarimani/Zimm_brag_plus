#!/usr/bin/env python3
"""BONDED (intra-arm) twist of the STAR, R and L separately (user 2026-09-24: "non bonded and bonded separately"),
the counterpart of the inter-arm (non-bonded) cholesteric analysis in star_cholesteric.py.  From the dumped frames
(positions, tangents and registries; last half of every run's dumps, J = 9 also the pivot control from sweep 2000 on):
  (a) backbone torsion tau_i = dihedral of the four consecutive residues (i-1, i, i+1, i+2) of one arm, for
      stretches where all four are helical of ONE hand: <tau> for R and for L separately (mirror images if the arm
      supercoils with a chiral sense; the bend potential has theta0 = 180 deg and no dihedral term, so any chirality
      here is emergent) and the chirality-weighted <chi tau>; also the coil stretches (achiral reference);
  (b) registry rotation d_alpha_i between consecutive same-handed residues (the angle from m_i to m_{i+1} about the
      bond, right-handed): the twist term imposes +25.3 deg for R.R and -25.3 for L.L at theta0 = 45 (kappa = 6.6
      kT/rad^2, sigma = 22 deg), so this is a check of the term, not a discovery;
  (c) the distributions of tau (helical R, helical L, coil) at the J with the most helical stretches.
  cd runs/t45/star50 && ~/venv/bin/python ../../star_bondedtwist.py     -> bondedtwist_star.png + table
"""
import glob
import os
import re
import sys
import warnings

import numpy as np
import matplotlib
warnings.filterwarnings("ignore", category=RuntimeWarning)
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def read_kv(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
    return d


def frames(path):
    out = []
    with open(path) as fh:
        while True:
            l = fh.readline()
            if not l:
                break
            n = int(l); hdr = fh.readline()
            sw = int(hdr.split("sweep=")[1].split()[0]) if "sweep=" in hdr else -1
            sp = []; pos = np.empty((n, 3)); reg = np.zeros((n, 3))
            for i in range(n):
                t = fh.readline().split(); sp.append(t[0]); pos[i] = (float(t[1]), float(t[2]), float(t[3]))
                if len(t) >= 14:
                    reg[i] = (float(t[11]), float(t[12]), float(t[13]))
            out.append((sw, np.array(sp), pos, reg))
    return out


def dihedral(p0, p1, p2, p3):
    """signed dihedral of the four points (deg), right-handed about p1 -> p2"""
    b0 = p1 - p0; b1 = p2 - p1; b2 = p3 - p2
    n1 = np.cross(b0, b1); n2 = np.cross(b1, b2)
    b1n = b1 / np.maximum(np.linalg.norm(b1, axis=1), 1e-12)[:, None]
    m1 = np.cross(n1, b1n)
    return np.degrees(np.arctan2((m1 * n2).sum(1), (n1 * n2).sum(1)))


def signed_angle_about(a, b, axis):
    """angle from a to b about axis (deg), right-handed, both projected perpendicular to axis"""
    ax = axis / np.maximum(np.linalg.norm(axis, axis=1), 1e-12)[:, None]
    a = a - (a * ax).sum(1)[:, None] * ax; b = b - (b * ax).sum(1)[:, None] * ax
    return np.degrees(np.arctan2((np.cross(a, b) * ax).sum(1), (a * b).sum(1)))


INPUTS = sorted(glob.glob("inputs/J*_s*.dat"))
if not INPUTS:
    sys.exit("run from runs/<geom>/star<n>")
kv = read_kv(INPUTS[0]); NARM, F = int(kv["N"]), int(kv["n_arms"])
GEOM = os.path.basename(os.path.dirname(os.getcwd()))

groups = {}
def add(J, fr_list, tag):
    g = groups.setdefault(J, dict(tau_R=[], tau_L=[], tau_C=[], da_R=[], da_L=[], nfr=0, runs=[]))
    for sw, sp, pos, reg in fr_list:
        m = sp != "X"; s = sp[m]; p = pos[m].reshape(F, NARM, 3); rg = reg[m].reshape(F, NARM, 3)
        st = np.where(s == "R", 1, np.where(s == "L", -1, 0)).reshape(F, NARM)
        # torsion of (i-1, i, i+1, i+2) for i = 1 .. N-3
        tau = dihedral(p[:, :-3].reshape(-1, 3), p[:, 1:-2].reshape(-1, 3), p[:, 2:-1].reshape(-1, 3), p[:, 3:].reshape(-1, 3)).reshape(F, NARM - 3)
        q = st[:, :-3]; allR = (q == 1) & (st[:, 1:-2] == 1) & (st[:, 2:-1] == 1) & (st[:, 3:] == 1)
        allL = (q == -1) & (st[:, 1:-2] == -1) & (st[:, 2:-1] == -1) & (st[:, 3:] == -1)
        allC = (q == 0) & (st[:, 1:-2] == 0) & (st[:, 2:-1] == 0) & (st[:, 3:] == 0)
        g["tau_R"].append(tau[allR]); g["tau_L"].append(tau[allL]); g["tau_C"].append(tau[allC])
        # registry rotation between consecutive same-handed residues about the bond
        bond = (p[:, 1:] - p[:, :-1]).reshape(-1, 3)
        da = signed_angle_about(rg[:, :-1].reshape(-1, 3), rg[:, 1:].reshape(-1, 3), bond).reshape(F, NARM - 1)
        rr = (st[:, :-1] == 1) & (st[:, 1:] == 1); ll = (st[:, :-1] == -1) & (st[:, 1:] == -1)
        g["da_R"].append(da[rr]); g["da_L"].append(da[ll]); g["nfr"] += 1
    g["runs"].append(f"{tag}×{len(fr_list)}")

for f in sorted(glob.glob("out/J*_s*_conf.xyz")):
    fr = frames(f)
    if not fr:
        continue
    J = float(re.match(r"J([\d.]+)_s", os.path.basename(f)).group(1))
    add(J, fr[len(fr) // 2:] if len(fr) > 1 else fr, os.path.basename(f)[:-9])
pc = "ctrl_pivot/out/J9_piv50_conf.xyz"
if os.path.isfile(pc):
    fr = [x for x in frames(pc) if x[0] >= 2000]
    if fr:
        add(9.0, fr, "pivot-control")
if not groups:
    sys.exit("no frames yet")


def mse(v, nmin=20): v = np.asarray(v, float); return (v.mean(), v.std(ddof=1) / np.sqrt(len(v))) if len(v) >= nmin else (np.nan, np.nan)


res = {}
for J, g in sorted(groups.items()):
    tR = np.concatenate(g["tau_R"]) if g["tau_R"] else np.array([]); tL = np.concatenate(g["tau_L"]) if g["tau_L"] else np.array([]); tC = np.concatenate(g["tau_C"]) if g["tau_C"] else np.array([])
    dR = np.concatenate(g["da_R"]) if g["da_R"] else np.array([]); dL = np.concatenate(g["da_L"]) if g["da_L"] else np.array([])
    chi_tau = np.concatenate([tR, -tL]) if len(tR) + len(tL) else np.array([])
    res[J] = dict(tR=mse(tR), tL=mse(tL), tC=mse(tC), chi=mse(chi_tau), dR=mse(dR), dL=mse(dL), sR=(dR.std() if len(dR) else np.nan), sL=(dL.std() if len(dL) else np.nan),
                  nR=len(tR), nL=len(tL), nC=len(tC), hR=tR, hL=tL, hC=tC, nfr=g["nfr"])

Js = sorted(res)
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 7.5, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
INK, MUTED, RED, BLUE = "#0b0b0b", "#898781", "#c0392f", "#1c5cab"
fig, (a, b, c) = plt.subplots(1, 3, figsize=(15.0, 4.4)); X = np.arange(len(Js))
a.axhline(0, color=MUTED, lw=0.8)
a.errorbar(X - 0.15, [res[J]["tR"][0] for J in Js], yerr=[res[J]["tR"][1] for J in Js], fmt="o", ms=5, color=RED, ecolor=RED, elinewidth=0.9, capsize=2, label="R helical stretches")
a.errorbar(X, [res[J]["tL"][0] for J in Js], yerr=[res[J]["tL"][1] for J in Js], fmt="s", ms=5, color=BLUE, ecolor=BLUE, elinewidth=0.9, capsize=2, label="L helical stretches")
a.errorbar(X + 0.15, [res[J]["tC"][0] for J in Js], yerr=[res[J]["tC"][1] for J in Js], fmt="^", ms=5, mfc="white", color=MUTED, ecolor=MUTED, elinewidth=0.9, capsize=2, label="coil stretches")
a.plot(X, [res[J]["chi"][0] for J in Js], "k_", ms=14, mew=1.5, label="χ-weighted (R − L)/2")
a.set_xticks(X); a.set_xticklabels([f"{J:g}" for J in Js]); a.set_xlabel("coupling J [k_BT]"); a.set_ylabel(r"mean backbone torsion $\langle \tau \rangle$ of 4 consecutive residues [deg]")
a.text(0.03, 0.05, "(a)  bonded: supercoiling of the arm; mirror images = chiral", transform=a.transAxes, fontsize=8, color=INK); a.legend(loc="upper right", handlelength=1.5)
b.axhline(25.3, color=RED, lw=0.8, ls=":"); b.axhline(-25.3, color=BLUE, lw=0.8, ls=":"); b.axhline(0, color=MUTED, lw=0.8)
b.errorbar(X - 0.1, [res[J]["dR"][0] for J in Js], yerr=[res[J]["sR"] for J in Js], fmt="o", ms=5, color=RED, ecolor=RED, elinewidth=0.9, capsize=2, label="R·R junctions (bar = s.d.)")
b.errorbar(X + 0.1, [res[J]["dL"][0] for J in Js], yerr=[res[J]["sL"] for J in Js], fmt="s", ms=5, color=BLUE, ecolor=BLUE, elinewidth=0.9, capsize=2, label="L·L junctions (bar = s.d.)")
b.set_xticks(X); b.set_xticklabels([f"{J:g}" for J in Js]); b.set_xlabel("coupling J [k_BT]"); b.set_ylabel(r"registry rotation per residue $\Delta\alpha$ [deg]"); b.set_ylim(-60, 60)
b.text(0.03, 0.05, "(b)  bonded: the twist term's ±25.3° (dotted), σ = 22° expected", transform=b.transAxes, fontsize=8, color=INK); b.legend(loc="upper right", handlelength=1.5)
Jh = max(Js, key=lambda J: res[J]["nR"] + res[J]["nL"]); bins = np.arange(-180, 180.1, 10)
for key, lab, colr, ls in (("hR", "R helical", RED, "-"), ("hL", "L helical", BLUE, "-"), ("hC", "coil", MUTED, "--")):
    v = res[Jh][key]
    if len(v) > 20:
        h, _ = np.histogram(v, bins, density=True); c.step(bins, np.r_[h, h[-1]], where="post", color=colr, lw=1.4, ls=ls, label=f"{lab}  (n = {len(v)}, mean {v.mean():+.1f}°)")
c.axvline(0, color=MUTED, lw=0.8); c.set_xlim(-180, 180); c.set_xlabel(r"backbone torsion $\tau$ [deg]"); c.set_ylabel("probability density [1/deg]")
c.text(0.03, 0.05, f"(c)  J = {Jh:g}", transform=c.transAxes, fontsize=8, color=INK); c.legend(loc="upper right", handlelength=1.5)
fig.tight_layout(); fig.savefig("bondedtwist_star.png", dpi=170); fig.savefig("bondedtwist_star.pdf"); plt.close(fig)

print(f"{GEOM} star f = {F}: BONDED (intra-arm) twist, R and L separately; frames = last half of each run (+ pivot control at J = 9)")
print(f"{'J':>4} {'frames':>6} {'<tau> R':>8} {'+-':>5} {'<tau> L':>8} {'+-':>5} {'<tau> coil':>10} {'+-':>5} {'chi-weighted':>12} {'+-':>5} {'n_R/n_L/n_C':>16} {'<da> R.R':>8} {'sd':>5} {'<da> L.L':>8} {'sd':>5}")
for J in Js:
    r = res[J]
    print(f"{J:>4g} {r['nfr']:>6} {r['tR'][0]:>8.2f} {r['tR'][1]:>5.2f} {r['tL'][0]:>8.2f} {r['tL'][1]:>5.2f} {r['tC'][0]:>10.2f} {r['tC'][1]:>5.2f} {r['chi'][0]:>12.2f} {r['chi'][1]:>5.2f} "
          f"{r['nR']:>5}/{r['nL']:<5}/{r['nC']:<5} {r['dR'][0]:>8.1f} {r['sR']:>5.1f} {r['dL'][0]:>8.1f} {r['sL']:>5.1f}")
