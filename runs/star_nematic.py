#!/usr/bin/env python3
"""Orientational order of the STAR's helical rods (user 2026-09-24: "nematic order as well"), from the dumped
frames of runs/<geom>/star<n> (tangent columns of the extended xyz; last half of every run's dumps, J = 9 also the
pivot control from sweep 2000 on):
  (a) radial order  S_rad(r) = <P2(u . r_hat)> per shell, helical residues (solid) and all residues (dashed):
      1 = every rod points along the radius (rigid-rod star), 0 = isotropic;
  (b) local nematic order between helices of DIFFERENT arms closer than R_LOC: S_loc(r) = <P2(u_i . u_j)> per
      shell of the pair midpoint: 1 = parallel bundles, 0 = uncorrelated, -1/2 = perpendicular;
  (c) chiral crossing of CONTACTING helix pairs of different arms (r < 3 a, the table cutoff): the crossing angle
      psi from the table's convention (sin psi = (u_i x u_j) . r_hat / (sin bA sin bB)); <chi sin psi> per J for
      same-handed pairs with chi = +1 (R.R) / -1 (L.L) -- the mirror symmetry makes R.R and L.L twist oppositely,
      so they must be combined with opposite signs, not averaged -- and <sin psi> for opposite-handed pairs (zero by
      symmetry, a noise reference); 0 = achiral arrangement, plus <P2(u_i . u_j)> of the same contacts.
Also printed: the GLOBAL nematic order S (largest eigenvalue of Q = 3/2 <u u> - 1/2 over all helical rods), which
must be ~0 for a radially symmetric star.
  cd runs/t45/star50 && ~/venv/bin/python ../../star_nematic.py     -> nematic_star.png + table
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

RC, R_LOC = 3.0, 5.0        # contact cutoff (table r_max) and local-order neighbourhood [a]


def read_kv(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
    return d


def frames(path):
    """(sweep, species, pos, tangent) per frame; columns: species pos(3) orientation(4) tangent(3) ..."""
    out = []
    with open(path) as fh:
        while True:
            l = fh.readline()
            if not l:
                break
            n = int(l); hdr = fh.readline()
            sw = int(hdr.split("sweep=")[1].split()[0]) if "sweep=" in hdr else -1
            sp = []; pos = np.empty((n, 3)); tan = np.zeros((n, 3))
            for i in range(n):
                t = fh.readline().split(); sp.append(t[0]); pos[i] = (float(t[1]), float(t[2]), float(t[3]))
                if len(t) >= 11:
                    tan[i] = (float(t[8]), float(t[9]), float(t[10]))
            out.append((sw, np.array(sp), pos, tan))
    return out


def P2(c): return 1.5 * c * c - 0.5


INPUTS = sorted(glob.glob("inputs/J*_s*.dat"))
if not INPUTS:
    sys.exit("run from runs/<geom>/star<n>")
kv = read_kv(INPUTS[0])
NARM, F, R0 = int(kv["N"]), int(kv["n_arms"]), float(kv["core_radius"]) + 0.5
BC = float(kv["bond_len_CC"])
GEOM = os.path.basename(os.path.dirname(os.getcwd()))
edges = np.geomspace(R0 - 0.5, R0 + NARM * BC + 5, 24); rc = np.sqrt(edges[:-1] * edges[1:]); NB = len(rc)

groups = {}
def add(J, fr_list, tag):
    g = groups.setdefault(J, dict(rad_h=np.zeros(NB), n_h=np.zeros(NB), rad_a=np.zeros(NB), n_a=np.zeros(NB), loc=np.zeros(NB), n_loc=np.zeros(NB),
                                  sin_same=[], sin_opp=[], p2_same=[], p2_opp=[], S=[], nfr=0, runs=[]))
    for sw, sp, pos, tan in fr_list:
        m = sp != "X"; s = sp[m]; p = pos[m]; u = tan[m]
        nu = np.linalg.norm(u, axis=1); ok = nu > 1e-9; u[ok] /= nu[ok][:, None]
        r = np.linalg.norm(p, axis=1); rhat = p / np.maximum(r, 1e-12)[:, None]
        hel = (s != "C") & ok; arm = np.arange(len(s)) // NARM; sign = np.where(s == "R", 1, np.where(s == "L", -1, 0))
        b = np.clip(np.searchsorted(edges, r, side="right") - 1, 0, NB); ins = b < NB
        c = (u * rhat).sum(1)
        g["rad_a"] += np.bincount(b[ins & ok], weights=P2(c[ins & ok]), minlength=NB); g["n_a"] += np.bincount(b[ins & ok], minlength=NB)
        g["rad_h"] += np.bincount(b[ins & hel], weights=P2(c[ins & hel]), minlength=NB); g["n_h"] += np.bincount(b[ins & hel], minlength=NB)
        # global nematic order of the helical rods
        if hel.sum() > 1:
            Q = 1.5 * (u[hel].T @ u[hel]) / hel.sum() - 0.5 * np.eye(3); g["S"].append(np.linalg.eigvalsh(Q).max())
        # helix pairs of different arms within R_LOC (brute force in chunks)
        H = np.where(hel)[0]; PH = p[H]; UH = u[H]; AH = arm[H]; SH = sign[H]
        for i0 in range(0, len(H), 400):
            i1 = min(len(H), i0 + 400)
            d = PH[i0:i1, None, :] - PH[None, :, :]; d2 = (d * d).sum(-1)
            pair = (d2 < R_LOC * R_LOC) & (AH[i0:i1, None] != AH[None, :])
            pair &= np.arange(i0, i1)[:, None] < np.arange(len(H))[None, :]          # each pair once (i < j)
            ii, jj = np.where(pair)
            if len(ii) == 0:
                continue
            ii += i0
            cij = (UH[ii] * UH[jj]).sum(1); mid = 0.5 * (PH[ii] + PH[jj]); rm = np.linalg.norm(mid, axis=1)
            bm = np.clip(np.searchsorted(edges, rm, side="right") - 1, 0, NB); okm = bm < NB
            g["loc"] += np.bincount(bm[okm], weights=P2(cij[okm]), minlength=NB); g["n_loc"] += np.bincount(bm[okm], minlength=NB)
            # contacts: chiral crossing angle
            dist = np.sqrt(d2[ii - i0, jj]); ct = dist < RC
            if ct.any():
                rh = (PH[jj[ct]] - PH[ii[ct]]) / dist[ct][:, None]
                e1 = (rh * UH[ii[ct]]).sum(1); e2 = (rh * UH[jj[ct]]).sum(1)
                cr = np.cross(UH[ii[ct]], UH[jj[ct]]); num = (cr * rh).sum(1); den = cij[ct] - e1 * e2
                psi = np.arctan2(num, den); same = SH[ii[ct]] == SH[jj[ct]]; chi = SH[ii[ct]]
                g["sin_same"].extend(chi[same] * np.sin(psi[same])); g["sin_opp"].extend(np.sin(psi[~same]))   # same-handed: chi-weighted (R.R +, L.L -)
                g["p2_same"].extend(P2(cij[ct][same])); g["p2_opp"].extend(P2(cij[ct][~same]))
        g["nfr"] += 1
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

Js = sorted(groups); cmap = plt.get_cmap("viridis"); col = {J: cmap(0.1 + 0.8 * i / max(1, len(Js) - 1)) for i, J in enumerate(Js)}
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 7.5, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
INK, MUTED = "#0b0b0b", "#898781"
fig, (a, b, c) = plt.subplots(1, 3, figsize=(15.0, 4.4))
for J in Js:
    g = groups[J]
    mh = g["n_h"] > 20; ma = g["n_a"] > 20
    a.plot(rc[ma], g["rad_a"][ma] / g["n_a"][ma], "--", color=col[J], lw=1.0)
    a.plot(rc[mh], g["rad_h"][mh] / g["n_h"][mh], "-", color=col[J], lw=1.6, label=f"J = {J:g}")
    ml = g["n_loc"] > 20
    b.plot(rc[ml], g["loc"][ml] / g["n_loc"][ml], "-", color=col[J], lw=1.6, label=f"J = {J:g}")
a.plot([], [], "-", color=INK, lw=1.6, label="helical residues"); a.plot([], [], "--", color=INK, lw=1.0, label="all residues")
a.axhline(0, color=MUTED, lw=0.8); a.axhline(1, color=MUTED, lw=0.8); a.axvline(R0, color=MUTED, lw=0.7)
a.set_xscale("log"); a.set_ylim(-0.55, 1.05); a.set_xlabel("distance from the core centre r [a]"); a.set_ylabel(r"radial order $\langle P_2(\hat u \cdot \hat r) \rangle$")
a.text(0.03, 0.05, "(a)  1 = radial rods, 0 = isotropic", transform=a.transAxes, fontsize=8, color=INK); a.legend(loc="lower right", ncol=2, handlelength=1.5)
b.axhline(0, color=MUTED, lw=0.8); b.axvline(R0, color=MUTED, lw=0.7)
b.set_xscale("log"); b.set_ylim(-0.55, 1.05); b.set_xlabel("distance of the pair midpoint from the core r [a]"); b.set_ylabel(r"local nematic order $\langle P_2(\hat u_i \cdot \hat u_j) \rangle$, pairs < %g a" % R_LOC)
b.text(0.03, 0.05, "(b)  helix pairs of different arms", transform=b.transAxes, fontsize=8, color=INK); b.legend(loc="upper right", ncol=2, handlelength=1.5)
X = np.arange(len(Js))
def mse(v): v = np.asarray(v); return (v.mean(), v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else (np.nan, np.nan)
ss = [mse(groups[J]["sin_same"]) for J in Js]; so = [mse(groups[J]["sin_opp"]) for J in Js]
c.axhline(0, color=MUTED, lw=0.8)
c.errorbar(X - 0.12, [x[0] for x in ss], yerr=[x[1] for x in ss], fmt="o", ms=5, color="#c0392f", ecolor="#c0392f", elinewidth=0.9, capsize=2, label="same-handed contacts, χ-weighted (R·R +, L·L −)")
c.errorbar(X + 0.12, [x[0] for x in so], yerr=[x[1] for x in so], fmt="s", ms=5, mfc="white", color="#1c5cab", ecolor="#1c5cab", elinewidth=0.9, capsize=2, label="opposite-handed contacts (R·L)")
c.set_xticks(X); c.set_xticklabels([f"{J:g}" for J in Js]); c.set_xlabel("coupling J [k_BT]"); c.set_ylabel(r"$\langle \chi \sin\psi \rangle$ of contacting helix pairs (< 3 a)")
c.set_ylim(-0.5, 0.5); c.text(0.03, 0.05, "(c)  0 = no preferred twist sense between arms", transform=c.transAxes, fontsize=8, color=INK); c.legend(loc="lower right", handlelength=1.5)
for i, J in enumerate(Js):
    c.text(X[i], 0.44, f"{len(groups[J]['sin_same'])}/{len(groups[J]['sin_opp'])}", ha="center", fontsize=7, color=MUTED)
fig.tight_layout(); fig.savefig("nematic_star.png", dpi=170); fig.savefig("nematic_star.pdf"); plt.close(fig)

print(f"{GEOM} star f = {F}: orientational order of the helical rods; frames = last half of each run (+ pivot control at J = 9)")
print(f"{'J':>4} {'frames':>6} {'S_global':>8} {'S_rad r<8':>9} {'S_rad 8-30':>10} {'S_rad >30':>9} {'S_loc r<8':>9} {'S_loc 8-30':>10} {'n_pairs':>7} {'<P2> contacts same/opp':>22} {'<chi sin psi> same':>18} {'+-':>5} {'opp':>7} {'+-':>5} {'n_same/n_opp':>13}")
for J in Js:
    g = groups[J]
    def band(num, den, r1, r2):
        m = (rc >= r1) & (rc < r2); return num[m].sum() / den[m].sum() if den[m].sum() > 0 else np.nan
    p2s, p2o = mse(g["p2_same"]), mse(g["p2_opp"]); s_s, s_o = mse(g["sin_same"]), mse(g["sin_opp"])
    print(f"{J:>4g} {g['nfr']:>6} {np.mean(g['S']):>8.3f} {band(g['rad_h'], g['n_h'], 0, 8):>9.3f} {band(g['rad_h'], g['n_h'], 8, 30):>10.3f} {band(g['rad_h'], g['n_h'], 30, 200):>9.3f} "
          f"{band(g['loc'], g['n_loc'], 0, 8):>9.3f} {band(g['loc'], g['n_loc'], 8, 30):>10.3f} {int(g['n_loc'].sum()):>7} {p2s[0]:>10.3f} / {p2o[0]:<9.3f} {s_s[0]:>14.3f} {s_s[1]:>5.3f} {s_o[0]:>7.3f} {s_o[1]:>5.3f} {len(g['sin_same']):>6}/{len(g['sin_opp'])}")
