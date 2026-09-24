#!/usr/bin/env python3
"""Pair correlation of the ARM handedness (user 2026-09-24: "pair correlatio"): does an arm share its hand with
its neighbours?  From the dumped frames of runs/<geom>/star<n> (last half of every run's dumps; J = 9 also the pivot
control from sweep 2000 on).  Arm k gets s_k = sign(n_R - n_L) of its helical residues (arms with n_R = n_L or no
helix are left out) and m_k = (n_R - n_L)/(n_R + n_L).
  (a) C(gamma) = <s_k s_l> for arm pairs binned by the angular separation gamma of their GRAFT sites (Fibonacci
      points, nearest neighbours ~30 deg apart for 50 arms), with the 1/sqrt(n_pairs) error of independent arms;
      C = 0 is the independent-arm null, C > 0 = like hands together, C < 0 = alternating hands.
  (b) P(same hand) for arm pairs that are IN CONTACT in that frame (any helix-helix residue pair closer than the
      pair cutoff 3 a) against pairs that are not, per J.  This is the direct test of chirality transfer through the
      helix-helix table.
Frames are treated as independent samples (they are 2000 sweeps apart).
  cd runs/t45/star50 && ~/venv/bin/python ../../star_handpair.py     -> handpair_star.png + table
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

RC = 3.0       # helix-helix pair cutoff [a] (table r_max)


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
            sp = []; pos = np.empty((n, 3))
            for i in range(n):
                t = fh.readline().split(); sp.append(t[0]); pos[i] = (float(t[1]), float(t[2]), float(t[3]))
            out.append((sw, np.array(sp), pos))
    return out


INPUTS = sorted(glob.glob("inputs/J*_s*.dat"))
if not INPUTS:
    sys.exit("run from runs/<geom>/star<n>")
kv = read_kv(INPUTS[0])
NARM, F, R0 = int(kv["N"]), int(kv["n_arms"]), float(kv["core_radius"]) + 0.5
GEOM = os.path.basename(os.path.dirname(os.getcwd()))

# graft geometry (Chain::init): Fibonacci points on the sphere of radius R0
ga = np.pi * (3.0 - np.sqrt(5.0)); k = np.arange(F)
z = 1.0 - 2.0 * (k + 0.5) / F; rr = np.sqrt(np.maximum(0.0, 1.0 - z * z)); ph = ga * k
U = np.stack([rr * np.cos(ph), rr * np.sin(ph), z], 1)
iu = np.triu_indices(F, 1)
GAMMA = np.degrees(np.arccos(np.clip((U @ U.T)[iu], -1, 1)))          # angular separation of every arm pair
GBINS = np.array([0, 40, 60, 90, 120, 150, 180.001]); gc = 0.5 * (GBINS[:-1] + GBINS[1:])

groups = {}
def add(J, fr_list, tag):
    g = groups.setdefault(J, dict(ss=[], mm=[], contact=[], runs=[], nfr=0))
    for sw, sp, pos in fr_list:
        m = sp != "X"; s = sp[m]; p = pos[m]
        sign = np.where(s == "R", 1, np.where(s == "L", -1, 0)).reshape(F, NARM); P = p.reshape(F, NARM, 3)
        nR = (sign > 0).sum(1); nL = (sign < 0).sum(1); nh = nR + nL
        sk = np.sign(nR - nL); mk = np.where(nh > 0, (nR - nL) / np.maximum(nh, 1), 0.0)
        ss = np.outer(sk, sk)[iu]; mm = np.outer(mk, mk)[iu]; valid = (np.outer(sk != 0, sk != 0))[iu]
        # helix-helix contacts between arms: any pair of helical residues closer than RC
        hel = sign != 0; ct = np.zeros(len(iu[0]), bool)
        for q, (a, b) in enumerate(zip(*iu)):
            if not (valid[q] and hel[a].any() and hel[b].any()):
                continue
            d2 = ((P[a][hel[a]][:, None, :] - P[b][hel[b]][None, :, :]) ** 2).sum(-1)
            ct[q] = (d2 < RC * RC).any()
        g["ss"].append(np.where(valid, ss, np.nan)); g["mm"].append(np.where(valid, mm, np.nan)); g["contact"].append(ct); g["nfr"] += 1
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

res = {}
for J, g in sorted(groups.items()):
    ss = np.concatenate(g["ss"]); mm = np.concatenate(g["mm"]); ct = np.concatenate(g["contact"]); gam = np.tile(GAMMA, g["nfr"])
    ok = np.isfinite(ss)
    C = np.full(len(gc), np.nan); E = np.full(len(gc), np.nan); n_b = np.zeros(len(gc), int)
    for i in range(len(gc)):
        m = ok & (gam >= GBINS[i]) & (gam < GBINS[i + 1])
        n_b[i] = m.sum()
        if n_b[i] > 0:
            C[i] = ss[m].mean(); E[i] = 1.0 / np.sqrt(n_b[i])
    same_ct = ss[ok & ct]; same_no = ss[ok & ~ct]
    res[J] = dict(C=C, E=E, n=n_b, nfr=g["nfr"],
                  p_ct=(0.5 * (1 + same_ct.mean()) if len(same_ct) else np.nan), e_ct=(0.5 / np.sqrt(len(same_ct)) if len(same_ct) else np.nan), n_ct=len(same_ct),
                  p_no=(0.5 * (1 + same_no.mean()) if len(same_no) else np.nan), e_no=(0.5 / np.sqrt(len(same_no)) if len(same_no) else np.nan), n_no=len(same_no),
                  Cm=np.nanmean(mm) / max(1e-12, np.nanmean(np.abs(mm))), Cnn=C[0], Enn=E[0],
                  arms_ok=ok.mean())

Js = sorted(res); cmap = plt.get_cmap("viridis"); col = {J: cmap(0.1 + 0.8 * i / max(1, len(Js) - 1)) for i, J in enumerate(Js)}
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 7.5, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
INK, MUTED = "#0b0b0b", "#898781"
fig, (a, b) = plt.subplots(1, 2, figsize=(11.0, 4.3))
a.axhline(0, color=MUTED, lw=0.8)
for i, J in enumerate(Js):
    off = (i - 0.5 * (len(Js) - 1)) * 2.0
    a.errorbar(gc + off, res[J]["C"], yerr=res[J]["E"], fmt="o-", ms=4, lw=1.0, color=col[J], ecolor=col[J], elinewidth=0.8, capsize=2, label=f"J = {J:g}")
a.set_xlabel("angular separation of the graft sites γ [deg]  (nearest neighbours ≈ 30°)"); a.set_ylabel(r"$C(\gamma) = \langle s_k s_l \rangle$ of the arm handedness")
a.set_xlim(0, 180); a.text(0.03, 0.05, "(a)  0 = independent arms; bars = 1/√n_pairs", transform=a.transAxes, fontsize=8, color=INK)
a.legend(loc="upper right", ncol=2, handlelength=1.5)
b.axhline(0.5, color=MUTED, lw=0.8)
X = np.arange(len(Js))
b.errorbar(X - 0.12, [res[J]["p_ct"] for J in Js], yerr=[res[J]["e_ct"] for J in Js], fmt="o", ms=5, color="#c0392f", ecolor="#c0392f", elinewidth=0.9, capsize=2, label="arm pairs in helix–helix contact (< 3 a)")
b.errorbar(X + 0.12, [res[J]["p_no"] for J in Js], yerr=[res[J]["e_no"] for J in Js], fmt="s", ms=5, mfc="white", color="#1c5cab", ecolor="#1c5cab", elinewidth=0.9, capsize=2, label="arm pairs not in contact")
b.set_xticks(X); b.set_xticklabels([f"{J:g}" for J in Js]); b.set_xlabel("coupling J [k_BT]"); b.set_ylabel("P(same hand)")
b.set_ylim(0.0, 1.0); b.text(0.03, 0.05, "(b)  1/2 = independent arms", transform=b.transAxes, fontsize=8, color=INK)
b.legend(loc="lower left", handlelength=1.5)
for i, J in enumerate(Js):
    b.text(X[i], 0.93, f"{res[J]['n_ct']}", ha="center", fontsize=7, color="#c0392f")
b.text(0.5, 0.97, "number of contacting pairs", transform=b.transAxes, ha="center", va="top", fontsize=7, color=MUTED)
fig.tight_layout(); fig.savefig("handpair_star.png", dpi=170); fig.savefig("handpair_star.pdf"); plt.close(fig)

print(f"{GEOM} star f = {F}: arm-handedness pair correlation; frames = last half of each run (+ pivot control at J = 9)")
print(f"{'J':>4} {'frames':>6} {'C(gamma<40)':>11} {'+-':>5} {'n_nn':>5} {'C(all m-weighted)':>17} {'P_same contact':>14} {'+-':>5} {'n_ct':>5} {'P_same no contact':>17} {'+-':>5} {'n_no':>6} {'arms with a hand':>16}")
for J in Js:
    r = res[J]
    print(f"{J:>4g} {r['nfr']:>6} {r['Cnn']:>11.3f} {r['Enn']:>5.3f} {r['n'][0]:>5} {r['Cm']:>17.3f} {r['p_ct']:>14.3f} {r['e_ct']:>5.3f} {r['n_ct']:>5} {r['p_no']:>17.3f} {r['e_no']:>5.3f} {r['n_no']:>6} {r['arms_ok']:>16.2f}")
