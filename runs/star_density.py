#!/usr/bin/env python3
"""Radial structure of the STAR (user 2026-09-24: "plot the density of monomers as a function of radial distance and
compare with stiff arm star polymers and soft ones"; "also handedness per r"): from the dumped frames of
runs/<geom>/star<n>,
  (a) residue number density rho(r) around the core centre, one curve per J (sequential colours, low J = coil-rich
      "soft" arms, high J = helical "stiff" arms), against
        - the RIGID-ROD star, absolute: rho = f / (4 pi b_H r^2) for R0 < r < R0 + N b_H (every arm a straight radial rod),
        - flexible-arm scalings drawn as slope guides through the softest MC profile at r = 10 a: Daoud-Cotton good
          solvent rho ~ r^-4/3 and ideal (theta) arms rho ~ r^-1 (prefactors of order one are not predicted);
  (b) mean radial position <r_i> of residue i along the arm (rod: R0 + i b_H);
  (c) local helicity theta(r), the helix fraction of the residues in each shell;
  (d) local handedness <|m(r)|>, m = (n_R - n_L)/(n_R + n_L) of the helical residues in each shell, per frame, averaged
      over frames (solid), against the INDEPENDENT-ARM null of the same frames (dashed): every arm mirrored with
      probability 1/2, 50 random mirrorings per frame.  Solid above dashed = arms of like handedness sit together.
Frames: the last half of every run's dumps (all of them when a run has only one) plus, for J = 9, the pivot control
(ctrl_pivot, its frames from sweep 2000 on).  One MC residue is a block of 7 fine monomers.
  cd runs/t45/star50 && ~/venv/bin/python ../../star_density.py     -> density_star.png + table
"""
import glob
import os
import re
import sys

import warnings

import numpy as np
import matplotlib
warnings.filterwarnings("ignore", category=RuntimeWarning)   # nanmean of shells without data
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
    """extended xyz written by the MC: (sweep, species array, positions array) per frame; the core X first"""
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
BH, BC = float(kv["bond_len_HH"]), float(kv["bond_len_CC"])
GEOM = os.path.basename(os.path.dirname(os.getcwd()))
NULL_DRAWS = 50
rng = np.random.default_rng(1)

groups = {}          # J -> dict(frames=[(r, sign, arm)], runs=[...])
def add(J, fr_list, tag):
    g = groups.setdefault(J, dict(frames=[], runs=[]))
    for sw, sp, pos in fr_list:
        core = pos[sp == "X"][0] if (sp == "X").any() else np.zeros(3)
        m = sp != "X"; p = pos[m] - core; s = sp[m]
        r = np.sqrt((p * p).sum(1)); sign = np.where(s == "R", 1, np.where(s == "L", -1, 0))
        g["frames"].append((r, sign, np.arange(len(r)) // NARM))
    g["runs"].append(f"{tag}×{len(fr_list)}")

for f in sorted(glob.glob("out/J*_s*_conf.xyz")):
    fr = frames(f)
    if not fr:
        continue
    J = float(re.match(r"J([\d.]+)_s", os.path.basename(f)).group(1))
    keep = fr[len(fr) // 2:] if len(fr) > 1 else fr
    add(J, keep, os.path.basename(f)[:-9])
pc = "ctrl_pivot/out/J9_piv50_conf.xyz"
if os.path.isfile(pc):
    fr = [x for x in frames(pc) if x[0] >= 2000]
    if fr:
        add(9.0, fr, "pivot-control")
if not groups:
    sys.exit("no frames yet")

# ---------------------------------------------------------------- binning
edges = np.geomspace(R0 - 0.5, R0 + NARM * BC + 5, 36); rc = np.sqrt(edges[:-1] * edges[1:]); vol = 4 * np.pi / 3 * (edges[1:] ** 3 - edges[:-1] ** 3)
NB = len(rc)
prof = {}
for J, g in sorted(groups.items()):
    nfr = len(g["frames"]); n = np.zeros(NB); nh = np.zeros(NB); ri = np.zeros(NARM)
    am = np.full((nfr, NB), np.nan); am0 = np.full((nfr, NB), np.nan)
    for k, (r, sign, arm) in enumerate(g["frames"]):
        b = np.clip(np.searchsorted(edges, r, side="right") - 1, 0, NB)          # NB = outside the last edge
        inside = b < NB
        n += np.bincount(b[inside], minlength=NB); nh += np.bincount(b[inside], weights=(sign[inside] != 0), minlength=NB)
        ri += np.array([r[np.arange(len(r)) % NARM == i].mean() for i in range(NARM)]) / nfr
        hel = np.bincount(b[inside], weights=np.abs(sign[inside]), minlength=NB)
        ok = hel >= 20
        s_shell = np.bincount(b[inside], weights=sign[inside], minlength=NB)
        am[k, ok] = np.abs(s_shell[ok]) / hel[ok]
        # independent-arm null: mirror every arm with probability 1/2
        flips = rng.choice([-1, 1], size=(NULL_DRAWS, F))
        acc = np.zeros(NB)
        for d in range(NULL_DRAWS):
            acc += np.abs(np.bincount(b[inside], weights=sign[inside] * flips[d][arm[inside]], minlength=NB))
        am0[k, ok] = acc[ok] / NULL_DRAWS / hel[ok]
    rho = n / (nfr * vol); th = np.where(n > 0, nh / np.maximum(n, 1), np.nan)
    prof[J] = dict(rho=rho, th=th, ri=ri, n=n, nfr=nfr, am=np.nanmean(am, 0), am0=np.nanmean(am0, 0), amn=np.sum(np.isfinite(am), 0))

Js = sorted(prof); cmap = plt.get_cmap("viridis"); col = {J: cmap(0.1 + 0.8 * i / max(1, len(Js) - 1)) for i, J in enumerate(Js)}
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 7.5, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
INK, MUTED = "#0b0b0b", "#898781"
fig, axs = plt.subplots(2, 2, figsize=(11.0, 8.6)); a, b, c, d = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]

# (a) density
rr = np.geomspace(R0, R0 + NARM * BH, 200)
a.plot(rr, F / (4 * np.pi * BH * rr ** 2), color=INK, lw=1.6, label=f"rigid radial rods, f/(4π b r²), b = {BH:.2f} a")
J0 = Js[0]; i10 = int(np.argmin(abs(rc - 10.0))); rho10 = prof[J0]["rho"][i10]
if np.isfinite(rho10) and rho10 > 0:
    a.plot(rr, rho10 * (rr / rc[i10]) ** (-4 / 3), color=MUTED, lw=1.2, ls="--", label="flexible arms, good solvent (Daoud–Cotton) ∝ r⁻⁴ᐟ³")
    a.plot(rr, rho10 * (rr / rc[i10]) ** (-1.0), color=MUTED, lw=1.2, ls=":", label="ideal arms (θ solvent) ∝ r⁻¹")
for J in Js:
    m = prof[J]["n"] > 0
    a.plot(rc[m], prof[J]["rho"][m], "-", color=col[J], lw=1.6, label=f"MC J = {J:g}")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("distance from the core centre r [a]"); a.set_ylabel(r"residue number density $\rho(r)$ [a$^{-3}$]")
a.axvline(R0, color=MUTED, lw=0.7); a.axvline(R0 + NARM * BH, color=MUTED, lw=0.7)
a.text(0.97, 0.95, f"(a)  θ₀ = {GEOM[1:]}°, f = {F} arms × N = {NARM}", transform=a.transAxes, fontsize=8, color=INK, ha="right", va="top")
a.legend(loc="lower left", handlelength=1.8)

# (b) arm extension profile
b.plot(np.arange(NARM), R0 + np.arange(NARM) * BH, color=INK, lw=1.6, label="rigid radial rod")
for J in Js:
    b.plot(np.arange(NARM), prof[J]["ri"], "-", color=col[J], lw=1.6, label=f"J = {J:g}")
b.set_xlabel("residue index i along the arm (0 = graft)"); b.set_ylabel(r"$\langle r_i \rangle$ [a]")
b.text(0.03, 0.92, "(b)  mean radial position of residue i", transform=b.transAxes, fontsize=8, color=INK)
b.legend(loc="lower right", ncol=2, handlelength=1.8)

# (c) local helicity
for J in Js:
    m = prof[J]["n"] > 20
    c.plot(rc[m], prof[J]["th"][m], "-", color=col[J], lw=1.6, label=f"J = {J:g}")
c.set_xscale("log"); c.set_xlabel("distance from the core centre r [a]"); c.set_ylabel(r"local helicity $\theta(r)$"); c.set_ylim(-0.02, 1.02)
c.axvline(R0, color=MUTED, lw=0.7)
c.text(0.03, 0.05, "(c)  helix fraction per shell", transform=c.transAxes, fontsize=8, color=INK)
c.legend(loc="center right", ncol=2, handlelength=1.8)

# (d) local handedness against the independent-arm null
for J in Js:
    m = prof[J]["amn"] > 0
    d.plot(rc[m], prof[J]["am0"][m], "--", color=col[J], lw=1.0)
    d.plot(rc[m], prof[J]["am"][m], "-", color=col[J], lw=1.6, label=f"J = {J:g}")
d.plot([], [], "-", color=INK, lw=1.6, label="MC shell |m|"); d.plot([], [], "--", color=INK, lw=1.0, label="same frames, arms mirrored at random")
d.set_xscale("log"); d.set_xlabel("distance from the core centre r [a]"); d.set_ylabel(r"shell handedness $\langle |m(r)| \rangle$")
d.set_ylim(0, max(0.35, float(np.nanmax([np.nanmax(prof[J]["am"]) for J in Js])) * 1.1))
d.axvline(R0, color=MUTED, lw=0.7)
d.text(0.03, 0.03, "(d)  |n_R − n_L| / (n_R + n_L) per shell and frame", transform=d.transAxes, fontsize=8, color=INK, va="bottom")
d.legend(loc="upper right", ncol=2, handlelength=1.8)
fig.tight_layout(); fig.savefig("density_star.png", dpi=170); fig.savefig("density_star.pdf"); plt.close(fig)

# ---------------------------------------------------------------- numbers
def rho_at(J, r): i = int(np.argmin(abs(rc - r))); return prof[J]["rho"][i]
def slope(J, r1, r2):
    m = (rc >= r1) & (rc <= r2) & (prof[J]["n"] > 0)
    return np.polyfit(np.log(rc[m]), np.log(prof[J]["rho"][m]), 1)[0] if m.sum() >= 3 else float("nan")
def band(J, key, r1, r2):
    m = (rc >= r1) & (rc < r2) & (prof[J]["amn"] > 0); v = prof[J][key][m]
    return float(np.nanmean(v)) if m.any() else float("nan")
print(f"{GEOM} star f = {F}, N = {NARM}, core wall R0 = {R0:g} a; rod star: rho = f/(4 pi b r^2) with b_H = {BH:.3f} a; frames from the last half of each run")
print(f"{'J':>4} {'frames':>6} {'rho(5a)':>8} {'rho(10a)':>8} {'rho(20a)':>8} {'rho(40a)':>8} {'slope 8-40a':>11} {'<r_tip>':>8} {'tip/L_rod':>9} {'th r<8a':>8} {'th r>8a':>8} {'|m| r<10a':>9} {'null':>6} {'|m| 10-40a':>10} {'null':>6}   runs")
for J in Js:
    g = prof[J]; L = R0 + (NARM - 1) * BH
    mi = (rc < 8) & (g["n"] > 0); mo = (rc >= 8) & (g["n"] > 0)
    inner = np.nansum(g["th"][mi] * g["n"][mi]) / g["n"][mi].sum(); outer = np.nansum(g["th"][mo] * g["n"][mo]) / g["n"][mo].sum()
    print(f"{J:>4g} {g['nfr']:>6} {rho_at(J, 5):>8.4f} {rho_at(J, 10):>8.4f} {rho_at(J, 20):>8.4f} {rho_at(J, 40):>8.4f} {slope(J, 8, 40):>11.2f} {g['ri'][-1]:>8.1f} {g['ri'][-1] / L:>9.2f} {inner:>8.3f} {outer:>8.3f} "
          f"{band(J, 'am', 0, 10):>9.3f} {band(J, 'am0', 0, 10):>6.3f} {band(J, 'am', 10, 40):>10.3f} {band(J, 'am0', 10, 40):>6.3f}   {', '.join(groups[J]['runs'])}")
print(f"rod star at 5 / 10 / 20 / 40 a: " + " / ".join(f"{F / (4 * np.pi * BH * r * r):.4f}" for r in (5, 10, 20, 40)) + "  (slope -2); Daoud-Cotton good solvent slope -4/3, theta -1")
