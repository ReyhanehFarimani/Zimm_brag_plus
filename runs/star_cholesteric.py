#!/usr/bin/env python3
"""Cholesteric (twist) order of the STAR's helical rods (user 2026-09-24: "compute the cholesteric order as well").
A cholesteric arrangement means neighbouring rods are rotated about their connecting vector with a PREFERRED sense,
by an angle that grows with the separation: psi(r) = q r, pitch P = 2 pi / |q|.  From the dumped frames of
runs/<geom>/star<n> (tangents; last half of every run's dumps, J = 9 also the pivot control from sweep 2000 on),
for helix pairs on DIFFERENT arms closer than R_MAX:
  psi     = dihedral of the two axes about r_hat_ij (table convention: sin psi ~ (u_i x u_j).r_hat, cos psi ~
            u_i.u_j - (r.u_i)(r.u_j)); the arms are polar (graft -> tip) so psi itself is meaningful, and the
            APOLAR twist psi_a = psi wrapped into (-90, 90] deg is the nematic-director version;
  R and L SEPARATELY (user 2026-09-24: "did u consider r and l separately?"): by mirror symmetry an R.R pair twists
  with the opposite sense of an L.L pair and the star has as many L arms as R arms, so a plain average over
  same-handed pairs cancels the signal; R.L pairs carry none (their ensemble is mapped onto itself by the mirror).
  Every twist quantity below is therefore weighted by the pair chirality chi = +1 (R.R), -1 (L.L), and R.R / L.L
  are also reported on their own (they must come out as mirror images):
  (a) twist correlation T_chi(r) = <chi sin psi cos psi>, the pseudoscalar chiral correlation of apolar nematogens
      (0 = achiral arrangement), with R.R and L.L separately in the table;
  (b) mean chirality-weighted apolar twist <chi psi_a>(r) of SIDE-BY-SIDE pairs (|r_hat.u| < 0.5 for both rods,
      where the dihedral is a clean twist angle), straight-line fit through the origin over r <= R_FIT: slope
      q = twist rate of the R.R pairs (L.L: -q), pitch P = 360 / |q|; q within 2 sigma of zero -> lower bound on P;
  (c) distributions of psi_a for close side-by-side pairs (r < 4 a): R.R, L.L and R.L at the J with the most pairs;
      a cholesteric tendency skews R.R and L.L to opposite sides.
  cd runs/t45/star50 && ~/venv/bin/python ../../star_cholesteric.py     -> cholesteric_star.png + table
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

R_MAX, R_FIT, SIDE = 15.0, 8.0, 0.5
R_CLOSE = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0      # 'close pair' cutoff for the table and panel (c); user 2026-09-24: also 'under 2a'
SUFFIX = "" if R_CLOSE == 4.0 else f"_r{R_CLOSE:g}"
REDGES = np.arange(1.0, R_MAX + 0.001, 1.0); rc = 0.5 * (REDGES[:-1] + REDGES[1:]); NB = len(rc)


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
            sp = []; pos = np.empty((n, 3)); tan = np.zeros((n, 3))
            for i in range(n):
                t = fh.readline().split(); sp.append(t[0]); pos[i] = (float(t[1]), float(t[2]), float(t[3]))
                if len(t) >= 11:
                    tan[i] = (float(t[8]), float(t[9]), float(t[10]))
            out.append((sw, np.array(sp), pos, tan))
    return out


INPUTS = sorted(glob.glob("inputs/J*_s*.dat"))
if not INPUTS:
    sys.exit("run from runs/<geom>/star<n>")
kv = read_kv(INPUTS[0]); NARM, F = int(kv["N"]), int(kv["n_arms"])
GEOM = os.path.basename(os.path.dirname(os.getcwd()))

groups = {}
def add(J, fr_list, tag):
    g = groups.setdefault(J, dict(r=[], psi=[], side=[], chi=[], run=[], nfr=0, runs=[]))
    rid = len(g["runs"])
    for sw, sp, pos, tan in fr_list:
        m = sp != "X"; s = sp[m]; p = pos[m]; u = tan[m]
        nu = np.linalg.norm(u, axis=1); ok = nu > 1e-9; u[ok] /= nu[ok][:, None]
        hel = (s != "C") & ok; arm = np.arange(len(s)) // NARM; sign = np.where(s == "R", 1, -1)
        H = np.where(hel)[0]; PH = p[H]; UH = u[H]; AH = arm[H]; SH = sign[H]
        for i0 in range(0, len(H), 400):
            i1 = min(len(H), i0 + 400)
            d = PH[None, :, :] - PH[i0:i1, None, :]; d2 = (d * d).sum(-1)              # r_ij = p_j - p_i
            pair = (d2 < R_MAX * R_MAX) & (AH[i0:i1, None] != AH[None, :]) & (np.arange(i0, i1)[:, None] < np.arange(len(H))[None, :])
            ii, jj = np.where(pair)
            if len(ii) == 0:
                continue
            dist = np.sqrt(d2[ii, jj]); rh = d[ii, jj] / dist[:, None]; ii += i0
            e1 = (rh * UH[ii]).sum(1); e2 = (rh * UH[jj]).sum(1); cij = (UH[ii] * UH[jj]).sum(1)
            num = (np.cross(UH[ii], UH[jj]) * rh).sum(1); den = cij - e1 * e2
            g["r"].append(dist); g["psi"].append(np.arctan2(num, den)); g["side"].append((np.abs(e1) < SIDE) & (np.abs(e2) < SIDE))
            g["chi"].append(np.where(SH[ii] == SH[jj], SH[ii], 0))                        # +1 R.R, -1 L.L, 0 R.L
            g["run"].append(np.full(len(ii), rid))
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


def mse(v): v = np.asarray(v); return (v.mean(), v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else (np.nan, np.nan)


res = {}
for J, g in sorted(groups.items()):
    r = np.concatenate(g["r"]); psi = np.concatenate(g["psi"]); side = np.concatenate(g["side"]); chi = np.concatenate(g["chi"]); run = np.concatenate(g["run"])
    psi_a = (psi + np.pi / 2) % np.pi - np.pi / 2                         # apolar twist in (-90, 90]
    tw = np.sin(psi) * np.cos(psi); sh = chi != 0
    T = np.full(NB, np.nan); Te = np.full(NB, np.nan); A = np.full(NB, np.nan); Ae = np.full(NB, np.nan); n_side = np.zeros(NB, int)
    b = np.clip(np.searchsorted(REDGES, r, side="right") - 1, 0, NB)
    for k in range(NB):
        mk = (b == k) & sh
        if mk.sum() > 5:
            T[k], Te[k] = mse(chi[mk] * tw[mk])
        ms = mk & side; n_side[k] = ms.sum()
        if n_side[k] > 5:
            A[k], Ae[k] = mse(np.degrees(chi[ms] * psi_a[ms]))
    def qfit(sel):
        """weighted fit psi = q r through the origin over r <= R_FIT for the selected side-by-side pairs"""
        num = den = 0.0
        for k in range(NB):
            if rc[k] > R_FIT:
                continue
            m = (b == k) & side & sel
            if m.sum() > 5:
                mu, se = mse(np.degrees(psi_a[m])); w = 1.0 / se ** 2
                num += w * rc[k] * mu; den += w * rc[k] ** 2
        return (num / den, 1.0 / np.sqrt(den)) if den > 0 else (np.nan, np.nan)
    q_rr, qe_rr = qfit(chi > 0); q_ll, qe_ll = qfit(chi < 0)
    mf = np.isfinite(A) & (rc <= R_FIT) & (Ae > 0)
    if mf.sum() >= 2:
        w = 1.0 / Ae[mf] ** 2; q = (w * rc[mf] * A[mf]).sum() / (w * rc[mf] ** 2).sum(); qe = 1.0 / np.sqrt((w * rc[mf] ** 2).sum())
    else:
        q = qe = np.nan
    close = side & (r < R_CLOSE)
    per_run = [np.mean(chi[m] * tw[m]) for k in range(len(g["runs"])) for m in [(run == k) & (r < R_CLOSE) & sh] if m.sum() > 20]
    T_run = (np.mean(per_run), np.std(per_run, ddof=1) / np.sqrt(len(per_run)) if len(per_run) > 1 else np.nan, len(per_run))
    res[J] = dict(T_run=T_run, per_run=per_run, T=T, Te=Te, A=A, Ae=Ae, n_side=n_side, q=q, qe=qe, q_rr=q_rr, qe_rr=qe_rr, q_ll=q_ll, qe_ll=qe_ll,
                  hist_rr=np.degrees(psi_a[close & (chi > 0)]), hist_ll=np.degrees(psi_a[close & (chi < 0)]), hist_rl=np.degrees(psi_a[close & (chi == 0)]),
                  n_close=int(close.sum()), T_close=mse(chi[(r < R_CLOSE) & sh] * tw[(r < R_CLOSE) & sh]),
                  T_rr=mse(tw[(r < R_CLOSE) & (chi > 0)]), T_ll=mse(tw[(r < R_CLOSE) & (chi < 0)]), T_rl=mse(tw[(r < R_CLOSE) & (chi == 0)]),
                  n_rr=int(((r < R_CLOSE) & (chi > 0)).sum()), n_ll=int(((r < R_CLOSE) & (chi < 0)).sum()), n_pairs=len(r), nfr=g["nfr"])

Js = sorted(res); cmap = plt.get_cmap("viridis"); col = {J: cmap(0.1 + 0.8 * i / max(1, len(Js) - 1)) for i, J in enumerate(Js)}
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 7.5, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
INK, MUTED = "#0b0b0b", "#898781"
fig, (a, b, c) = plt.subplots(1, 3, figsize=(15.0, 4.4))
a.axhline(0, color=MUTED, lw=0.8)
for J in Js:
    m = np.isfinite(res[J]["T"])
    a.fill_between(rc[m], (res[J]["T"] - res[J]["Te"])[m], (res[J]["T"] + res[J]["Te"])[m], color=col[J], alpha=0.15, lw=0)
    a.plot(rc[m], res[J]["T"][m], "-", color=col[J], lw=1.6, label=f"J = {J:g}")
a.set_xlabel("separation of the two helical residues r [a]"); a.set_ylabel(r"twist correlation $T_\chi(r) = \langle \chi \sin\psi \cos\psi \rangle$")
a.set_ylim(-0.12, 0.12); a.text(0.03, 0.05, "(a)  χ = +1 R·R, −1 L·L; 0 = no twist sense; bands = ±1 s.e.", transform=a.transAxes, fontsize=8, color=INK)
a.legend(loc="upper right", ncol=2, handlelength=1.5)
b.axhline(0, color=MUTED, lw=0.8)
for J in Js:
    m = np.isfinite(res[J]["A"])
    b.errorbar(rc[m], res[J]["A"][m], yerr=res[J]["Ae"][m], fmt="o-", ms=3.5, lw=1.0, color=col[J], ecolor=col[J], elinewidth=0.8, capsize=2,
               label=f"J = {J:g}: q = {res[J]['q']:+.2f} ± {res[J]['qe']:.2f} °/a" if np.isfinite(res[J]["q"]) else f"J = {J:g}")
b.set_xlabel("separation r [a]  (side-by-side same-handed pairs, |r̂·û| < %g)" % SIDE); b.set_ylabel(r"chirality-weighted apolar twist $\langle \chi\,\psi_a \rangle$ [deg]")
b.set_ylim(-25, 25); b.text(0.03, 0.05, f"(b)  cholesteric: ψ = q r for R·R (−q for L·L); fit over r ≤ {R_FIT:g} a", transform=b.transAxes, fontsize=8, color=INK)
b.legend(loc="upper right", handlelength=1.5, fontsize=7)
bins = np.arange(-90, 90.1, 10)
Jh = max(Js, key=lambda J: res[J]["n_close"])
for key, lab, colr, ls in (("hist_rr", "R·R", "#c0392f", "-"), ("hist_ll", "L·L", "#1c5cab", "-"), ("hist_rl", "R·L", MUTED, "--")):
    v = res[Jh][key]
    if len(v) > 10:
        h, _ = np.histogram(v, bins, density=True)
        c.step(bins, np.r_[h, h[-1]], where="post", color=colr, lw=1.4, ls=ls, label=f"{lab}  (n = {len(v)}, mean {v.mean():+.1f}°)")
c.axvline(0, color=MUTED, lw=0.8); c.set_xlim(-90, 90); c.set_xlabel(r"apolar twist $\psi_a$ of close side-by-side pairs (r < %g a) [deg]" % R_CLOSE); c.set_ylabel("probability density [1/deg]")
c.text(0.03, 0.05, f"(c)  J = {Jh:g}; a cholesteric tendency skews R·R and L·L to opposite sides", transform=c.transAxes, fontsize=8, color=INK)
c.legend(loc="upper right", handlelength=1.5)
fig.tight_layout(); fig.savefig(f"cholesteric_star{SUFFIX}.png", dpi=170); fig.savefig(f"cholesteric_star{SUFFIX}.pdf"); plt.close(fig)

print(f"{GEOM} star f = {F}: cholesteric (twist) order between helices of different arms, R.R (+) and L.L (-) weighted; frames = last half of each run (+ pivot control at J = 9)")
for J in sorted(res):
    print(f"   J = {J:g}: T_chi(r<{R_CLOSE:g}a) per run = " + ", ".join(f"{v:+.4f}" for v in res[J]["per_run"]) + f"   ({', '.join(groups[J]['runs'])})")
print("errors: +- = pair counting (optimistic: pairs of one frame and frames of one run are correlated); 'per run' = mean and s.e.m. over the runs' own values")
print(f"{'J':>4} {'frames':>6} {'T_chi(r<' + f'{R_CLOSE:g}' + ')':>10} {'+-':>6} {'per run':>8} {'+-':>6} {'n':>2} {'T_RR':>7} {'+-':>6} {'T_LL':>7} {'+-':>6} {'T_RL':>7} {'n_RR/n_LL':>10} {'q_chi[deg/a]':>12} {'+-':>5} {'q_RR':>6} {'q_LL':>6} {'pitch [a]':>18}")
for J in Js:
    r_ = res[J]; q, qe = r_["q"], r_["qe"]
    if np.isfinite(q) and abs(q) > 2 * qe:
        pitch = f"{360.0 / abs(q):.0f}"
    elif np.isfinite(qe):
        pitch = f"> {360.0 / (abs(q) + 2 * qe):.0f} (2σ bound)"
    else:
        pitch = "–"
    print(f"{J:>4g} {r_['nfr']:>6} {r_['T_close'][0]:>10.4f} {r_['T_close'][1]:>6.4f} {r_['T_run'][0]:>8.4f} {r_['T_run'][1]:>6.4f} {r_['T_run'][2]:>2} {r_['T_rr'][0]:>7.4f} {r_['T_rr'][1]:>6.4f} {r_['T_ll'][0]:>7.4f} {r_['T_ll'][1]:>6.4f} {r_['T_rl'][0]:>7.4f} "
          f"{r_['n_rr']:>4}/{r_['n_ll']:<5} {q:>12.2f} {qe:>5.2f} {r_['q_rr']:>6.2f} {r_['q_ll']:>6.2f} {pitch:>18}")
