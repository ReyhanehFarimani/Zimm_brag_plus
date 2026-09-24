#!/usr/bin/env python3
"""Live page for the STAR J scan (user 2026-09-23: "make an artifact and show me the new results every 30 min").
Runs star_plot.py and star_relax.py in runs/<geom>/star<n>, then writes <out>/index.html + <out>/img/*.png for the
Artifact tool: status strip (finished / running / queued, sweep rate, ETA), the two figures, the per-J table that
star_plot.py prints, the control runs (ctrl_pivot, ctrl_allR, ctrl_linear50) and a dated log of the campaign.

  ~/venv/bin/python runs/star_page.py --out <dir> [--dir runs/t45/star50] [--no-refresh]
"""
import argparse
import datetime as dt
import glob
import html
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

LOG = [  # (time, text) -- the campaign log shown at the bottom, newest first
    ("12:15", "Box RESTARTED from the last dumped frame of run 1 (production sweep 6000, states + positions + registries, new seed; the reader reproduces the logged energies of that frame exactly) for 20k more sweeps. Run 1 (~6700 sweeps) is archived in prev_run1/. New input keys restart_file / restart_frame."),
    ("04:50", "Box relaunched with dense-system move steps: max_disp 0.3 (was 1.0), pivot angle 0.3 rad (was π), skin 2. The first launch accepted 3.5 % of position moves and 0.3 % of pivots and sat frozen: helicity 0.91 at J = 9 against 1.00, creeping up 0.002 per 100 sweeps. Its 250 production sweeps are archived in prev_disp1/."),
    ("03:02", "Periodic box launched: 500 free chains × 50 residues, box 63 a (residue density 0.1 = monomer density 0.7), J = 9, 8, 7, 6 on 4 workers. New code: `box` input with minimum-image non-bonded interactions; a cell grid for pivots; the neighbour list re-references a bead before an uncovered trial move (the all-pairs fallback was 95 % of the cost of large systems)."),
    ("02:45", "Harness on a dense periodic box exposed a PRE-EXISTING pivot bug: the pivot bead's own pair energies (its axis turns with the rotated bond) were never recomputed; fixed. In the dilute star this was 1 in 588 pivots by 0.002 kT (the mismatch seen earlier); in the box 16 kT. A second periodic-only piece: intra-tail pairs change through the minimum image when a chain is longer than half the box; counted now."),
    ("02:35", "User: J < 5 star runs stopped (inputs in inputs_dropped/), remaining star seeds (J5_s2, all J ≥ 5 seed 3) HELD; the freed cores go to the box."),
    ("02:55", "Bonded twist (bondedtwist_star.png): the backbone torsion of helical stretches is zero for R and for L within 1–2° at J ≥ 6.5 (no chiral supercoiling of the arms); the registry rotation is +25.3° (R·R) / −25.3° (L·L) with σ = 22°, exactly the twist term."),
    ("02:45", "Cholesteric order redone with R and L separated (χ-weighted; the first version averaged R·R with L·L and cancelled any signal by symmetry). Pair-counting errors now show a positive T_χ at J = 6.5, 7.5, 8 (0.02–0.035) but zero at J = 7 and 9 and an inconsistent R·R / L·L mirror check at J = 6: the frames of one run are correlated, so the run-to-run scatter is the honest error; undecided until more frames and the third seeds arrive."),
    ("02:30", "Cholesteric order (cholesteric_star.png) first version, R·R and L·L averaged together: zero by construction."),
    ("02:15", "Arm-handedness pair correlation (handpair_star.png): nearest-neighbour graft pairs and contacting pairs are uncorrelated at every J (P(same hand) = 0.47–0.54, all within 2σ of 1/2). Nematic order figure added (nematic_star.png)."),
    ("00:45", "Shell handedness added to the radial figure (panel d): |m| per shell against the independent-arm null of the same frames; with 2 frames per J it scatters on both sides of the null, no radial handedness segregation."),
    ("00:30", "Radial density added (density_star.png): helical arms decay with slope −1.5 to −1.7 between 8 and 40 a, between the rigid-rod star (−2) and the Daoud–Cotton flexible star (−4/3); arm tips sit at 0.56–0.59 of the rod length; the shells nearest the core are more helical than the rest of the arm (J = 6: 0.99 vs 0.88)."),
    ("23:03", "Batch stopped and relaunched with pivots (n_pivot = 50, one attempt per arm per sweep). The no-pivot partial data are archived in prev_nopivot."),
    ("22:58", "Control: the J = 9 star WITH pivots reaches helicity 0.995 by sweep 400 and 0.998 by 550 (1D theory 0.9997); the no-pivot production runs sat at 0.970."),
    ("22:53", "Linear 50-residue chains at J = 9 with the star protocol reach helicity 0.98 within 150–700 sweeps with or without pivots: the freeze is specific to the crowded star."),
    ("22:53", "Transition plot now averages only the LAST half of each run's production rows."),
    ("22:40", "Trajectory frames: the missing helicity sits in short coil gaps inside the arms (R·C·R and R·C..C·L), never at the graft or the tip; walls per arm decay only slowly."),
    ("22:20", "Helix-pair table: same- and opposite-handed contacts differ by < 0.02 kT wherever the table is trusted (r ≥ 1.9 a); no attraction anywhere (deepest value −0.03 kT). No mechanism for chirality transfer between arms."),
    ("21:14", "Scan launched: 13 J values × 3 seeds, 2000 + 30000 sweeps, 20 workers (no pivots)."),
]


def read_kv(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
    return d


def obs_rows(f):
    try:
        return [[float(x) for x in l.split()] for l in open(f) if not l.startswith("#") and l.strip()]
    except (OSError, ValueError):
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--dir", default=os.path.join(HERE, "t45", "star50"))
    ap.add_argument("--no-refresh", action="store_true")
    a = ap.parse_args()
    d = os.path.abspath(a.dir); geom = os.path.basename(os.path.dirname(d))
    inputs = sorted(glob.glob(os.path.join(d, "inputs", "J*_s*.dat")))
    if not inputs:
        sys.exit(f"no inputs in {d}")
    kv = read_kv(inputs[0]); NEQ, NSW, NARM, NARMS = int(kv["n_equil"]), int(kv["n_sweeps"]), int(kv["N"]), int(kv["n_arms"])
    CORE, NPIV = float(kv["core_radius"]), int(kv.get("n_pivot", 0))

    # ---- figures + table (star_plot.py prints the per-J table)
    table_txt = ""
    if not a.no_refresh:
        r = subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_plot.py")], cwd=d, capture_output=True, text=True, timeout=900)
        table_txt = r.stdout
        subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_relax.py")], cwd=d, capture_output=True, text=True, timeout=900)
        subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_density.py")], cwd=d, capture_output=True, text=True, timeout=900)
        subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_handpair.py")], cwd=d, capture_output=True, text=True, timeout=900)
        subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_nematic.py")], cwd=d, capture_output=True, text=True, timeout=900)
        subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_cholesteric.py")], cwd=d, capture_output=True, text=True, timeout=900)
        subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_bondedtwist.py")], cwd=d, capture_output=True, text=True, timeout=900)
    else:
        r = subprocess.run([PY, os.path.join(HERE, "star_plot.py")], cwd=d, capture_output=True, text=True, timeout=900); table_txt = r.stdout

    # ---- run status
    now = time.time(); rows = []; rates = []; n_fin = n_run = n_q = n_held = 0; slowest_eta = 0.0
    rate_J = {}; remaining = []; queued = []
    for f in inputs:
        b = os.path.basename(f)[:-4]; J = float(re.match(r"J([\d.]+)_s", b).group(1)); s = int(b.split("_s")[1])
        log = os.path.join(d, "logs", b + ".log"); lock = os.path.join(d, "logs", b + ".lock")
        done = os.path.isfile(log) and "summary" in open(log).read()
        ob = obs_rows(os.path.join(d, "out", b + "_obs.dat"))
        last = int(ob[-1][0]) if ob else -1
        if done:
            n_fin += 1; state = "finished"; prog = NSW
        elif os.path.isdir(lock) and subprocess.run(["pgrep", "-f", f"zimm inputs/{b}[.]dat"], capture_output=True).returncode != 0:
            n_held += 1; state = "held"; prog = 0                                     # lock placeholder, no process (HELD_SEEDS.txt)
        elif os.path.isdir(lock):
            n_run += 1; state = "running"; prog = NEQ + last + 1 if last >= 0 else 0     # 0 = still equilibrating (progress unknown)
            t0 = os.path.getmtime(lock); el = max(1.0, now - t0)
            if last >= 200:
                rate = el / prog; rates.append(rate); rate_J.setdefault(J, rate); rem = (NEQ + NSW - prog) * rate
            else:                                    # still equilibrating: the rate is at least el / NEQ, the remaining time at least this
                rate = el / NEQ; rate_J.setdefault(J, rate); rem = NSW * rate
            remaining.append(rem); slowest_eta = max(slowest_eta, rem)
        else:
            n_q += 1; state = "queued"; prog = 0; queued.append((s, -J))
        rows.append((J, s, state, prog, last))
    med_rate = sorted(rates)[len(rates) // 2] if rates else float("nan")
    # list-scheduling estimate of the whole scan: the queued runs (xargs order: seed-major, J descending) take the
    # slots as the running ones free them, each lasting (NEQ + NSW) x the rate measured for the same J (else the slowest known)
    slots = sorted(remaining) + [0.0] * max(0, 20 - len(remaining)); slots = sorted(slots)[:20] if slots else [0.0] * 20
    fallback = max(rate_J.values()) if rate_J else float("nan")
    for s_, negJ in sorted(queued):
        slots.sort(); slots[0] += (NEQ + NSW) * rate_J.get(-negJ, fallback)
    eta_batch = max(slots) if (queued or remaining) and fallback == fallback else float("nan")
    def hm(sec): return "–" if not sec or sec != sec else f"{sec / 3600:.1f} h"

    # ---- per-J table from star_plot.py
    trs = []
    for l in table_txt.splitlines():
        m = re.match(r"\s*([\d.]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)(.*)", l)
        if not m:
            continue
        J, n, th, tht, am, amt, u4, u4t, rg, enb, etw, rest = m.groups()
        part = re.search(r"partial: (\d+)k", rest)
        trs.append(f"<tr><td>{float(J):g}</td><td>{n}</td><td>{float(th):.3f}</td><td class='th'>{float(tht):.3f}</td>"
                   f"<td>{float(am):.3f}</td><td class='th'>{float(amt):.3f}</td><td>{float(u4):.2f}</td><td class='th'>{float(u4t):.2f}</td>"
                   f"<td>{float(rg):.0f}</td><td>{float(enb):.1f}</td><td>{float(etw):.0f}</td>"
                   f"<td class='mono'>{'≥ ' + part.group(1) + 'k, running' if part else 'finished'}</td></tr>")
    theory_line = next((l for l in table_txt.splitlines() if l.startswith("1D theory")), "")

    # ---- controls
    ctrl = []
    p = os.path.join(d, "ctrl_pivot", "out", "J9_piv50_obs.dat")
    if os.path.isfile(p) and obs_rows(p):
        ob = obs_rows(p); ctrl.append(("Star, J = 9, pivots on, all-coil start", f"sweep {int(ob[-1][0])}: helicity {ob[-1][3]:.4f}", "1D theory 0.9997; the no-pivot runs plateaued at 0.970"))
    p = os.path.join(d, "ctrl_allR", "out", "J9_allR_obs.dat")
    if os.path.isfile(p) and obs_rows(p):
        ob = obs_rows(p); ctrl.append(("Star, J = 9, no pivots, all-R start", f"sweep {int(ob[-1][0])}: helicity {ob[-1][3]:.4f}", "random initial registries re-nucleated it within one sweep; same slow healing as production"))
    lin = {}
    for f in sorted(glob.glob(os.path.join(d, "ctrl_linear50", "out", "*_obs.dat"))):
        ob = obs_rows(f)
        if len(ob) < 40:
            continue
        key = "with pivots" if "piv1" in f else "no pivots"
        th = [r_[3] for r_ in ob]; h = len(th) // 2
        first = next((int(r_[0]) for r_ in ob if r_[3] >= 0.98), None)
        lin.setdefault(key, []).append((sum(th[h:]) / len(th[h:]), first))
    for key, v in lin.items():
        ctrl.append((f"Linear chain of {NARM}, J = 9, {key}, {len(v)} seeds",
                     f"helicity {sum(x[0] for x in v) / len(v):.4f} over the last half",
                     "reaches 0.98 after " + ", ".join(str(x[1]) for x in v) + " sweeps"))

    # ---- periodic box (runs/box_gen.sh) as a second system on the page
    box_html = ""
    bd = os.path.join(HERE, "t45", "box500")
    if glob.glob(os.path.join(bd, "inputs", "J*_s*.dat")):
        bkv = read_kv(sorted(glob.glob(os.path.join(bd, "inputs", "J*_s*.dat")))[0])
        r = subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, "star_plot.py")], cwd=bd, capture_output=True, text=True, timeout=900)
        btab = []
        for l in r.stdout.splitlines():
            m = re.match(r"\s*([\d.]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)(.*)", l)
            if m:
                J, n, th, tht, am, amt, u4, u4t, rg, enb, etw, rest = m.groups(); part = re.search(r"partial: (\d+)k", rest)
                btab.append(f"<tr><td>{float(J):g}</td><td>{n}</td><td>{float(th):.3f}</td><td class='th'>{float(tht):.3f}</td><td>{float(am):.3f}</td><td class='th'>{float(amt):.3f}</td>"
                            f"<td>{float(u4):.2f}</td><td class='th'>{float(u4t):.2f}</td><td>{float(rg):.0f}</td><td>{float(enb):.1f}</td><td>{float(etw):.0f}</td><td class='mono'>{'≥ ' + part.group(1) + 'k, running' if part else 'finished'}</td></tr>")
        bstat = []
        for f in sorted(glob.glob(os.path.join(bd, "inputs", "J*_s*.dat"))):
            b = os.path.basename(f)[:-4]; ob = obs_rows(os.path.join(bd, "out", b + "_obs.dat")); last = int(ob[-1][0]) if ob else -1
            done = os.path.isfile(os.path.join(bd, "logs", b + ".log")) and "summary" in open(os.path.join(bd, "logs", b + ".log")).read()
            bstat.append(f'<span class="cell {"finished" if done else "running"}">{b.split("_")[0][1:]}<small>{"done" if done else (f"{int(bkv["n_equil"]) + last + 1} sw" if last >= 0 else "equil.")}</small></span>')
        fig = ""
        src = os.path.join(bd, "transitions_star.png")
        if os.path.isfile(src) and glob.glob(os.path.join(bd, "out", "J*_obs.dat")):
            shutil.copyfile(src, os.path.join(a.out, "img", "transitions_box.png"))
            fig = f'<figure class="fig"><img src="img/transitions_box.png?v={int(os.path.getmtime(src))}" alt="box transition"><figcaption><span>Same panels for the periodic box; blue = {bkv["n_arms"]} independent 1D chains with registries (handedness ceiling of independent chains ≈ 0.036).</span><span class="mono">t45/box500/transitions_star.png</span></figcaption></figure>'
        rho = int(bkv["n_arms"]) * int(bkv["N"]) / float(bkv["box"]) ** 3
        box_html = f'''
  <section>
    <h2>Periodic box: {bkv["n_arms"]} free chains × {bkv["N"]} residues</h2>
    <p class="muted" style="margin:0 0 8px">Cubic box of {float(bkv["box"]):.1f} a, residue density {rho:.3f} a⁻³ (monomer density {7 * rho:.2f} with 7 monomers per residue), minimum-image non-bonded interactions,
    one pivot and one flip per chain per sweep, {int(bkv["n_equil"])} + {int(bkv["n_sweeps"]) // 1000}k sweeps, J = 9, 8, 7, 6, one seed each; max_disp {float(bkv.get("max_disp", 1)):g}, pivot angle {float(bkv.get("max_rot", 3.14)):g} rad.{" Restarted 12:15 from the last frame of the previous run (its sweep 6000); sweeps count from 0 again." if bkv.get("restart_file") else " Relaunched 04:50."}</p>
    <div class="grid">{"".join(bstat)}</div>
    {fig}
    <div class="tbl"><table><thead><tr><th>J</th><th>seeds</th><th>θ</th><th class="th">theory</th><th>⟨|m|⟩</th><th class="th">theory</th><th>U₄</th><th class="th">theory</th><th>⟨R_g²⟩ per chain [a²]</th><th>⟨E_nb⟩</th><th>⟨E_twist⟩</th><th>status</th></tr></thead>
    <tbody>{"".join(btab) or "<tr><td colspan='12' class='muted'>no run has 20 production rows yet</td></tr>"}</tbody></table></div>
  </section>'''

    os.makedirs(os.path.join(a.out, "img"), exist_ok=True)
    # ---- snapshots (tools/render_box.py, tools/ovito_preset.py): every png in <dir>/snapshots of the box and the star
    snaps = []
    for label, sd in (("periodic box", os.path.join(HERE, "t45", "box500", "snapshots")), ("star", os.path.join(d, "snapshots"))):
        for f in sorted(glob.glob(os.path.join(sd, "*.jpg")) + glob.glob(os.path.join(sd, "*.png"))):
            name = os.path.basename(f)[:-4]; ext = f[-4:]; shutil.copyfile(f, os.path.join(a.out, "img", "snap_" + name + ext))
            m = re.match(r"J([\d.]+)_(box|slab|star)_sweep(\d+)", name)
            cap = (f"{label}, J = {m.group(1)}, production sweep {m.group(3)}" + (", slab of 3 a through the middle" if m.group(2) == "slab" else (", whole box, positions wrapped" if m.group(2) == "box" else ""))) if m else name
            snaps.append(f'<figure class="snap"><img src="img/snap_{name}{ext}?v={int(os.path.getmtime(f))}" alt="{html.escape(cap)}"><figcaption>{cap}</figcaption></figure>')
    snap_html = (f'''
  <section>
    <h2>Snapshots</h2>
    <p class="muted" style="margin:0 0 8px">Rendered with the house OVITO preset: right-handed helices red, left-handed amber, coil blue, the registry direction as a dark arrow on every rod (tools/render_box.py, tools/ovito_preset.py).</p>
    <div class="snaps">{"".join(snaps)}</div>
  </section>''' if snaps else "")

    # ---- assemble
    os.makedirs(os.path.join(a.out, "img"), exist_ok=True)
    figs = []
    for name, cap in (("transitions_star", "Helicity θ, χ_θ, handedness ⟨|m|⟩ and Binder U₄ of the whole star against J (top); ⟨R_g²⟩, ⟨E_nb⟩, ⟨E_twist⟩ (bottom). Blue: independent arms, 1D theory with registries; dashed: one arm alone; grey: no registries. Filled = finished, hollow = running; averages over the last half of each run."),
                      ("relax_star", "Relaxation check: block means (250 sweeps) of the helicity and of the signed handedness against production sweep for every run with data; dotted = the 1D theory per J; solid / dashed / dotted lines = seed 1 / 2 / 3."),
                      ("density_star", "Radial structure from the dumped frames (last half of each run; J = 9 also uses the pivot control). (a) Residue number density around the core centre per J, against the rigid-rod star f/(4π b r²) in absolute terms and the flexible-arm scalings (Daoud–Cotton good solvent ∝ r⁻⁴ᐟ³, ideal arms ∝ r⁻¹) drawn as slope guides through the softest MC profile at 10 a. (b) Mean radial position of residue i along the arm against the rod line. (c) Local helicity per shell: the helix fraction as a function of distance from the core. (d) Local handedness: |n_R − n_L| / (n_R + n_L) of the helical residues in each shell and frame, averaged over frames (solid), against the same frames with every arm mirrored at random (dashed) — solid above dashed would mean arms of like handedness sit together. One residue is a block of 7 fine monomers."),
                      ("handpair_star", "Pair correlation of the arm handedness. (a) ⟨s_k s_l⟩ of the arm handedness signs against the angular separation of the two graft sites (nearest neighbours ≈ 30°), bars = 1/√n_pairs, 0 = independent arms. (b) Probability that two arms share a hand when they are in helix–helix contact in that frame (< 3 a) against pairs that are not; 1/2 = independent arms; the small numbers are the counts of contacting pairs."),
                      ("nematic_star", "Orientational order of the helical rods. (a) Radial order ⟨P₂(û·r̂)⟩ per shell, helical residues solid, all residues dashed: 1 = rigid radial rods, 0 = isotropic. (b) Local nematic order ⟨P₂(û_i·û_j)⟩ between helices of different arms closer than 5 a, per shell of the pair midpoint: 1 = parallel bundles. (c) Mean sine of the crossing angle ψ of contacting helix pairs (< 3 a), same-handed and opposite-handed separately: 0 = no preferred twist sense between neighbouring arms; counts of same/opposite pairs above."),
                      ("cholesteric_star", "Cholesteric (twist) order between helices of different arms, R and L treated separately: by mirror symmetry R·R pairs twist opposite to L·L pairs and the star has as many L as R arms, so every twist quantity is weighted by the pair chirality χ = +1 (R·R), −1 (L·L); R·L pairs carry no signal. (a) Twist correlation T_χ(r) = ⟨χ sin ψ cos ψ⟩ against the separation of the two residues; 0 = no preferred twist sense, bands ± 1 s.e. by pair counting (optimistic: pairs of one frame and frames of one run are correlated). (b) Chirality-weighted apolar twist angle of side-by-side pairs against separation with a straight-line fit through the origin over r ≤ 8 a: slope q = twist rate of R·R pairs (L·L: −q), pitch = 360°/|q|. (c) Distributions of the apolar twist angle of close side-by-side pairs (r < 4 a) for R·R, L·L and R·L at the J with the most pairs; a cholesteric tendency skews R·R and L·L to opposite sides."),
                      ("bondedtwist_star", "Bonded (intra-arm) twist, R and L separately, the counterpart of the inter-arm cholesteric figure. (a) Mean backbone torsion of four consecutive residues for R helical, L helical and coil stretches, and the χ-weighted combination; the bend potential has no dihedral term, so a nonzero mirror-image pair would be emergent supercoiling of the arm. (b) Registry rotation between consecutive same-handed residues: the twist term imposes ±25.3° with σ = 22° at θ₀ = 45°, a check of the term. (c) Torsion distributions at the J with the most helical stretches.")):
        src = os.path.join(d, name + ".png")
        if os.path.isfile(src):
            shutil.copyfile(src, os.path.join(a.out, "img", name + ".png")); v = int(os.path.getmtime(src))
            figs.append(f'<figure class="fig"><img src="img/{name}.png?v={v}" alt="{html.escape(cap)}"><figcaption><span>{cap}</span><span class="mono">{geom}/star{NARMS}/{name}.png</span></figcaption></figure>')
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    grid = "".join(f'<span class="cell {st}" title="J = {J:g}, seed {s}: {st}, {(str(prog) + " sweeps") if prog else ("equilibrating" if st == "running" else ("held for later" if st == "held" else ""))}">{J:g}<small>s{s}</small></span>' for J, s, st, prog, last in rows)
    ctrl_html = "".join(f"<tr><td>{html.escape(k)}</td><td class='mono'>{html.escape(v)}</td><td class='muted'>{html.escape(n)}</td></tr>" for k, v, n in ctrl) or "<tr><td colspan='3' class='muted'>no control data yet</td></tr>"
    log_html = "".join(f"<li><span class='mono'>{t}</span><span>{html.escape(x)}</span></li>" for t, x in LOG)
    css = """
:root { --bg:#f6f5f1; --ink:#1f2023; --muted:#6a6d75; --rule:#dad8d1; --panel:#ffffff; --accent:#7f1310; --ok:#2f6b3a; --okbg:#e3efe2; --runbg:#f5e4e2; --qbg:#ecebe6; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme:dark; --bg:#17181b; --ink:#e8e6e0; --muted:#9a9da6; --rule:#33353b; --panel:#1f2024; --accent:#e07a5f; --ok:#8fcf98; --okbg:#1f3323; --runbg:#3a2320; --qbg:#26272b; } }
:root[data-theme="dark"] { color-scheme:dark; --bg:#17181b; --ink:#e8e6e0; --muted:#9a9da6; --rule:#33353b; --panel:#1f2024; --accent:#e07a5f; --ok:#8fcf98; --okbg:#1f3323; --runbg:#3a2320; --qbg:#26272b; }
body { background:var(--bg); color:var(--ink); font-family:"Source Sans 3","Helvetica Neue",Arial,sans-serif; font-size:15px; line-height:1.5; padding-inline:16px; padding-block:24px 48px; }
.wrap { max-width:1120px; margin:0 auto; display:flex; flex-direction:column; gap:24px; }
h1 { font-family:"Source Serif 4",Georgia,serif; font-weight:600; font-size:1.9rem; margin:0; text-wrap:balance; }
h2 { font-family:"Source Serif 4",Georgia,serif; font-weight:600; font-size:1.2rem; margin:0 0 8px; }
.sub { color:var(--muted); margin:6px 0 0; max-width:78ch; }
.status { display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:14px 22px; border-block:1px solid var(--rule); padding-block:14px; }
.status .k { font-size:.72rem; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); }
.status .v { font-family:"JetBrains Mono",Menlo,monospace; font-size:1.05rem; font-variant-numeric:tabular-nums; }
.grid { display:flex; flex-wrap:wrap; gap:4px; }
.cell { font-family:"JetBrains Mono",Menlo,monospace; font-size:.74rem; padding:2px 6px; border:1px solid var(--rule); background:var(--qbg); color:var(--muted); }
.cell small { margin-left:3px; opacity:.7; }
.cell.running { background:var(--runbg); color:var(--accent); border-color:var(--accent); }
.cell.finished { background:var(--okbg); color:var(--ok); border-color:var(--ok); }
.cell.held { border-style:dashed; }
.fig { margin:0; background:var(--panel); border:1px solid var(--rule); padding:10px; display:flex; flex-direction:column; gap:8px; }
.fig img { width:100%; max-width:100%; height:auto; display:block; }
figcaption { display:flex; justify-content:space-between; gap:12px; flex-wrap:wrap; font-size:.85rem; color:var(--muted); }
.mono { font-family:"JetBrains Mono",Menlo,monospace; font-size:.78rem; }
.muted { color:var(--muted); }
.tbl { overflow-x:auto; }
table { border-collapse:collapse; width:100%; font-variant-numeric:tabular-nums; font-size:.9rem; }
th, td { text-align:right; padding:5px 9px; border-bottom:1px solid var(--rule); white-space:nowrap; }
th { font-size:.72rem; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); font-weight:600; }
td.th { color:var(--muted); }
.ctrl td:first-child, .ctrl td:last-child { text-align:left; white-space:normal; }
.ctrl td:first-child { min-width:16ch; }
.snaps { display:grid; grid-template-columns:repeat(auto-fit, minmax(300px, 1fr)); gap:12px; }
.snap { margin:0; background:var(--panel); border:1px solid var(--rule); padding:8px; display:flex; flex-direction:column; gap:6px; }
.snap img { width:100%; max-width:100%; height:auto; display:block; }
.snap figcaption { font-size:.82rem; color:var(--muted); }
ul.log { list-style:none; padding:0; margin:0; display:flex; flex-direction:column; gap:6px; }
ul.log li { display:grid; grid-template-columns:5ch 1fr; gap:12px; font-size:.9rem; }
"""
    page = f"""<title>Star50 J Scan</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,500;8..60,600&family=Source+Sans+3:wght@400;600&family=JetBrains+Mono:wght@400;500&display=swap">
<style>{css}</style>
<div class="wrap">
  <header>
    <h1>Star50 J Scan</h1>
    <p class="sub">{NARMS} arms × {NARM} residues on a core of radius {CORE:g} a (θ₀ = {geom[1:]}°), six-argument helix–helix table plus registry twist,
    E<sub>HH</sub> = −J, E<sub>RL</sub> = +J, E<sub>CH</sub> = 0. J = {min(r[0] for r in rows):g} … {max(r[0] for r in rows):g} ({len(set(r[0] for r in rows))} values, seeds {", ".join(str(x) for x in sorted(set(r[1] for r in rows if r[2] != "held")))} running or done, seed 3 held), {NEQ} + {NSW // 1000}k sweeps per run,
    {NPIV} pivot attempts per sweep (one per arm). J &lt; 5 was dropped at 02:35 (inputs_dropped/). Blue curves: the 1D transfer matrix with the registries integrated out for {NARMS}
    independent arms; for independent arms the star's handedness cannot exceed 0.112, so anything above it would be chirality transfer between arms.</p>
  </header>
  <div class="status">
    <div><div class="k">Finished</div><div class="v">{n_fin} / {len(inputs)}</div></div>
    <div><div class="k">Running</div><div class="v">{n_run}</div></div>
    <div><div class="k">Queued / held</div><div class="v">{n_q} / {n_held}</div></div>
    <div><div class="k">Sweep time (median)</div><div class="v">{f"{med_rate:.2f} s" if med_rate == med_rate else "–"}</div></div>
    <div><div class="k">Running runs end within</div><div class="v">{hm(slowest_eta)}</div></div>
    <div><div class="k">Scan ends in (estimate)</div><div class="v">{hm(eta_batch)}</div></div>
    <div><div class="k">Last update</div><div class="v">{stamp}</div></div>
  </div>
  <div class="grid">{grid}</div>
  {"".join(figs)}
  <section>
    <h2>Averages over the last half of each run</h2>
    <p class="muted" style="margin:0 0 8px">{html.escape(theory_line)}. Grey columns: independent arms with registries.</p>
    <div class="tbl"><table><thead><tr><th>J</th><th>seeds</th><th>θ</th><th class="th">theory</th><th>⟨|m|⟩</th><th class="th">theory</th><th>U₄</th><th class="th">theory</th><th>⟨R_g²⟩ [a²]</th><th>⟨E_nb⟩</th><th>⟨E_twist⟩</th><th>status</th></tr></thead>
    <tbody>{"".join(trs) or "<tr><td colspan='12' class='muted'>no run has 20 production rows yet (the first 2000 sweeps are equilibration)</td></tr>"}</tbody></table></div>
  </section>
  {box_html}
  {snap_html}
  <section>
    <h2>Control runs</h2>
    <div class="tbl"><table class="ctrl"><tbody>{ctrl_html}</tbody></table></div>
  </section>
  <section>
    <h2>Log</h2>
    <ul class="log">{log_html}</ul>
  </section>
</div>
"""
    open(os.path.join(a.out, "index.html"), "w").write(page)
    print(f"{a.out}/index.html  finished {n_fin}, running {n_run}, queued {n_q}, median {med_rate:.2f} s/sweep, wave ends in {hm(slowest_eta)}, scan in {hm(eta_batch)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
