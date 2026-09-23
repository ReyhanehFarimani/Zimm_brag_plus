#!/usr/bin/env python3
"""Live page for the J scans after the k = 2K / 1-4-exclusion fixes (user 2026-09-23: "plot the transition
curves for t45 and t100 while running").  Runs the two existing figure scripts on the FINISHED runs
(t100: jscan_eps/plot_transitions.py -> transitions.png; t45: jscan_sigc/analyze.py -> transitions_t45.png,
eps_t45.png), then writes <out>/index.html + <out>/img/*.png for the Artifact tool.

  python runs/jscan_page.py --out <dir> [--no-refresh]
"""
import argparse
import datetime as dt
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
SCANS = {"t45": dict(d=os.path.join(HERE, "t45", "jscan_db"), script=os.path.join(HERE, "jscan_db_plot.py"),
                    figs=[("transitions_db", "θ₀ = 45°: helicity, χ_θ, ⟨|m|⟩, Binder U₄; ⟨R_g²⟩, ⟨E_nb⟩, ⟨E_twist⟩ against J")]),
         "t100": dict(d=os.path.join(HERE, "t100", "jscan_db"), script=os.path.join(HERE, "jscan_db_plot.py"),
                     figs=[("transitions_db", "θ₀ = 100°: same panels")])}


def counts(d):
    inputs = len(glob.glob(os.path.join(d, "inputs", "J*.dat")))
    logs = glob.glob(os.path.join(d, "logs", "J*_db_s*.log"))
    done = sum(1 for l in logs if "summary" in open(l).read())
    return done, inputs


def refresh():
    for camp, s in SCANS.items():
        if not glob.glob(os.path.join(s["d"], "out", "J*_obs.dat")):
            continue
        subprocess.run(["nice", "-n", "10", PY, s["script"]], cwd=s["d"], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=900)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-refresh", action="store_true")
    a = ap.parse_args()
    if not a.no_refresh:
        refresh()
    os.makedirs(os.path.join(a.out, "img"), exist_ok=True)
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    status, figs = [], []
    all_done = True
    for camp, s in SCANS.items():
        done, n = counts(s["d"])
        all_done &= (done >= n) if n > 0 else True                 # a scan without inputs is not part of the campaign
        status.append(f'<div><div class="k">{camp} runs finished</div><div class="v">{f"{done} / {n}" if n else "not started"}</div></div>')
        for name, cap in s["figs"]:
            p = os.path.join(s["d"], name + ".png")
            if os.path.isfile(p) and glob.glob(os.path.join(s["d"], "out", "J*_obs.dat")):
                shutil.copyfile(p, os.path.join(a.out, "img", f"{name}_{camp}.png"))
                figs.append(f'<figure class="fig"><img src="img/{name}_{camp}.png?v={int(os.path.getmtime(p))}" alt="{cap}">'
                            f'<figcaption><span>{cap}</span><span class="mono">{camp}/jscan_db/{name} · filled = finished, hollow = still running</span></figcaption></figure>')
            else:
                figs.append(f'<figure class="fig"><div class="missing">{cap}: no finished runs yet</div></figure>')
    page = f"""<title>J Scan Transitions</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,500;8..60,600&family=Source+Sans+3:wght@400;600&family=JetBrains+Mono:wght@400;500&display=swap">
<style>
:root {{ --bg:#f6f5f1; --ink:#1f2023; --muted:#6a6d75; --rule:#dad8d1; --panel:#ffffff; --accent:#7f1310; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#17181b; --ink:#e8e6e0; --muted:#9a9da6; --rule:#33353b; --panel:#1f2024; --accent:#e07a5f; }} }}
:root[data-theme="dark"] {{ --bg:#17181b; --ink:#e8e6e0; --muted:#9a9da6; --rule:#33353b; --panel:#1f2024; --accent:#e07a5f; }}
body {{ background:var(--bg); color:var(--ink); font-family:"Source Sans 3","Helvetica Neue",Arial,sans-serif; font-size:15px; line-height:1.5; padding-inline:16px; padding-block:24px 48px; }}
.wrap {{ max-width:1120px; margin:0 auto; display:flex; flex-direction:column; gap:24px; }}
h1 {{ font-family:"Source Serif 4",Georgia,serif; font-weight:600; font-size:1.9rem; margin:0; text-wrap:balance; }}
.sub {{ color:var(--muted); margin:6px 0 0; max-width:70ch; }}
.status {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:14px 22px; border-block:1px solid var(--rule); padding-block:14px; }}
.status .k {{ font-size:.72rem; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); }}
.status .v {{ font-family:"JetBrains Mono",Menlo,monospace; font-size:1.05rem; font-variant-numeric:tabular-nums; }}
.fig {{ margin:0; background:var(--panel); border:1px solid var(--rule); padding:10px; display:flex; flex-direction:column; gap:8px; }}
.fig img {{ width:100%; max-width:100%; height:auto; display:block; }}
figcaption {{ display:flex; justify-content:space-between; gap:12px; flex-wrap:wrap; font-size:.85rem; color:var(--muted); }}
.mono {{ font-family:"JetBrains Mono",Menlo,monospace; font-size:.78rem; }}
.missing {{ color:var(--muted); padding:40px; text-align:center; }}
</style>
<div class="wrap">
  <header>
    <h1>J Scan Transitions</h1>
    <p class="sub">Zimm–Bragg chain, N = 200, E<sub>HH</sub> = −J, E<sub>RL</sub> = +J, E<sub>CH</sub> = 0, with the
    six-argument helix–helix TABLE (r, β<sub>A</sub>, β<sub>B</sub>, ψ, α<sub>A</sub>, α<sub>B</sub>) and the registry twist
    terms (R·R / L·L: Δα₀ = ±25.3°, κ = 6.6 k<sub>B</sub>T/rad²; R·L: 180°, κ = 4.9 at θ₀ = 45°) — no ε<sub>s</sub> series.
    J = 5.5 … 7.5 in steps of 0.25 across the crossovers with three seeds, plus 1, 2, 3, 3.75 … 5.25, 8, 9, 10 with one seed;
    1 M sweeps (exact 1D crossover J*(N = 200) = 4.68 at θ₀ = 45°; the explicit registries shift it up by ≈ 2). Bonded constants k = 2K of the fits, non-bonded
    pairs from 1-4 (both fixed 2026-09-23); restarted 19:31 with the R·L twist term. Blue curves: exact 1D chain without non-bonded interactions, N = 100 … 1600.
    Redrawn from the finished runs.</p>
  </header>
  <div class="status">{"".join(status)}<div><div class="k">Last update</div><div class="v">{now}</div></div>
    <div><div class="k">Campaign</div><div class="v">{"FINISHED" if all_done else "running"}</div></div></div>
  {"".join(figs)}
</div>
"""
    open(os.path.join(a.out, "index.html"), "w").write(page)
    d45, n45 = counts(SCANS["t45"]["d"]); d100, n100 = counts(SCANS["t100"]["d"])
    print(f"{a.out}/index.html  (t45 {d45}/{n45}, t100 {d100}/{n100}{'; FINISHED' if all_done else ''})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
