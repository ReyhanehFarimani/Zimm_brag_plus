#!/usr/bin/env python3
"""jscan_eps summary: what does the chiral amplitude eps_s change?  FINISHED runs only (log has a summary).
Per run: block-averaged means (samples are serially correlated).  Per (J, eps_s) cell: mean over seeds, error =
seed scatter when >= 2 seeds, else the single run's block error.  For every observable the eps_s dependence is
reported as the difference to the eps_s = 0 cell of the same J, in units of the combined error (z).
  m = (n_R - n_L) / N  (handedness excess);  contact = fraction of samples with E_nb != 0."""
import glob, os, re
import numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
N = 200

def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def blk(x, nb=20):
    b = x[: len(x) // nb * nb].reshape(nb, -1).mean(1); return x.mean(), b.std(ddof=1) / np.sqrt(nb)

OBS = ["helicity", "absm", "m2", "E_nb", "contact", "Rg2", "E_state"]
runs = {}
for log in sorted(glob.glob("logs/J*_es*_s*.log")):
    if "summary" not in open(log).read(): continue
    g_ = re.search(r"J([\d.]+)_es(\d+)_s(\d+)", log).groups(); J, es, s = float(g_[0]), int(g_[1]), int(g_[2])
    o = read_obs("out/" + os.path.basename(log)[:-4] + "_obs.dat")
    m = (o["n_R"] - o["n_L"]) / N
    o.update(absm=np.abs(m), m2=m * m, contact=(o["E_nb"] != 0).astype(float))
    runs.setdefault((J, es), []).append({k: blk(o[k]) for k in OBS})

cell = {}
for key, rs in runs.items():
    cell[key] = {}
    for k in OBS:
        v = np.array([r[k][0] for r in rs]); e = np.array([r[k][1] for r in rs])
        cell[key][k] = (v.mean(), v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else e[0], len(v))
Js = sorted({k[0] for k in cell}); ESs = sorted({k[1] for k in cell})
print(f"finished runs: {sum(len(r) for r in runs.values())}   cells: {len(cell)}   seeds per cell: "
      f"{min(c['helicity'][2] for c in cell.values())}-{max(c['helicity'][2] for c in cell.values())}")
for k in OBS:
    print(f"\n== {k}:  value +- err   [z vs eps_s = 0 of the same J]")
    print("   J  " + "".join(f"{'eps_s=' + str(e):>26s}" for e in ESs))
    for J in Js:
        row = f"  {J:3g} "
        for es in ESs:
            if (J, es) not in cell: row += f"{'--':>26s}"; continue
            v, e, _ = cell[(J, es)][k]; txt = f"{v:.4g} +- {e:.2g}"
            if es != 0 and (J, 0) in cell:
                v0, e0, _ = cell[(J, 0)][k]; d = np.hypot(e, e0)
                txt += f" [{(v - v0) / d:+.1f}]" if d > 0 else " [ 0 ]"
            row += f"{txt:>26s}"
        print(row)
# one number per observable: is there ANY eps_s dependence beyond noise?
print("\n== all (J, eps_s > 0) cells vs their eps_s = 0 cell: RMS z (1 = pure noise) and largest |z|")
for k in OBS:
    z = []
    for (J, es), c in cell.items():
        if es == 0 or (J, 0) not in cell: continue
        d = np.hypot(c[k][1], cell[(J, 0)][k][1])
        if d > 0: z.append((c[k][0] - cell[(J, 0)][k][0]) / d)
    z = np.array(z)
    if len(z): print(f"  {k:9s} RMS z = {np.sqrt((z * z).mean()):5.2f}   max |z| = {np.abs(z).max():5.2f}   (n = {len(z)})")
