#!/usr/bin/env python3
import os; os.chdir(os.path.dirname(os.path.abspath(__file__)))
# helix-helix pairs (|i-j| > 1) inside the fit cutoff (3 a = 3.693 code units) per dumped frame -- the only pairs eps_s acts on
import glob, re, numpy as np
HF_LEN = 1.0   # code units per a; the archived prev_hflen1.2311 runs used 1.2311
RC = 3.0 * HF_LEN; out = {}
for f in sorted(glob.glob("out/J*_es*_s*_conf.xyz")):
    J = int(re.search(r"J(\d+)_", f).group(1)); L = open(f).read().split("\n"); i = 0
    while i < len(L) and L[i].strip():
        n = int(L[i]); rows = [l.split() for l in L[i + 2:i + 2 + n]]; i += 2 + n
        if len(rows) < n: break
        hel = np.array([r[0] != "C" for r in rows]); p = np.array([[float(x) for x in r[1:4]] for r in rows])
        d = np.linalg.norm(p[:, None] - p[None], axis=2); iu = np.triu_indices(n, 2)
        close = d[iu] < RC; hh = hel[iu[0]] & hel[iu[1]]
        out.setdefault(J, []).append(((close & hh).sum(), close.sum(), hel.sum()))
print("  J  frames  <helical residues>  <pairs inside cutoff>  <of which helix-helix>  frames with >=1 helix-helix contact")
for J in sorted(out):
    a = np.array(out[J]); print(f"  {J}  {len(a):6d}  {a[:,2].mean():18.1f}  {a[:,1].mean():21.2f}  {a[:,0].mean():22.3f}  {100*(a[:,0]>0).mean():10.1f} %")
