#!/usr/bin/env python3
"""Time series of a periodic-box campaign (user 2026-09-24: "show me the time step of energy and helical fraction for
chains"): helicity and energies per residue against the CUMULATIVE production sweep, the restarted runs stitched to
their parents (restart_file -> the parent's last frame sets the offset; the parent's rows after that frame were not
continued and are dropped). Looks in <dir>/inputs and <dir>/prev_*/inputs.
  cd runs/t45/box500 && ~/venv/bin/python ../../box_series.py     -> series_box.png + table
"""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def read_kv(f):
    d = {}
    for l in open(f):
        l = l.split("#")[0]
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]; d[k] = v.split()[0]
    return d


def last_frame_sweep(conf):
    sw = -1
    try:
        with open(conf) as fh:
            for l in fh:
                if "sweep=" in l:
                    sw = int(l.split("sweep=")[1].split()[0])
    except OSError:
        pass
    return sw


def read_obs(f):
    rows = [[float(x) for x in l.split()] for l in open(f) if not l.startswith("#") and l.strip()]
    return np.array(rows) if rows else np.zeros((0, 19))


HERE = os.path.abspath(os.getcwd())
nodes = {}                                                    # conf realpath -> node
for inp in sorted(glob.glob("inputs/J*_s*.dat") + glob.glob("prev_*/inputs/J*_s*.dat")):
    kv = read_kv(inp); base = os.path.dirname(os.path.dirname(inp)) or "."
    J = float(re.match(r"J([\d.]+)_s", os.path.basename(inp)).group(1))
    prefix = os.path.join(base, kv["out_prefix"]); conf = os.path.realpath(prefix + "_conf.xyz")
    nodes[conf] = dict(J=J, name=os.path.basename(inp)[:-4], base=base, obs=prefix + "_obs.dat", conf=conf,
                       parent=os.path.realpath(kv["restart_file"]) if kv.get("restart_file") else None, N=int(kv["N"]) * int(kv["n_arms"]),
                       neq=int(kv["n_equil"]), offset=None, cut=None)
for n in nodes.values():                                      # the child's restart frame cuts the parent's history
    if n["parent"] in nodes:
        nodes[n["parent"]]["cut"] = last_frame_sweep(n["parent"])
def offset(n):
    if n["offset"] is None:
        n["offset"] = 0 if n["parent"] not in nodes else offset(nodes[n["parent"]]) + last_frame_sweep(n["parent"])
    return n["offset"]
for n in nodes.values():
    offset(n)

# keep, per J, only the lineage of the LIVE run (the one in inputs/): abandoned first attempts (e.g. prev_disp1) are dropped
live = {}
for n in nodes.values():
    if n["base"] == ".":
        chain = []; m = n
        while m is not None:
            chain.append(m["conf"]); m = nodes.get(m["parent"])
        live.setdefault(n["J"], set()).update(chain)
series = {}                                                   # J -> list of (sweep, rows, name)
for n in sorted(nodes.values(), key=lambda x: (x["J"], x["offset"])):
    if n["J"] in live and n["conf"] not in live[n["J"]]:
        continue
    r = read_obs(n["obs"])
    if len(r) == 0:
        continue
    if n["cut"] is not None:
        r = r[r[:, 0] <= n["cut"]]
    series.setdefault(n["J"], []).append((r[:, 0] + n["offset"], r, n["name"], n["N"]))
if not series:
    sys.exit("no rows yet")

Js = sorted(series); cmap = plt.get_cmap("viridis"); col = {J: cmap(0.1 + 0.8 * i / max(1, len(Js) - 1)) for i, J in enumerate(Js)}
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "legend.frameon": False, "legend.fontsize": 8, "axes.grid": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
INK, MUTED = "#0b0b0b", "#898781"
fig, ax = plt.subplots(2, 2, figsize=(12.0, 7.6))
COMP = ((4, "E_state"), (5, "E_bond"), (6, "E_bend"), (7, "E_nb"), (8, "E_twist"))
for J in Js:
    first = True
    for sw, r, name, N in series[J]:
        ax[0, 0].plot(sw, r[:, 3], "-", color=col[J], lw=1.2, label=f"J = {J:g}" if first else None)
        E = r[:, 4:9].sum(1) / N
        ax[0, 1].plot(sw, E, "-", color=col[J], lw=1.2, label=f"J = {J:g}" if first else None)
        ax[1, 0].plot(sw, r[:, 7] / N, "-", color=col[J], lw=1.2, label=f"J = {J:g}" if first else None)
        ax[1, 1].plot(sw, (r[:, 5] + r[:, 6] + r[:, 8]) / N, "-", color=col[J], lw=1.2, label=f"J = {J:g}" if first else None)
        for a in ax.flat:
            a.axvline(sw[0], color=MUTED, lw=0.6, ls=":")
        first = False
ax[0, 0].set_ylabel("helicity θ of all chains"); ax[0, 1].set_ylabel(r"total energy per residue [$k_BT$]")
ax[1, 0].set_ylabel(r"non-bonded energy per residue [$k_BT$]"); ax[1, 1].set_ylabel(r"bond + bend + twist energy per residue [$k_BT$]")
for k, a in enumerate(ax.flat):
    a.set_xlabel("cumulative production sweep"); a.text(0.02, 0.05, f"({'abcd'[k]})", transform=a.transAxes, fontsize=9, color=INK)
ax[0, 0].legend(loc="lower right", ncol=2)
ax[0, 0].text(0.02, 0.93, "dotted: restart / run boundary", transform=ax[0, 0].transAxes, fontsize=8, color=MUTED)
fig.suptitle(os.path.basename(HERE) + ": " + ", ".join(f"J = {J:g}: " + " → ".join(nm for _, _, nm, _ in series[J]) for J in Js), fontsize=8, color=MUTED, y=0.995)
fig.tight_layout(); fig.savefig("series_box.png", dpi=170); plt.close(fig)

print(f"{os.path.basename(HERE)}: series per J (cumulative production sweeps; runs stitched at their restart frames)")
for J in Js:
    tot = sum(len(r) for _, r, _, _ in series[J]); last_sw, last_r, _, N = series[J][-1]
    E = last_r[:, 4:9].sum(1) / N
    print(f"  J = {J:g}: {' -> '.join(nm for _, _, nm, _ in series[J])}; sweeps 0 .. {int(last_sw[-1])}, {tot} rows; theta first {series[J][0][1][0, 3]:.4f} last {last_r[-1, 3]:.4f}; "
          f"E/N first {series[J][0][1][:, 4:9].sum(1)[0] / N:.3f} last {E[-1]:.3f}; E_nb/N last {last_r[-1, 7] / N:.3f}")
