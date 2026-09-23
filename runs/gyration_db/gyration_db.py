#!/usr/bin/env python3
"""Gyration-tensor eigenvalues of the frozen-pattern N = 40 chains run by gen_gyration_db.sh:
lambda_1 >= lambda_2 >= lambda_3 per dumped frame, averaged after dropping the first 20 % of the
frames; the error of one run is the standard error over 10 time blocks, seeds are combined as
independent runs.  Reported per pattern for the tabulated potential (db) and the radial fit
(fit): the eigenvalues, Rg^2 = sum lambda_i, the shape ratios lambda_2/lambda_1 and
lambda_3/lambda_1, the run length in integrated autocorrelation times of Rg^2, and the
acceptance of pivot and twist moves.  The allL / allR pair is the mirror check of the table.

  ~/.conda/envs/sim_analysis/bin/python runs/gyration_db/gyration_db.py [--fig]
"""
import glob
import os
import re
import sys

import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
PATTERNS = ["allL", "allR", "halfL_C", "halfL_R"]
POTS = ["db", "fit"]


def read_frames(f):
    frames = []
    with open(f) as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            n = int(head)
            fh.readline()
            rows = [fh.readline().split() for _ in range(n)]
            if len(rows[-1]) < 4:
                break                                       # half-written last frame
            frames.append([[float(v) for v in r[1:4]] for r in rows])
            species = [r[0] for r in rows]
    return np.array(frames), np.array(species)


def contacts(P, species, rmax):
    """mean number per frame of helix-helix pairs j >= i + 2 (the non-bonded exclusion of the code)
    closer than rmax -- how often the pair table is actually consulted."""
    h = np.where(species != "C")[0]
    if len(h) < 3:
        return 0.0
    ii, jj = np.triu_indices(len(h), 2)
    d = np.linalg.norm(P[:, h[ii]] - P[:, h[jj]], axis=2)
    return float((d < rmax).sum(1).mean())


def read_obs(f):
    lines = open(f).read().splitlines()
    hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines
                     if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}


def tau_int(x):
    x = np.asarray(x, float); x = x - x.mean(); n = len(x)
    if n < 16 or x.std() == 0:
        return np.nan
    f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); ac /= ac[0]
    t = 0.5
    for M in range(1, n // 2):
        t += ac[M]
        if M > 6 * t:
            break
    return t


def summary_value(log, key):
    m = re.search(r"# " + re.escape(key) + r"\s*=\s*([-0-9.eE+]+)", open(log).read())
    return float(m.group(1)) if m else np.nan


def collect():
    runs = {}
    for f in sorted(glob.glob("out/*_conf.xyz")):
        m = re.match(r"out/(\w+?)_(db|fit)_s(\d+)_conf", f)
        if not m:
            continue
        pat, pot, seed = m.group(1), m.group(2), int(m.group(3))
        log = f"logs/{pat}_{pot}_s{seed}.log"
        if not os.path.exists(log) or "summary" not in open(log).read():
            continue                                        # unfinished
        P, species = read_frames(f)
        if len(P) < 50:
            continue
        P = P[len(P) // 5:]
        nc3, nc2 = contacts(P, species, 3.0), contacts(P, species, 2.0)
        c = P - P.mean(1, keepdims=True)
        T = np.einsum("fia,fib->fab", c, c) / P.shape[1]
        lam = np.sort(np.linalg.eigvalsh(T), axis=1)[:, ::-1]
        nb = 10
        blk = lam[: len(lam) // nb * nb].reshape(nb, -1, 3).mean(1)
        o = read_obs(f.replace("_conf.xyz", "_obs.dat"))
        dt = o["sweep"][1] - o["sweep"][0]
        runs.setdefault((pat, pot), []).append(dict(
            lam=lam.mean(0), err=blk.std(0, ddof=1) / np.sqrt(nb), nfr=len(lam),
            tau=tau_int(o["Rg2"]) * dt, nsw=o["sweep"][-1] - o["sweep"][0],
            acc_piv=summary_value(log, "acc_pivot"), acc_tw=summary_value(log, "acc_twist"),
            enb=summary_value(log, "<E_nb>"), nc3=nc3, nc2=nc2))
    res = {}
    for key, rs in runs.items():
        n = len(rs)
        res[key] = dict(lam=np.mean([r["lam"] for r in rs], 0),
                        err=np.sqrt(np.sum([r["err"] ** 2 for r in rs], 0)) / n,
                        n=n, tau=np.nanmean([r["tau"] for r in rs]), nsw=rs[0]["nsw"],
                        acc_piv=np.mean([r["acc_piv"] for r in rs]),
                        acc_tw=np.nanmean([r["acc_tw"] for r in rs]),
                        enb=np.mean([r["enb"] for r in rs]),
                        nc3=np.mean([r["nc3"] for r in rs]), nc2=np.mean([r["nc2"] for r in rs]))
    return res


def report(res):
    print(f"{'pattern':>8} {'pot':>4} {'seeds':>5} {'lam1':>8} {'+-':>5} {'lam2':>7} {'+-':>5} {'lam3':>7} {'+-':>5} "
          f"{'Rg2':>7} {'l2/l1':>6} {'l3/l1':>6} {'<E_nb>':>7} {'HH<3a':>6} {'HH<2a':>6} {'run/tau':>7} {'acc_piv':>7} {'acc_tw':>6}")
    for pat in PATTERNS:
        for pot in POTS:
            if (pat, pot) not in res:
                continue
            r = res[(pat, pot)]; lam, err = r["lam"], r["err"]
            print(f"{pat:>8} {pot:>4} {r['n']:>5} {lam[0]:>8.2f} {err[0]:>5.2f} {lam[1]:>7.2f} {err[1]:>5.2f} "
                  f"{lam[2]:>7.2f} {err[2]:>5.2f} {lam.sum():>7.2f} {lam[1] / lam[0]:>6.3f} {lam[2] / lam[0]:>6.3f} "
                  f"{r['enb']:>7.2f} {r['nc3']:>6.3f} {r['nc2']:>6.3f} {r['nsw'] / r['tau']:>7.1f} {r['acc_piv']:>7.3f} {r['acc_tw']:>6.3f}")
    print("\ndb - fit per pattern, in units of the combined error (z):")
    for pat in PATTERNS:
        if (pat, "db") in res and (pat, "fit") in res:
            a, b = res[(pat, "db")], res[(pat, "fit")]
            z = (a["lam"] - b["lam"]) / np.sqrt(a["err"] ** 2 + b["err"] ** 2)
            d = a["lam"] - b["lam"]
            print(f"{pat:>8}  dlam = {d[0]:+.2f} {d[1]:+.2f} {d[2]:+.2f}   z = {z[0]:+.1f} {z[1]:+.1f} {z[2]:+.1f}")
    if ("allL", "db") in res and ("allR", "db") in res:
        a, b = res[("allL", "db")], res[("allR", "db")]
        z = (a["lam"] - b["lam"]) / np.sqrt(a["err"] ** 2 + b["err"] ** 2)
        print(f"\nmirror check allL vs allR (db): z = {z[0]:+.1f} {z[1]:+.1f} {z[2]:+.1f}")


def figure(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4), sharex=True)
    x = np.arange(len(PATTERNS))
    for k, pot in enumerate(POTS):
        for i in range(3):
            y = [res[(p, pot)]["lam"][i] if (p, pot) in res else np.nan for p in PATTERNS]
            e = [res[(p, pot)]["err"][i] if (p, pot) in res else np.nan for p in PATTERNS]
            ax[i].errorbar(x + (k - 0.5) * 0.18, y, e, fmt="o" if pot == "db" else "s", ms=6, capsize=3,
                           color="C3" if pot == "db" else "C0", label={"db": "table (db)", "fit": "radial fit"}[pot])
    for i in range(3):
        ax[i].set_ylabel(rf"$\lambda_{i + 1}$ [$a^2$]")
        ax[i].set_xticks(x); ax[i].set_xticklabels(["all L", "all R", "L+C", "L+R"])
    ax[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig("gyration_db.pdf"); fig.savefig("gyration_db.png", dpi=150)
    print("wrote gyration_db.pdf/png")


if __name__ == "__main__":
    res = collect()
    if not res:
        print("no finished runs"); sys.exit(1)
    report(res)
    if "--fig" in sys.argv:
        figure(res)
