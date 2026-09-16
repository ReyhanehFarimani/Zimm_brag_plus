#!/usr/bin/env python3
"""Chirality zoom (E_helix = -11, J0 = -9..-6, sense-flip moves on): handedness order across the sense coupling.
<S^2>/N, sqrt(N)<|m|>, Binder U4, R|L wall count vs J0 with the exact no-steric references; P(m) histograms; m(t)."""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# exact reference (needs S2N -> use the no_steric exact_bend with J2 = 11)
sys.argv = ["x"]
ns = open("../../no_steric/analyze.py").read()
exec(ns.split("inp0 = read_input")[0])
exec(ns[ns.index("def triple_class"):ns.index("return dict(hel=hel, Es=Es, Eb=Eb, Re2=Re2, S2N=tot / N)") + len("return dict(hel=hel, Es=Es, Eb=Eb, Re2=Re2, S2N=tot / N)")])
inp0 = read_input("../sample_data.dat"); inp0["E_helix"] = "-11"; N = int(inp0["N"])

def read_obs(f):
    lines = open(f).read().splitlines(); hdr = next(l for l in lines if l.startswith("#")).split()[1:]
    rows = np.array([[float(x) for x in l.split()] for l in lines if not l.startswith("#") and len(l.split()) == len(hdr)])
    return {k: rows[:, i] for i, k in enumerate(hdr)}
def walls(f):
    out = []
    with open(f) as fh:
        while True:
            h = fh.readline()
            if not h: break
            n = int(h); fh.readline(); sp = np.array([int(fh.readline().split()[-1]) for _ in range(n)])
            out.append(np.sum(sp[:-1] * sp[1:] < 0))
    return np.mean(out) if out else np.nan
def tau_int(x):
    x = x - x.mean(); n = len(x)
    if x.std() == 0: return np.nan
    F = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(F * np.conj(F))[:n] / np.arange(n, 0, -1); ac /= ac[0]; t = 0.5
    for M in range(1, n // 2):
        t += ac[M]
        if M > 6 * t: break
    return t

cells = {}; series = {}
for f in sorted(glob.glob("out/J*_s*_obs.dat")):
    J = float(re.search(r"J(-?[\d.]+)_s", f).group(1)); o = read_obs(f)
    if len(o["helicity"]) < 40: continue
    m = (o["n_R"] - o["n_L"]) / N
    nb = 20; blk = lambda x: (x[: len(x) // nb * nb].reshape(nb, -1).mean(1).std(ddof=1) / np.sqrt(nb))
    m2 = (m**2).mean(); m4 = (m**4).mean()
    cells.setdefault(J, []).append(dict(hel=o["helicity"].mean(), S2N=(N * m2, N * blk(m**2)), am=(np.abs(m).mean(), blk(np.abs(m))),
                                        U4=1 - m4 / (3 * m2**2), tau=tau_int(m) * (o["sweep"][1] - o["sweep"][0]),
                                        w=walls(f.replace("_obs.dat", "_conf.xyz")), accf=o["acc_flip"][-1]))
    series.setdefault(J, []).append((o["sweep"], m))
Js = sorted(cells)
print(f"E_helix = -11, N = {N}, sterics + sense-flip moves; 3 seeds x 3e5 sweeps")
print(f"{'J0':>6} {'J_sense':>8} {'hel':>6} | {'S2/N':>8} {'+-':>6} {'exact-ns':>9} | {'sqrtN<|m|>':>10} {'U4':>7} {'walls':>7} {'exact-ns':>9} | {'tau(m)':>8} {'acc_flip':>8}")
exS = {}
for J in Js:
    d = dict(inp0); e = exact_bend(d, J, J2=11.0); exS[J] = e
    v = cells[J]; s2 = np.mean([x["S2N"][0] for x in v]); s2e = np.sqrt(np.sum([x["S2N"][1]**2 for x in v])) / len(v)
    am = np.mean([x["am"][0] for x in v]); u4 = np.mean([x["U4"] for x in v]); w = np.nanmean([x["w"] for x in v])
    # exact wall count: from S2N relation? use paper eq-9-like via TM is messy; estimate exact walls from the 9x9 TM pair probs
    print(f"{J:>6} {(11 - abs(J)) / 2:>8.2f} {np.mean([x['hel'] for x in v]):>6.3f} | {s2:>8.3f} {s2e:>6.3f} {e['S2N']:>9.3f} |"
          f" {np.sqrt(N) * am:>10.2f} {u4:>7.3f} {w:>7.1f} {'':>9} | {np.mean([x['tau'] for x in v]):>8.0f} {np.mean([x['accf'] for x in v]):>8.3f}")

INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.3), dpi=160, facecolor=BG)
for ax in axs:
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.set_xlabel("J0   (sense coupling = (11−|J0|)/2)", color=INK)
ax = axs[0]
Jline = np.linspace(min(Js) - 0.2, max(Js) + 0.2, 40)
ax.plot(Jline, [exact_bend(dict(inp0), Jv, J2=11.0)["S2N"] for Jv in Jline], "--", lw=1.2, color=INK2, label="exact, no sterics")
ax.errorbar(Js, [np.mean([x["S2N"][0] for x in cells[J]]) for J in Js],
            yerr=[np.sqrt(np.sum([x["S2N"][1]**2 for x in cells[J]])) / len(cells[J]) for J in Js], fmt="o-", ms=5, lw=1, color=GREEN, capsize=2.5, label="MC with sterics")
ax.set_yscale("log"); ax.set_ylabel("⟨S²⟩/N", color=INK); ax.set_title("Chiral order fluctuations", loc="left", fontsize=9.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
ax = axs[1]
ax.errorbar(Js, [np.sqrt(N) * np.mean([x["am"][0] for x in cells[J]]) for J in Js],
            yerr=[np.sqrt(N) * np.sqrt(np.sum([x["am"][1]**2 for x in cells[J]])) / len(cells[J]) for J in Js], fmt="o-", ms=5, lw=1, color=BLUE, capsize=2.5, label="√N·⟨|m|⟩")
ax.axhline(np.sqrt(2 / np.pi), color=INK2, lw=0.9, ls=":"); ax.text(Js[0], 0.85, "Gaussian if ⟨S²⟩/N=1", fontsize=7.5, color=INK2)
ax2 = ax.twinx(); ax2.plot(Js, [np.mean([x["U4"] for x in cells[J]]) for J in Js], "s--", ms=4, lw=1, color=ORANGE)
ax2.set_ylabel("U₄ (orange)", color=ORANGE); ax2.set_ylim(-0.2, 0.75); ax2.axhline(2/3, color=ORANGE, lw=0.7, ls=":")
ax.set_ylabel("√N·⟨|m|⟩", color=INK); ax.set_title("Handedness order parameter and Binder cumulant", loc="left", fontsize=9.5, color=INK)
ax = axs[2]
ax.errorbar(Js, [np.nanmean([x["w"] for x in cells[J]]) for J in Js], fmt="o-", ms=5, lw=1, color=ORANGE, capsize=2.5)
ax.set_yscale("log"); ax.set_ylabel("R|L walls per chain", color=INK); ax.set_title("Sense-reversal walls", loc="left", fontsize=9.5, color=INK)
fig.tight_layout(); fig.savefig("chirality.png")

fig, axs = plt.subplots(len(Js), 1, figsize=(11, 1.35 * len(Js) + 0.8), dpi=130, sharex=True, facecolor=BG)
for ax, J in zip(np.atleast_1d(axs), Js):
    ax.set_facecolor(BG); ax.grid(True, color=GRID, lw=0.6)
    for (s, m), col in zip(series[J], (BLUE, ORANGE, GREEN)):
        ax.plot(s / 1e5, m, lw=0.6, color=col)
    ax.axhline(0, color=INK2, lw=0.7); ax.set_ylim(-1.05, 1.05); ax.set_ylabel("m", color=INK)
    ax.text(0.005, 0.82, f"J0={J}", transform=ax.transAxes, fontsize=8.5, color=INK)
    ax.tick_params(colors=INK2, labelsize=7.5)
np.atleast_1d(axs)[-1].set_xlabel("sweeps [×10⁵]", color=INK)
np.atleast_1d(axs)[0].set_title("m(t), three seeds per panel — sense-flip moves keep the handedness sector moving", loc="left", fontsize=9.5, color=INK)
fig.tight_layout(); fig.savefig("chirality_traces.png")
print("wrote chirality.png, chirality_traces.png")
