#!/bin/bash
# Reproduce Fig. 2 of Chauhan/Mueller/Daoulas (Macromolecules 2025): <Re^2>_id / <Re^2>_0 vs number of bonds M
# for J = 0.5 and J = 2, M = 100 ... 1500 (their sizes), ideal chain (no excluded volume).
# Same mapping as runs/paper_ideal (bond spin = left bead, bend_key = pair, n_states = 2), rigid limit kappa = 1e4,
# sampled with n_hinge = N/10 hinge moves per sweep.
cd "$(dirname "$0")"
BIN=${BIN:-../../zimm}
mkdir -p inputs out logs
J0_for() {  # J kappa -> J0
python3 - "$1" "$2" <<'PY'
import sys, numpy as np
J, k = float(sys.argv[1]), float(sys.argv[2])
t = np.linspace(0, np.pi, 400001); w = np.sin(t) * np.exp(-0.5 * k * (t - np.pi)**2)
print(f"{J - 0.5*np.log(0.5*np.trapezoid(w, t)):.6f}")
PY
}
KAPPA=${KAPPA:-10000}
mk() {  # M J seed nsweeps
N=$(( $1 + 1 )); J0=$(J0_for $2 $KAPPA)
cat > inputs/M$1_J$2_s$3.dat <<INP
# M = $1 bonds, target J = $2, kappa = $KAPPA -> J0 = $J0
N = $N
init = rod
init_state = random
n_states = 2
bend_key = pair
seed = $(( 7919*$1 + 100*$3 + ${2/./} ))
bond_len_CC = 1.0
k_bond_CC = 1000.0
bond_len_CH = 1.0
k_bond_CH = 1000.0
bond_len_HH = 1.0
k_bond_HH = 1000.0
bond_len_RL = 1.0
k_bond_RL = 1000.0
kappa_HHH = $KAPPA
theta0_HHH = 3.141592653589793
kappa_CCC = 0
kappa_CHC = 0
kappa_RLR = 0
J0 = $J0
J1 = 2.0
J2 = $J0
kT = 1.0
n_equil = $(( $4 / 4 ))
n_sweeps = $4
max_disp = 0.06
n_hinge = $(( N / 10 ))
log_every = 20
dump_every = 1000
out_prefix = out/M$1_J$2_s$3
INP
}
for M in 100 200 300 400 500 600 800 1000 1200 1500; do
  for s in 1 2 3 4 5 6 7 8 9 10; do
    mk $M 0.5 $s 200000
    mk $M 2.0 $s 400000     # J = 2: few, long-lived walls and large Re^2 fluctuations -> longer run
  done
done
# largest first
ls inputs/*.dat | sort -t M -k2 -n -r | xargs -P ${NPROC:-10} -I{} sh -c 'b=$(basename {} .dat); '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
