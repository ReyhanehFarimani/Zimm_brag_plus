#!/bin/bash
# Ideal-chain (no excluded volume) test against Chauhan, Mueller, Daoulas, Macromolecules 58, 5408 (2025):
# freely-jointed chain with reversible hinges. Mapping onto our model:
#   * bond i carries the spin of bead i (bend_key = pair): the hinge at bead j is keyed by (state[j-1], state[j])
#   * HH hinge: very stiff harmonic bend (kappa_HHH, theta0 = pi = straight)  ~ rigid;  RL hinges: free (kappa = 0)
#   * sampled with the spin-flip + hinge-resampling move (n_hinge), whose acceptance carries only the free-energy
#     cost of a wall, so kappa can be made as stiff as we like
#   * E_HH = -J0, E_RL = +J0 (J2 = J0); two-state model (n_states = 2): no coil, R <-> L flipped directly
#     (a coil intermediate is entropically favoured once the helix hinges are stiff, and R->C->L would be far too slow)
#   * effective Ising coupling J = J0 + 1/2 ln(w_r/w_f), w_r/w_f = 1/2 int_0^pi sin(t) exp(-kappa (t-pi)^2/2) dt
#     (paper: J = J0 - ln(4pi)/2 for perfectly rigid hinges)
#   * from the rod start the wall count and Re^2 need ~5e6 sweeps to relax (walls annihilate while all hinges
#     are still straight, and are re-created only at rate ~exp(-4 J0)); hence the long n_equil
cd "$(dirname "$0")"
BIN=${BIN:-../../zimm}
N=${N:-101}          # M = N-1 = 100 bonds, as in the paper's figures
mkdir -p inputs out logs
J0_for() {  # J kappa -> J0
python3 - "$1" "$2" <<'PY'
import sys, numpy as np
J, k = float(sys.argv[1]), float(sys.argv[2])
t = np.linspace(0, np.pi, 400001); w = np.sin(t) * np.exp(-0.5 * k * (t - np.pi)**2)
ratio = 0.5 * np.trapezoid(w, t)
print(f"{J - 0.5*np.log(ratio):.6f}")
PY
}
mk() {  # tag J kappa seed
J0=$(J0_for $2 $3)
cat > inputs/$1.dat <<INP
# target J = $2, kappa = $3  ->  J0 = $J0
N = $N
init = rod
init_state = random
n_states = 2
bend_key = pair
seed = $4
bond_len_CC = 1.0
k_bond_CC = 1000.0
bond_len_CH = 1.0
k_bond_CH = 1000.0
bond_len_HH = 1.0
k_bond_HH = 1000.0
bond_len_RL = 1.0
k_bond_RL = 1000.0
kappa_HHH = $3
theta0_HHH = 3.141592653589793
kappa_CCC = 0.0
kappa_CHC = 0.0
kappa_RLR = 0.0
J0 = $J0
J1 = 2.0
J2 = $J0
kT = 1.0
n_equil = 500000
n_sweeps = 5000000
max_disp = 0.06
n_hinge = 10
log_every = 500
dump_every = 500
out_prefix = out/$1
INP
}
for k in 100 400 1600 10000; do for s in 1 2 3 4 5; do mk J0.5_k${k}_s$s 0.5 $k $((100*k + s)); done; done
ls inputs/*.dat | xargs -P ${NPROC:-10} -I{} sh -c 'b=$(basename {} .dat); '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
