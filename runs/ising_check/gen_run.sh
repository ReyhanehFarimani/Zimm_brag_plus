#!/bin/bash
# 1D Ising check of the state part. All bond parameters equal, so the spins decouple from geometry
# and the state statistics are an exact 3-state (spin-1) transfer-matrix problem.
#   scan A: kT scan at J0=1, J1=1, J2=1.5  -> helicity and <E_state> vs kT
#   scan B: correlation function <s_i s_{i+r}> at kT=1, J0=1, J1=2.5, J2=1 (coil nearly suppressed:
#           2-state Ising limit with J_eff = (J0+J2)/2), 5 independent seeds
cd "$(dirname "$0")"
BIN=${BIN:-../../zimm_dev}
N=${N:-1024}
mkdir -p inputs out logs
mk() {  # name kT J0 J1 J2 dump_every [seed] [n_equil]
cat > inputs/$1.dat <<INP
N = $N
init = rod
seed = ${7:-4242}
bond_len_CC = 1.0
k_bond_CC = 100.0
bond_len_CH = 1.0
k_bond_CH = 100.0
bond_len_HH = 1.0
k_bond_HH = 100.0
bond_len_RL = 1.0
k_bond_RL = 100.0
J0 = $3
J1 = $4
J2 = $5
kT = $2
n_equil = ${8:-50000}
n_sweeps = 1000000
max_disp = 0.25
log_every = 50
dump_every = $6
out_prefix = out/$1
INP
}
for kT in 0.3 0.4 0.5 0.6 0.8 1.0 1.25 1.5 2.0 3.0; do mk A_kT$kT $kT 1.0 1.0 1.5 0; done
for s in 11 12 13 14 15; do mk B_s$s 1.0 1.0 2.5 1.0 1000 $s 200000; done
ls inputs/*.dat | xargs -P ${NPROC:-5} -I{} sh -c 'b=$(basename {} .dat); '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
