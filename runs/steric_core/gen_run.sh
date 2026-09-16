#!/bin/bash
# Steric-core test: a coil of x soft beads (Gaussian core A e^{-(r/sigma)^2}, A = 9.01, sigma = 1.5, coil bonds
# from the user's parameter set) vs a chain of 7x WCA beads (sigma = 1, eps = 1, harmonic bonds r0 = 1, k = 100).
# One soft bead is meant to represent 7 fine monomers; if the mapping holds the gyration-tensor eigenvalues of
# the two chains agree up to one common length-scale factor. Non-bonded pairs only (|i-j| > 1), no bending,
# states frozen in coil (J1 = 50). Sampled with local moves + pivots.
cd "$(dirname "$0")"
mkdir -p inputs out logs
mk() {  # tag N nb_type A sigma bond_len k_bond max_disp n_sweeps seed
cat > inputs/$1.dat <<INP
N = $2
init = walk
seed = ${10}
J0 = 0
J1 = 50
J2 = 0
nb_type = $3
nb_A = $4
nb_sigma = $5
kT = 1.0
n_equil = $(( $9 / 5 ))
n_sweeps = $9
max_disp = $8
n_pivot = 4
max_rot = 3.14159265358979
log_every = $(( $9 / 2000 ))
dump_every = $(( $9 / 1000 ))
out_prefix = out/$1
INP
for c in CC CH HH RL; do printf "bond_len_$c = $6\nk_bond_$c = $7\n" >> inputs/$1.dat; done
for t in CCC CCH CHC CHH CRL HCH RCL HHH RRL RLR; do printf "kappa_$t = 0\n" >> inputs/$1.dat; done
}
for x in 8 16 32 64; do
  for s in 1 2 3; do
    mk gauss_x${x}_s$s $x   gauss 9.01 1.5 2.46 0.987 1.5 200000 $((500 + 10*x + s))
    mk wca_x${x}_s$s $((7*x)) wca 1.0 1.0 1.0 100.0 0.2 $((x <= 16 ? 200000 : 100000)) $((900 + 10*x + s))
  done
done
# biggest WCA chains first; skip finished runs (resumable); NPROC limits concurrency (user wants 4 cores free)
ls inputs/*.dat | sort -t x -k2 -n -r | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; ../../zimm {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
