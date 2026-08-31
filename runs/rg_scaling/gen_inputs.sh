#!/bin/bash
# Generate inputs: N = 2^5 .. 2^10, 5 seeds each; everything else identical.
# Run length scales as N^2 (Rouse scaling of the Rg autocorrelation time for local moves):
#   tau_int(Rg2) ~ 1e5 sweeps at N=256  ->  n_sweeps = 3e6 (N/256)^2 ~ 30 tau, with a 2e7 floor for small N.
cd "$(dirname "$0")"
mkdir -p inputs out logs
for N in 32 64 128 256 512 1024; do
for s in 1 2 3 4 5; do
seed=$((1000*N + s))
nsw=$(( 3000000 * (N/256) * (N/256) )); [ $nsw -lt 20000000 ] && nsw=20000000
neq=$(( nsw / 100 ))
lev=$(( nsw / 20000 ))
cat > inputs/N${N}_s${s}.dat <<INP
# Rg scaling study: bonds only (harmonic), all coil, random-walk start
N           = ${N}
init        = walk
seed        = ${seed}

bond_len_CC = 1.0
k_bond_CC   = 100.0
bond_len_CH = 1.0
k_bond_CH   = 100.0
bond_len_HH = 1.0
k_bond_HH   = 100.0
bond_len_RL = 1.0
k_bond_RL   = 100.0

kT          = 1.0
n_equil     = ${neq}
n_sweeps    = ${nsw}
max_disp    = 0.25
log_every   = ${lev}
dump_every  = 0
out_prefix  = out/N${N}_s${s}
INP
done; done
ls inputs | wc -l
