#!/bin/bash
# t100 single-coupling scan with the measured (nb_hh = fit) helix-helix potential, per user 2026-09-17.
# State energies reduced to ONE coupling J:
#   E_CC = 0,  E_CH = 0 (J1 = 0),  E_HH = -J (J0 = J),  E_RL = +J (J2 = J),  E_helix = 0
#   J     in {1, 2, 3, 4, 5, 6}
#   eps_s in {0, 1, 3, 4, 7}   chiral AMPLITUDE only; the chiral range (sigma_c, hard-coded) is not changed
#   3 seeds                                                                              (90 runs)
# Bonds and bends are taken verbatim from ../sample_data.dat.  UNITS: the parser reads theta0_* in
# RADIANS unless the value carries a "deg" suffix; sample_data.dat writes "180 deg" -- keep the suffix.
# hf_clamp = 1: the chiral sector acts inside the fitted domain only (|e| <= 0.5, r in [1.5, 2.5] a); without it
# eps_s >= 3 opens extrapolation wells of -33 .. -136 kT.  Needs the binary with the relabelling fix (psi is
# invariant under 1 <-> 2) -- before it E(1,2) != E(2,1) for eps_s != 0.
# Cost (measured 2026-09-17, N = 200, neighbour list on): 0.2 - 0.4 ms/sweep, i.e. ~5 min per run.
#   GEN_ONLY=1 : write the inputs and stop.   SEEDS="4 5 6" : add replicas later (existing runs are skipped).
#   Workers claim a run with an atomic mkdir lock, so the script can be started again with a larger
#   NPROC to add workers to a live batch.  After a kill, remove logs/*.lock of unfinished runs to resume.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
JS=${JS:-"1 2 3 4 5 6"}
ESS=${ESS:-"0 1 3 4 7"}
SEEDS=${SEEDS:-"1 2 3"}
NEQ=${NEQ:-100000}
NSW=${NSW:-1000000}
mkdir -p inputs out logs
for J in $JS; do for ES in $ESS; do for s in $SEEDS; do
b=J${J}_es${ES}_s${s}
[ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^J1 .*/J1          = 0/" \
    -e "s/^J2 .*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = fit/' \
    -e "s/^n_equil.*/n_equil     = $NEQ/" \
    -e "s/^n_sweeps.*/n_sweeps    = $NSW/" \
    -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((700000 + 10000*J + 100*ES + s))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = 0
hf_theta0 = 100
hf_eps_s = $ES
hf_len = 1.2311
hf_clamp = 1
n_flip = 4
EOT
done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
# dispatch replica-first (a complete (J, eps_s) grid per seed), the expensive high-J runs first within a seed
ls inputs/J*.dat | sort -t_ -k3,3 -k1,1r | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
