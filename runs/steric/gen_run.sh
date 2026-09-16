#!/bin/bash
# Same J-scan as runs/no_steric (user's bonds and bends, theta0 = 180 deg, N = 200, J0 = J2 = J, J1 = 0, 3 seeds)
# with the steric interaction switched on: Gaussian core (A = 9.01, sigma = 1.5) for coil-coil and coil-helix,
# repulsive Gay-Berne rods (L = 2.5, D = 2) for helix-helix. Sampling: local moves + 4 pivots per sweep
# (the hinge move is not valid with sterics).
cd "$(dirname "$0")"
BIN=${BIN:-../../zimm_dev}
THETA0=${THETA0:-180}
mkdir -p inputs out logs
for J in 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0; do
for s in 1 2 3; do
sed -e "s/^J0 .*/J0          = $J/" -e "s/^J2 .*/J2          = $J/" \
    -e "s/^\(theta0_[A-Z]*\) .*/\1 = $THETA0 deg/" \
    -e "s/^init .*/init        = walk/" -e "s/^max_disp .*/max_disp    = 1.0/" \
    -e "s/^n_equil .*/n_equil     = 30000/" -e "s/^n_sweeps .*/n_sweeps    = 200000/" \
    -e "s/^log_every .*/log_every   = 100/" -e "s/^dump_every .*/dump_every  = 1000/" \
    -e "s/^seed .*/seed        = $((5000 + 100*s))${J/./}/" \
    -e "s#^out_prefix .*#out_prefix  = out/J${J}_s${s}#" data_sample.dat > inputs/J${J}_s${s}.dat
cat >> inputs/J${J}_s${s}.dat <<INP
nb_type = gauss
nb_A = 9.01
nb_sigma = 1.5
nb_hh = gb
gb_eps0 = 1.0
rod_L = 2.5
rod_r = 1.0
n_pivot = 4
INP
done; done
ls inputs/*.dat | sort -t J -k2 -n -r | xargs -P ${NPROC:-5} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
