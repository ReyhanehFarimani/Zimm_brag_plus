#!/bin/bash
# J-scan of the no-steric system (data_sample.dat): N = 200, coil allowed (J1 = 0), soft state-dependent bonds,
# state-dependent bending (kappa per triple class) with a common preferred angle THETA0 (degrees, valence angle:
# 180 = straight). J0 = J2 = J is varied; 3 seeds per J. Everything else from data_sample.dat.
# Sampling: walk start, max_disp = 1.5, plus 20 hinge moves per sweep (valid: no non-bonded interactions) so that
# Re^2 of the stiff high-J chains decorrelates in ~50 sweeps instead of ~1e6.
cd "$(dirname "$0")"
THETA0=${THETA0:-180}          # degrees, valence angle: 180 = straight
mkdir -p inputs out logs
for J in 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0; do
for s in 1 2 3; do
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^J2 .*/J2          = $J/" \
    -e "s/^\(theta0_[A-Z]*\) .*/\1 = $THETA0 deg/" \
    -e "s/^init .*/init        = walk/" \
    -e "s/^max_disp .*/max_disp    = 1.5\nn_hinge     = 20/" \
    -e "s/^n_sweeps .*/n_sweeps    = 3000000/" \
    -e "s/^n_equil .*/n_equil     = 200000/" \
    -e "s/^seed .*/seed        = $((3000 + 100*s))${J/./}/" \
    -e "s#^out_prefix .*#out_prefix  = out/J${J}_s${s}#" \
    data_sample.dat > inputs/J${J}_s${s}.dat
done; done
# skip runs that already finished (lets an interrupted batch be resumed); NPROC limits concurrency
ls inputs/*.dat | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; ../../zimm {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
