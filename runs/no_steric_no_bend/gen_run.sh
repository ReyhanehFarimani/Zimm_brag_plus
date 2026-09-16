#!/bin/bash
# J-scan of the user's no-steric / no-bend system (input_1.dat): N = 200, coil allowed (J1 = 0),
# state-dependent SOFT bonds (per-class bond_len / k_bond), all bending constants zero.
# J0 = J2 = J is varied; everything else is taken from input_1.dat.
cd "$(dirname "$0")"
mkdir -p inputs out logs
for J in 0.0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5; do
for s in 1 2 3; do
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^J2 .*/J2          = $J/" \
    -e "s/^n_equil .*/n_equil     = 500000/" \
    -e "s/^n_sweeps .*/n_sweeps    = 10000000/" \
    -e "s/^log_every .*/log_every   = 1000/" \
    -e "s/^dump_every .*/dump_every  = 10000/" \
    -e "s/^init .*/init        = walk/" \
    -e "s/^max_disp .*/max_disp    = 1.5/" \
    -e "s/^seed .*/seed        = $((1000 + 100*s))${J/./}/" \
    -e "s#^out_prefix .*#out_prefix  = out/J${J}_s${s}#" \
    input_1.dat > inputs/J${J}_s${s}.dat
done; done
# skip runs that already finished (lets an interrupted batch be resumed); NPROC limits concurrency
ls inputs/*.dat | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; ../../zimm {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
