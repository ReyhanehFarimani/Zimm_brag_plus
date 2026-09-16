#!/bin/bash
# Chirality zoom: E_helix = -11, J0 = -9 .. -6 (sense coupling (11-|J0|)/2 = 1 .. 2.5 kT), fully helical phase.
# Uses the validated sense-flip moves (n_flip: single-site R<->L + non-merging domain flips) so the handedness
# sector equilibrates; sterics on; 3 seeds.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm_dev}
mkdir -p inputs out logs
for J in -9 -8.5 -8 -7.5 -7 -6.5 -6; do for s in 1 2 3; do
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^seed .*/seed        = $((9800 + 100*s))${J//[-.]/}/" \
    -e "s/^n_sweeps .*/n_sweeps    = 300000/" -e "s/^dump_every .*/dump_every  = 1000/" \
    -e "s#^out_prefix .*#out_prefix  = out/J${J}_s${s}#" ../sample_data.dat > inputs/J${J}_s${s}.dat
printf "E_helix = -11\nn_flip = 20\n" >> inputs/J${J}_s${s}.dat
done; done
ls inputs/*.dat | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
