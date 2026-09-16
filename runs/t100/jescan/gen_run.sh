#!/bin/bash
# 2D scan of the t100 model (N = 200, sterics on): J0 = -6 .. -12 (E_HH = -J0 = +6 .. +12) x E_helix = -9 .. -12.
cd "$(dirname "$0")"
mkdir -p inputs out logs
for J in -6 -7 -8 -9 -10 -11 -12; do for E in -9 -10 -11 -12; do for s in 1 2; do
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^seed .*/seed        = $((9000 + 100*s))${J//-/}${E//-/}/" \
    -e "s/^n_sweeps .*/n_sweeps    = 200000/" -e "s/^dump_every .*/dump_every  = 0/" \
    -e "s#^out_prefix .*#out_prefix  = out/J${J}_E${E}_s${s}#" ../sample_data.dat > inputs/J${J}_E${E}_s${s}.dat
echo "E_helix = $E" >> inputs/J${J}_E${E}_s${s}.dat
done; done; done
ls inputs/*.dat | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; ../../../zimm {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
