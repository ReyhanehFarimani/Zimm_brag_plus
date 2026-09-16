#!/bin/bash
# Boundary zoom: fine J0 steps (0.25) across the helix-coil boundary for each E_helix, 3 seeds, N = 200, sterics on.
cd "$(dirname "$0")"
mkdir -p inputs out logs
pairs() { cat <<'P'
-9  -8.0 -7.75 -7.5 -7.25 -7.0
-10 -9.0 -8.75 -8.5 -8.25 -8.0
-11 -11.0 -10.75 -10.5 -10.25 -10.0
-12 -13.0 -12.75 -12.5 -12.25 -12.0
P
}
pairs | while read E J0s; do for J in $J0s; do for s in 1 2 3; do
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^seed .*/seed        = $((9500 + 100*s))${J//[-.]/}${E//-/}/" \
    -e "s/^n_sweeps .*/n_sweeps    = 200000/" -e "s/^dump_every .*/dump_every  = 0/" \
    -e "s#^out_prefix .*#out_prefix  = out/J${J}_E${E}_s${s}#" ../sample_data.dat > inputs/J${J}_E${E}_s${s}.dat
echo "E_helix = $E" >> inputs/J${J}_E${E}_s${s}.dat
done; done; done
ls inputs/*.dat | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; ../../../zimm {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
