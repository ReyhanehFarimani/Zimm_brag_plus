#!/bin/bash
# Finite-size check of the t100 crossover with sterics: N = 50 and 100 (N = 200 from ../out), E_helix near E*.
cd "$(dirname "$0")"
mkdir -p inputs out logs
for N in 50 100; do for E in -10.4 -10.5 -10.6 -10.7 -10.8 -11.0; do for s in 1 2; do
sed -e "s/^N .*/N           = $N/" -e "s/^seed .*/seed        = $((8000 + 100*s + N))${E//[-.]/}/" -e "s#^out_prefix .*#out_prefix  = out/N${N}_E${E}_s${s}#" \
    -e "s/^n_sweeps .*/n_sweeps    = 300000/" -e "s/^dump_every .*/dump_every  = 0/" ../sample_data.dat > inputs/N${N}_E${E}_s${s}.dat
echo "E_helix = $E" >> inputs/N${N}_E${E}_s${s}.dat
done; done; done
ls inputs/*.dat | sort -t N -k2 -n -r | xargs -P ${NPROC:-8} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; ../../../zimm {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
