#!/bin/bash
# t100 parameter set (sample_data.dat) scanned over the helix ground energy E_helix, 3 seeds each.
# SIGN: in the code +E_helix is added per R/L residue, so a helix GROUND energy of 8..12 below coil is E_helix = -8..-12.
cd "$(dirname "$0")"
BIN=${BIN:-../../zimm}
EH=${EH:-"-8 -10 -10.5 -10.75 -11 -11.25 -11.5 -12"}
mkdir -p inputs out logs
for E in $EH; do for s in 1 2 3; do
sed -e "s/^seed .*/seed        = $((7000 + 100*s))${E//[-.]/}/" -e "s#^out_prefix .*#out_prefix  = out/E${E}_s${s}#" sample_data.dat > inputs/E${E}_s${s}.dat
echo "E_helix = $E" >> inputs/E${E}_s${s}.dat
done; done
ls inputs/*.dat | xargs -P ${NPROC:-5} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
