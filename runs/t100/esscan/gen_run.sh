#!/bin/bash
# (E_helix, J2, eps_s) grid of the measured (nb_hh = fit) helix-helix
# potential at the t100 parameter set: looking for a handedness-selection
# transition driven by the chiral force (order parameter m = (nR-nL)/nH).
#   E_helix in {-8, -10, -12}   (helix-coil axis)
#   J2      in {-15, -5, +5}    (local R/L junction coupling; <0 favors mixing)
#   eps_s   in {0, 1, 2, 4, 6}  (chiral amplification; F_r repulsive < 7.8)
# 2 seeds each; the original es*_s* runs (E -12, J2 -15) stay as extra seeds.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
mkdir -p inputs out logs
n=0
for E in -8 -10 -12; do for J in -15 -5 5; do for ES in 0 1 2 4 6; do for s in 1 2; do
n=$((n + 1))
b=E${E}_J${J}_es${ES}_s${s}
sed -e "s/^J2.*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = fit/' \
    -e 's/^dump_every.*/dump_every  = 10000/' \
    -e "s/^seed .*/seed        = $((300000 + 1000*n + 7*s))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = $E
hf_theta0 = 100
hf_eps_s = $ES
hf_len = 1.2311
n_flip = 4
EOT
done; done; done; done
ls inputs/E*.dat inputs/es*.dat | xargs -P ${NPROC:-23} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log 2>/dev/null | grep -cv smoke) / $(ls inputs/E*.dat inputs/es*.dat | wc -l)"
