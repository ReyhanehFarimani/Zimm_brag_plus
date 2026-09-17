#!/bin/bash
# eps_s scan of the measured (nb_hh = fit) helix-helix potential at the t100
# parameter set with E_helix = -12, J2 = -15 (R/L junctions favored locally):
# does the nonlocal chiral interaction select a handedness anyway?
# eps_s in {0, 1, 2, 4, 6} (F_r stays repulsive for eps_s < 7.8), 2 seeds.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
mkdir -p inputs out logs
for ES in 0 1 2 4 6; do for s in 1 2; do
b=es${ES}_s${s}
sed -e 's/^J2.*/J2          = -15/' \
    -e 's/^nb_hh.*/nb_hh = fit/' \
    -e 's/^dump_every.*/dump_every  = 10000/' \
    -e "s/^seed .*/seed        = $((9100 + 10*s + ES))77/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = -12
hf_theta0 = 100
hf_eps_s = $ES
hf_len = 1.2311
n_flip = 4
EOT
done; done
ls inputs/es*.dat | xargs -P ${NPROC:-5} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/es*.log 2>/dev/null | wc -l) / $(ls inputs/es*.dat | wc -l)"
