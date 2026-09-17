#!/bin/bash
# (E_helix, J2, eps_s) grid + chain-length ladder of the measured
# (nb_hh = fit) helix-helix potential: handedness-selection transition
# search.  Long runs (1M sweeps) per user 2026-09-17.
#   grid  : N = 200, E in {-8, -10, -12}, J2 in {-15, -5, +5},
#           eps_s in {0, 1, 2, 4, 6}, 2 seeds                  (90 runs)
#   ladder: N in {100, 400}, E = -12, J2 in {-15, +5},
#           eps_s in {0, 1, 2, 4, 6}, 2 seeds                  (40 runs)
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
mkdir -p inputs out logs
gen () { # N E J ES s tag
local N=$1 E=$2 J=$3 ES=$4 s=$5 b=$6 n=$7
sed -e "s/^N  .*/N           = $N/" \
    -e "s/^J2.*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = fit/' \
    -e 's/^n_equil.*/n_equil     = 100000/' \
    -e 's/^n_sweeps.*/n_sweeps    = 1000000/' \
    -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((500000 + 1000*n + 7*s))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = $E
hf_theta0 = 100
hf_eps_s = $ES
hf_len = 1.2311
n_flip = 4
EOT
}
n=0
for E in -8 -10 -12; do for J in -15 -5 5; do for ES in 0 1 2 4 6; do for s in 1 2; do
n=$((n+1)); gen 200 $E $J $ES $s E${E}_J${J}_es${ES}_s${s} $n
done; done; done; done
for N in 100 400; do for J in -15 5; do for ES in 0 1 2 4 6; do for s in 1 2; do
n=$((n+1)); gen $N -12 $J $ES $s N${N}_J${J}_es${ES}_s${s} $n
done; done; done; done
ls inputs/E*.dat inputs/N*.dat | xargs -P ${NPROC:-23} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log 2>/dev/null | grep -cv smoke) / $(ls inputs/E*.dat inputs/N*.dat | wc -l)"
