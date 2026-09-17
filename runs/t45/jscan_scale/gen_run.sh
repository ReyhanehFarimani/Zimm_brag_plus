#!/bin/bash
# t45: amplify the WHOLE measured helix-helix potential, achiral and chiral sector by the same factor S, per user 2026-09-17:
#   U = S * ( A Ia + C Ic )        ->  hf_scale = S, hf_eps_s = 1  (as measured shape, default range sigma_c = 0.51 a)
# to be compared at the same J with the chiral-only amplification U = A Ia + S * C Ic (hf_scale = 1, hf_eps_s = S) of the
# eps_s scans, and with eps_s = 0.  S in {2 2.9} are values that exist in the eps_s scan at this range.  Because both
# sectors grow together the potential keeps its shape and stays purely repulsive for every S.
#   J in {3 4.5 5 5.5}: the helix-coil midpoint and the handedness crossover at N = 200;  3 seeds.
# Bonds/bends from ../sample_data.dat; all lengths in a (hf_len = 1).
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
JS=${JS:-"3 4.5 5 5.5"}; SCALES=${SCALES:-"2 2.9"}; SEEDS=${SEEDS:-"1 2 3"}; NEQ=${NEQ:-100000}; NSW=${NSW:-1000000}
mkdir -p inputs out logs; n=0
for J in $JS; do for S in $SCALES; do for s in $SEEDS; do
n=$((n+1)); b=J${J}_scale${S}_s${s}; [ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" -e "s/^J1 .*/J1          = 0/" -e "s/^J2 .*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = fit/' -e "s/^n_equil.*/n_equil     = $NEQ/" -e "s/^n_sweeps.*/n_sweeps    = $NSW/" -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((4700000 + n))/" -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
printf "E_helix = 0\nhf_theta0 = 45\nhf_len = 1.0\nhf_clamp = 1\nhf_eps_s = 1\nhf_scale = $S\nn_flip = 4\n" >> inputs/$b.dat
done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*.dat | LC_ALL=C sort -t_ -k3,3 -k1,1r | xargs -P ${NPROC:-20} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
