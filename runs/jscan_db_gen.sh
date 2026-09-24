#!/bin/bash
# J scan with the TABLE interaction (user 2026-09-23: "no e_s, just the table interaction we already had"):
# nb_hh = db (six-argument helix-helix table helix_pair_db_t<theta0>_win.bin + registry twist term), one series.
#   E_HH = -J, E_CH = 0, E_RL = +J, E_helix = 0;  J in {1 .. 6, half steps};  3 seeds;  N = 200, 1M sweeps
# Bonds, bends, coil core from runs/<geom>/sample_data.dat (k = 2K, generated from the fits; nb_min_sep = 3 default).
#   usage: cd runs/<t45|t100>/jscan_db && [NPROC=12] [SEEDS="1 2 3"] [JS="..."] ../../jscan_db_gen.sh
# run from inside runs/<geom>/jscan_db (no cd: the script lives in runs/)
GEOM=$(basename "$(cd .. && pwd)"); TH=${GEOM#t}; [ -f ../sample_data.dat ] || { echo "run from runs/<geom>/jscan_db"; exit 1; }
BIN=${BIN:-../../../zimm}
DB=/home/reyhaneh/Documents/Zimm_brag_plus/helix_pair_db_${GEOM}_win.bin
JS=${JS:-"1 2 3 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 8 9 10"}   # user 2026-09-23: 3 seeds for 5.5-7.5 (the transition), 1 seed elsewhere
SEEDS=${SEEDS:-"1 2 3"}
NEQ=${NEQ:-100000}; NSW=${NSW:-1000000}
mkdir -p inputs out logs
for J in $JS; do for s in $SEEDS; do
b=J${J}_db_s${s}
[ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" -e "s/^J1 .*/J1          = 0/" -e "s/^J2 .*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = db/' -e "s/^n_equil.*/n_equil     = $NEQ/" -e "s/^n_sweeps.*/n_sweeps    = $NSW/" \
    -e 's/^dump_every.*/dump_every  = 50000/' -e "s/^seed .*/seed        = $((${TH}00000 + 1000*${J/./} + s))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = 0
db_file = $DB
hf_theta0 = $TH
hf_len = 1.0
n_twist = 200
twist_step = 30
n_flip = 4
nl_skin = 1.0
EOT
done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*.dat | LC_ALL=C sort -t_ -k3,3 -k1,1r | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
