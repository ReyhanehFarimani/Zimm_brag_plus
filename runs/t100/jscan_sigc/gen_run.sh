#!/bin/bash
# UNITS: bonds, coil core and the fitted helix-helix potential are ALL in bead diameters a, so hf_len = 1.
# Until 2026-09-17 these scans used hf_len = 1.2311 (= rod_L / L_block), which stretched the helix-helix potential by
# 23 % relative to the chain; that data is archived in prev_hflen1.2311/ and superseded.
# t100 chiral-RANGE scan, per user 2026-09-17: two chiral widths above the campaign median sigma_c = 0.70 a, each with
# FIVE eps_s values evenly spaced up to that width's own limit (the eps_s at which a radially attractive
# configuration first appears; 40M-configuration scan, clamped potential, all orientations):
#   hf_sig_c = 0.80 a : first attraction at eps_s =  9.2  ->  eps_s in {1.75, 3.5, 5.25, 7, 8.75}
#   hf_sig_c = 0.93 a : first attraction at eps_s = 10.4  ->  eps_s in {2, 4, 6, 8, 10}
#   (0.93 a = the achiral width sigma_a; above it the limit collapses, 6.1 at 1.10 a.  The sampled limits are upper
#    estimates -- they dropped from 9.4 / 10.5 at 6M configurations -- hence the ~5% margin.)
#   J in {1..6} (E_CH = 0, E_HH = -J, E_RL = +J, E_helix = 0),  3 seeds                        (180 runs)
#   eps_s = 0 is independent of the width: use the jscan_eps runs.
# Bonds/bends verbatim from ../sample_data.dat (theta0 carries "deg").  Needs the binary with hf_sig_c.
# Dispatch: replica-first, the largest eps_s first (largest possible effect), high J first within that.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
JS=${JS:-"1 2 3 4 5 6"}
SEEDS=${SEEDS:-"1 2 3"}
NEQ=${NEQ:-100000}
NSW=${NSW:-1000000}
ess_of () { case $1 in 0.80) echo "1.75 3.5 5.25 7 8.75";; 0.93) echo "2 4 6 8 10";; esac; }
mkdir -p inputs out logs
n=0
for SC in 0.80 0.93; do for J in $JS; do for ES in $(ess_of $SC); do for s in 1 2 3 4 5 6 7 8 9; do
n=$((n+1))                                   # the seed depends on (SC, J, ES, s) only, not on SEEDS
case " $SEEDS " in *" $s "*) ;; *) continue;; esac
b=J${J}_es${ES}_sc${SC}_s${s}
[ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^J1 .*/J1          = 0/" \
    -e "s/^J2 .*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = fit/' \
    -e "s/^n_equil.*/n_equil     = $NEQ/" \
    -e "s/^n_sweeps.*/n_sweeps    = $NSW/" \
    -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((900000 + n))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = 0
hf_theta0 = 100
hf_eps_s = $ES
hf_len = 1.0
hf_clamp = 1
hf_sig_c = $SC
n_flip = 4
EOT
done; done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
# LC_ALL=C: in the de_DE locale of this machine "." is a thousands separator and sort -n reads 8.75 as 875
ls inputs/J*.dat | LC_ALL=C sort -t_ -k4,4 -k2.3,2nr -k1,1r | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
