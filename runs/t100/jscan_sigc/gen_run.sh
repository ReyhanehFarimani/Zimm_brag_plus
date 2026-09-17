#!/bin/bash
# t100 chiral-RANGE scan, per user 2026-09-17: the jscan_eps scan repeated at three chiral widths above the
# campaign median sigma_c = 0.70 a:   hf_sig_c in {0.78, 0.85, 0.93} a.
#   VALID range: the t100 chiral width rests on 5 families (median 0.70, spread 0.26 a), and 0.93 a is the achiral
#   width sigma_a -- above it the chiral force outlives the repulsion at large r and the eps_s limit collapses
#   (first attraction at eps_s = 6.1 for 1.10 a).  Inside [0.70, 0.93] the limit RISES (8.2 -> 10.5), so
#   eps_s <= 8 stays purely repulsive at all three widths.
#   J in {1..6} (E_CH = 0, E_HH = -J, E_RL = +J, E_helix = 0),  eps_s in {1, 3, 4, 7, 8},  3 seeds   (270 runs)
#   eps_s = 0 is independent of the width: use the jscan_eps runs.
# Bonds/bends verbatim from ../sample_data.dat (theta0 carries "deg").  Needs the binary with hf_sig_c.
# Dispatch: replica-first, the largest eps_s first (largest possible effect), high J first within that.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
JS=${JS:-"1 2 3 4 5 6"}
ESS=${ESS:-"1 3 4 7 8"}
SCS=${SCS:-"0.78 0.85 0.93"}
SEEDS=${SEEDS:-"1 2 3"}
NEQ=${NEQ:-100000}
NSW=${NSW:-1000000}
mkdir -p inputs out logs
k=0
for SC in $SCS; do k=$((k+1)); for J in $JS; do for ES in $ESS; do for s in $SEEDS; do
b=J${J}_es${ES}_sc${SC}_s${s}
[ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^J1 .*/J1          = 0/" \
    -e "s/^J2 .*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = fit/' \
    -e "s/^n_equil.*/n_equil     = $NEQ/" \
    -e "s/^n_sweeps.*/n_sweeps    = $NSW/" \
    -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((800000 + 100000*k + 10000*J + 100*ES + s))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = 0
hf_theta0 = 100
hf_eps_s = $ES
hf_len = 1.2311
hf_clamp = 1
hf_sig_c = $SC
n_flip = 4
EOT
done; done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*.dat | sort -t_ -k4,4 -k2,2r -k1,1r | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
