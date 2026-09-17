#!/bin/bash
# 2026-09-17 17:50, user decision "C": the exploratory range (t45 0.87 a / t100 0.93 a, not supported by the family fits) is dropped for
# all runs not yet started, to shorten the queue; unstarted inputs moved to inputs_dropped/. SCS="..." restores it.
# theta0 = 45 scan, per user 2026-09-17: THREE chiral ranges, FIVE eps_s each up to that range's own limit, same J.
# The t45 potential has a ~3x stronger measured chiral sector than t100 (max 0.75 vs 0.24 kT per unit eps_s), so the
# limits (eps_s at which a radially attractive configuration first appears; clamped potential, all orientations,
# 12M configurations) are much lower:
#   hf_sig_c = 0.51 a  campaign family median (default)      first attraction at 3.1  ->  eps_s in {1, 1.5, 2, 2.5, 2.9}
#   hf_sig_c = 0.62 a  median + family spread (0.11 a)        first attraction at 4.5  ->  eps_s in {1, 1.8, 2.6, 3.4, 4.2}
#   hf_sig_c = 0.87 a  = the achiral width sigma_a            first attraction at 6.6  ->  eps_s in {1, 2.3, 3.6, 4.9, 6.2}
#   Every list starts at eps_s = 1 (AS MEASURED) and is evenly spaced up to ~6 % below the limit (sampled limits are
#   upper estimates).  0.87 a is 3 family spreads above the median: exploratory, not supported by the fits.  Below
#   the median the limit collapses (1.4 at 0.40 a, 0 at the weighted mean 0.23 a), so no smaller range is scanned.
#   + the eps_s = 0 reference, which does not depend on the range.
#   J in {1..6}: E_CC = 0, E_CH = 0 (J1 = 0), E_HH = -J (J0 = J), E_RL = +J (J2 = J), E_helix = 0;  3 seeds  (288 runs)
# Bonds, bends, coil core: ../sample_data.dat, GENERATED from the measured t45 fits by ../make_sample_data.py.
# UNITS: all lengths in a, so hf_len = 1 (NOT the 1.2311 of the t100 scans -- see the note in make_sample_data.py).
# Extra couplings (2026-09-17, to resolve the handedness crossover, exact J* = 4.68 at N = 200):
#   JS="4.5 5.5" SEED_BASE=4600000 ./gen_run.sh      (own seed block; J may be non-integer)
# Dispatch: replica-first, the largest eps_s first, high J first within that.  Lock-based workers: start the script
# again with another NPROC to add workers to a live batch.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
JS=${JS:-"1 2 3 4 5 6"}
SEEDS=${SEEDS:-"1 2 3"}
NEQ=${NEQ:-100000}
NSW=${NSW:-1000000}
ess_of () { case $1 in 0.51) echo "0 1 1.5 2 2.5 2.9";; 0.62) echo "1 1.8 2.6 3.4 4.2";; 0.87) echo "1 2.3 3.6 4.9 6.2";; esac; }
mkdir -p inputs out logs
n=0
for SC in ${SCS:-0.51 0.62}; do for J in $JS; do for ES in $(ess_of $SC); do for s in 1 2 3 4 5 6 7 8 9; do
n=$((n+1))                                   # the seed depends on (SC, J, ES, s) only, not on SEEDS
case " $SEEDS " in *" $s "*) ;; *) continue;; esac
b=J${J}_es${ES}_sc${SC}_s${s}; [ "$ES" = 0 ] && b=J${J}_es0_sc0_s${s}      # eps_s = 0: the range is irrelevant
[ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" \
    -e "s/^J1 .*/J1          = 0/" \
    -e "s/^J2 .*/J2          = $J/" \
    -e "s/^n_equil.*/n_equil     = $NEQ/" \
    -e "s/^n_sweeps.*/n_sweeps    = $NSW/" \
    -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((${SEED_BASE:-4500000} + n))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
E_helix = 0
hf_eps_s = $ES
hf_sig_c = $SC
EOT
done; done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
# LC_ALL=C: in the de_DE locale of this machine sort -n reads 2.9 as 29
ls inputs/J*.dat | LC_ALL=C sort -t_ -k4,4 -k2.3,2nr -k1,1r | xargs -P ${NPROC:-15} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
