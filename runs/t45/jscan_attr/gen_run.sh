#!/bin/bash
# 2026-09-17 19:05: the user rejected the two-knob variant (hf_eps_rl) as a manipulation of the potential. DEFAULT = range 0.87 a with ONE
# eps_s for every pair type, the potential exactly as fitted. The 0.51 / 0.62 a two-knob runs are in prev_two_knob_hf_eps_rl/.
# t45: an attraction between SAME-handed helices only -- opposite-handed (R.L) pairs keep the measured, purely repulsive
# potential (per user 2026-09-17: "until -3 kBT but no attraction between R and L").
#   hf_eps_s  = chiral amplitude of the same-handed pairs, chosen for the depth of their deepest well (at r = 1.5 a, one sign of psi)
#   hf_eps_rl = 1  -> R.L pairs as measured: radial force repulsive everywhere, no local pocket (checked: min F_r = +0.11 kT/a)
#     range 0.51 a: eps_s in {8, 8.5, 9.2, 10.7, 12.2}     ->  wells of about {-0.2, -0.5, -1, -2, -3} kT
#     range 0.62 a: eps_s in {8.4, 8.85, 9.6, 11.2, 12.75} ->  wells of about {-0.2, -0.5, -1, -2, -3} kT
#     range 0.87 a: ONE eps_s (no hf_eps_rl) in {7.2, 7.6, 8.3, 9.6, 10.5}  ->  wells of about {-0.2, -0.5, -1, -2, -2.7} kT, R.L purely
#                   repulsive by itself -- the only range where a single eps_s does that (user request 2026-09-17 19:00)
#   J in {2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5};  3 seeds                                                     (360 runs)
#   ONLY="_sc0.87_" restricts a worker pool to one range.
# Reference: the eps_s = 0 runs of ../jscan_sigc (same J).  This is BEYOND the measured interaction (eps_s = 1): a model
# study of what a same-handed attraction would do.  The earlier runs that amplified R.L as well are in prev_RL_amplified/.
# Bonds/bends from ../sample_data.dat, all lengths in a (hf_len = 1).  Needs the binary with hf_eps_rl.
# Dispatch: replica-first, then the crossover couplings first (high J), then the deepest well first.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
JS=${JS:-"2 2.5 3 3.5 4 4.5 5 5.5"}; SEEDS=${SEEDS:-"1 2 3"}; NEQ=${NEQ:-100000}; NSW=${NSW:-1000000}
ess_of () { case $1 in 0.51) echo "8 8.5 9.2 10.7 12.2";; 0.62) echo "8.4 8.85 9.6 11.2 12.75";; 0.87) echo "7.2 7.6 8.3 9.6 10.5";; esac; }
mkdir -p inputs out logs; n=0
for SC in ${SCS:-0.87}; do for J in $JS; do for ES in $(ess_of $SC); do for s in 1 2 3 4 5 6 7 8 9; do
n=$((n+1)); case " $SEEDS " in *" $s "*) ;; *) continue;; esac
b=J${J}_es${ES}_sc${SC}_s${s}; [ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" -e "s/^J1 .*/J1          = 0/" -e "s/^J2 .*/J2          = $J/" \
    -e "s/^n_equil.*/n_equil     = $NEQ/" -e "s/^n_sweeps.*/n_sweeps    = $NSW/" -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((5100000 + n))/" -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
# 0.51 / 0.62 a: two knobs (hf_eps_rl = 1 keeps R.L as measured).  0.87 a: ONE eps_s for every pair type -- at that range
# the window eps_s = 7.0 .. 10.5 gives same-handed wells of -0.03 .. -2.7 kT while R.L stays purely repulsive by itself.
RL="hf_eps_rl = 1"; [ "$SC" = 0.87 ] && RL="# single eps_s for all pair types"
printf "E_helix = 0\nhf_eps_s = $ES\n$RL\nhf_sig_c = $SC\n" >> inputs/$b.dat
done; done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*.dat | grep -E "${ONLY:-.}" | LC_ALL=C sort -t_ -k4,4 -k1.9,1nr -k2.3,2nr | xargs -P ${NPROC:-16} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
