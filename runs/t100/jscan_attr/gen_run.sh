#!/bin/bash
# t100: a WEAK attraction between same-handed helices, opposite-handed pairs still purely repulsive (per user 2026-09-17).
# In this model that is the chiral amplitude just above the point where the same-handed pair energy first turns negative
# (clamped potential, default range sigma_c = 0.70 a).  Deepest same-handed well, always at r = 1.5 a and one sign of psi:
#     eps_s in {15 16 18}  ->  about {-0.3, -0.6, -1.0} kT;   R.L deepest energy stays > 0 for all of them.
# NOTE: this is beyond the measured interaction (eps_s = 1): a model study of what a same-handed attraction would do.
#   J in {3.5 4 4.5 5} (helix-coil midpoint / handedness crossover), 3 seeds.  Reference: the eps_s = 0 runs of ../jscan_*.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}
JS=${JS:-"3.5 4 4.5 5"}; ESS=${ESS:-"15 16 18"}; SEEDS=${SEEDS:-"1 2 3"}; NEQ=${NEQ:-100000}; NSW=${NSW:-1000000}
mkdir -p inputs out logs; n=0
for J in $JS; do for ES in $ESS; do for s in 1 2 3 4 5 6 7 8 9; do
n=$((n+1)); case " $SEEDS " in *" $s "*) ;; *) continue;; esac
b=J${J}_es${ES}_sc0.70_s${s}; [ -f inputs/$b.dat ] && continue
sed -e "s/^J0 .*/J0          = $J/" -e "s/^J1 .*/J1          = 0/" -e "s/^J2 .*/J2          = $J/" -e 's/^nb_hh.*/nb_hh = fit/' \
    -e "s/^n_equil.*/n_equil     = $NEQ/" -e "s/^n_sweeps.*/n_sweeps    = $NSW/" -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((1400000 + n))/" -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
printf "E_helix = 0\nhf_theta0 = 100\nhf_len = 1.0\nhf_clamp = 1\nhf_eps_s = $ES\nhf_sig_c = 0.70\nn_flip = 4\n" >> inputs/$b.dat
done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*.dat | LC_ALL=C sort -t_ -k4,4 -k2.3,2nr -k1,1r | xargs -P ${NPROC:-20} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
