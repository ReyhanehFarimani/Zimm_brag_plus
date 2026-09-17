#!/bin/bash
# theta0 = 45, chain length N = 400: finite-size test of the cells that looked suspicious at N = 200 (per user 2026-09-17).
# At N = 200 (one seed) <|m|> rose with eps_s in the handedness-crossover region, e.g. J = 5, sigma_c = 0.62 a:
# 0.6856 / 0.6951 / 0.6994 at eps_s = 2.6 / 3.4 / 4.2.  A real collective effect of the chiral term must GROW with N;
# a fluctuation will not repeat.  The exact 1D crossover moves from J* = 4.68 (N = 200) to 4.98 (N = 400), hence
#   J in {4.5, 5, 5.5, 6};   (sigma_c, eps_s) in {(-, 0), (0.62, 1), (0.62, 2.6), (0.62, 4.2), (0.51, 2.9)};   3 seeds   (60 runs)
# Everything else as ../jscan_sigc (same sample_data.dat, hf_len = 1, 1M sweeps).  Cost: 2.6 x the N = 200 run.
cd "$(dirname "$0")"
BIN=${BIN:-../../../zimm}; NN=${NN:-400}
JS=${JS:-"4.5 5 5.5 6"}; SEEDS=${SEEDS:-"1 2 3"}; NEQ=${NEQ:-100000}; NSW=${NSW:-1000000}
CELLS=${CELLS:-"0.51:0 0.62:1 0.62:2.6 0.62:4.2 0.51:2.9"}
mkdir -p inputs out logs; n=0
for J in $JS; do for C in $CELLS; do SC=${C%%:*}; ES=${C##*:}; for s in 1 2 3 4 5 6 7 8 9; do
n=$((n+1)); case " $SEEDS " in *" $s "*) ;; *) continue;; esac
b=J${J}_es${ES}_sc${SC}_s${s}; [ "$ES" = 0 ] && b=J${J}_es0_sc0_s${s}
[ -f inputs/$b.dat ] && continue
sed -e "s/^N  .*/N           = $NN/" -e "s/^J0 .*/J0          = $J/" -e "s/^J1 .*/J1          = 0/" -e "s/^J2 .*/J2          = $J/" \
    -e "s/^n_equil.*/n_equil     = $NEQ/" -e "s/^n_sweeps.*/n_sweeps    = $NSW/" -e 's/^dump_every.*/dump_every  = 50000/' \
    -e "s/^seed .*/seed        = $((4900000 + n))/" -e "s#^out_prefix .*#out_prefix  = out/$b#" ../sample_data.dat > inputs/$b.dat
printf "E_helix = 0\nhf_eps_s = $ES\nhf_sig_c = $SC\n" >> inputs/$b.dat
done; done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*.dat | LC_ALL=C sort -t_ -k4,4 -k1,1r | xargs -P ${NPROC:-20} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
