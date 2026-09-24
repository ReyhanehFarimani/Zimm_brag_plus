#!/bin/bash
# STAR J scan with the TABLE interaction (user 2026-09-23: "run the star50 ... consider longer run ... use 20 cpu ...
# 3 sample each"): n_arms arms of N residues on a spherical core (core_radius), nb_hh = db (six-argument helix-helix
# table + registry twist term), E_HH = -J, E_CH = 0, E_RL = +J, E_helix = 0.  Same generator pattern as jscan_db_gen.sh:
# bonds, bends, coil core from runs/<geom>/sample_data.dat; pivots ON (NPIV, one per arm per sweep); n_twist = one per residue per sweep,
# n_flip = one per ARM per sweep (the linear scan has 4 flips per chain of 200; the flip picks a residue uniformly, so
# a star of 50 arms needs ~50 to give every arm a flip attempt per sweep -- a sampling move, not physics).
#   usage: cd runs/<geom>/star50 && [NPROC=20] [SEEDS="1 2 3"] [JS="..."] [NARMS=50] [NRES=50] [CORE=2.5] \
#          [NEQ=2000] [NSW=30000] [NPIV=50] ../../star_gen.sh
# Restartable: finished runs (summary in the log) are skipped, running ones are protected by logs/<b>.lock.
GEOM=$(basename "$(cd .. && pwd)"); TH=${GEOM#t}; [ -f ../sample_data.dat ] || { echo "run from runs/<geom>/star<n>"; exit 1; }
BIN=${BIN:-../../../zimm}
DB=/home/reyhaneh/Documents/Zimm_brag_plus/helix_pair_db_${GEOM}_win.bin
JS=${JS:-"1 2 3 4 4.5 5 5.5 6 6.5 7 7.5 8 9"}     # 13 values x 3 seeds = 39 runs: two waves on 20 CPUs
SEEDS=${SEEDS:-"1 2 3"}
NARMS=${NARMS:-50}; NRES=${NRES:-50}; CORE=${CORE:-2.5}
NPIV=${NPIV:-$NARMS}          # pivots: one attempt per arm per sweep (2026-09-23: without pivots the coil kinks between helical
                              # segments freeze in the crowded corona; helicity plateaued 3 % under the 1D theory at J = 9, with pivots it reaches it)
NEQ=${NEQ:-2000}; NSW=${NSW:-30000}                 # ~1.9 s/sweep for 50 x 50 -> ~17 h per run
mkdir -p inputs out logs
for J in $JS; do for s in $SEEDS; do
b=J${J}_s${s}
[ -f inputs/$b.dat ] && continue
sed -e "s/^N  .*/N           = $NRES             # residues PER ARM/" \
    -e "s/^J0 .*/J0          = $J/" -e "s/^J1 .*/J1          = 0/" -e "s/^J2 .*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = db/' -e "s/^n_equil.*/n_equil     = $NEQ/" -e "s/^n_sweeps.*/n_sweeps    = $NSW/" \
    -e 's/^log_every.*/log_every   = 50/' -e 's/^dump_every.*/dump_every  = 2000/' \
    -e "s/^seed .*/seed        = $((${TH}00000 + 100000 + 1000*${J/./} + s))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" \
    -e "s/^n_pivot.*/n_pivot = $NPIV                # one pivot attempt per arm per sweep (see star_gen.sh)/" \
    -e "s/^n_flip.*/n_flip = $NARMS                 # one flip attempt per arm per sweep/" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
# ---- star: $NARMS arms x $NRES residues on a core of radius $CORE a, table + twist terms (user 2026-09-23) ----
n_arms = $NARMS
core_radius = $CORE
E_helix = 0
db_file = $DB
n_twist = $((NARMS * NRES))
twist_step = 30
nl_skin = 1.0
EOT
done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*.dat | LC_ALL=C sort -t_ -k2,2 -k1,1r | xargs -P ${NPROC:-20} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
