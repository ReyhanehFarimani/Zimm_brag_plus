#!/bin/bash
# PERIODIC-BOX J scan with the TABLE interaction (user 2026-09-24: "run a box of density 0.7 of 500 50-segment chains
# for J 9 8 7 6 ... density is the density of monomers so for our model is 0.1"): NCH free chains of NRES residues
# (one residue = a block of 7 fine monomers, so monomer density 0.7 = residue density 0.1) in a cubic periodic box,
# nb_hh = db (six-argument helix-helix table + registry twist term), E_HH = -J, E_CH = 0, E_RL = +J, E_helix = 0.
# Same generator pattern as star_gen.sh: bonds, bends, coil core from runs/<geom>/sample_data.dat; pivots one per
# chain per sweep, flips one per chain per sweep, twists one per residue per sweep; nl_skin = 3 (see
# neighbour-list notes of 2026-09-24: with skin 1 the list was rebuilt ~600 times per sweep).
#   usage: cd runs/<geom>/box500 && [NPROC=4] [SEEDS="1"] [JS="9 8 7 6"] [NCH=500] [NRES=50] [RHO=0.1] \
#          [NEQ=1000] [NSW=10000] [MAXDISP=0.3] [MAXROT=0.3] [NLSKIN=3.0] [SEGLEN=3] [RESTART_DIR=prev_run2] [SUFFIX=r3] ../../box_gen.sh
# Restartable: finished runs (summary in the log) are skipped, running ones are protected by logs/<b>.lock.
GEOM=$(basename "$(cd .. && pwd)"); TH=${GEOM#t}; [ -f ../sample_data.dat ] || { echo "run from runs/<geom>/box<n>"; exit 1; }
BIN=${BIN:-../../../zimm}
DB=/home/reyhaneh/Documents/Zimm_brag_plus/helix_pair_db_${GEOM}_win.bin
JS=${JS:-"9 8 7 6"}
SEEDS=${SEEDS:-"1"}
NCH=${NCH:-500}; NRES=${NRES:-50}; RHO=${RHO:-0.1}
BOX=$(python3 -c "print(round(($NCH * $NRES / $RHO) ** (1.0 / 3.0), 3))")     # residue density RHO -> box side
NEQ=${NEQ:-1000}; NSW=${NSW:-10000}
# move steps for a DENSE box (2026-09-24 04:50: with the dilute-chain steps max_disp 1.0 / max_rot pi the box accepted
# 3.5 % of the position moves and 0.3 % of the pivots and sat frozen 8 % under the 1D helicity): smaller steps
MAXDISP=${MAXDISP:-0.3}; MAXROT=${MAXROT:-0.3}
# segment moves (2026-09-24, user: "segment moves ... flip segments ... torsion moves"): crankshaft rotations of
# 1..SEGLEN beads (one per residue per sweep), handedness flips of helical segments and bond-axis torsions of tails
# (one per chain per sweep each); skin 3 keeps a 3-bead crankshaft of +-0.5 rad inside the neighbour list
SEGLEN=${SEGLEN:-3}; SEGROT=${SEGROT:-0.5}; SEGFLIPLEN=${SEGFLIPLEN:-8}; TORROT=${TORROT:-0.5}
NLSKIN=${NLSKIN:-3.0}
RESTART_DIR=${RESTART_DIR:-}                        # set to a previous run dir (with out/J<J>_s1*_conf.xyz) to continue from its last frames
mkdir -p inputs out logs
for J in $JS; do for s in $SEEDS; do
b=J${J}_s${s}${SUFFIX:-}
[ -f inputs/$b.dat ] && continue
sed -e "s/^N  .*/N           = $NRES             # residues PER CHAIN/" \
    -e "s/^J0 .*/J0          = $J/" -e "s/^J1 .*/J1          = 0/" -e "s/^J2 .*/J2          = $J/" \
    -e 's/^nb_hh.*/nb_hh = db/' -e "s/^n_equil.*/n_equil     = $NEQ/" -e "s/^n_sweeps.*/n_sweeps    = $NSW/" \
    -e 's/^log_every.*/log_every   = 10/' -e 's/^dump_every.*/dump_every  = 1000/' \
    -e "s/^seed .*/seed        = $((${TH}00000 + 200000 + 1000*${J/./} + s))/" \
    -e "s#^out_prefix .*#out_prefix  = out/$b#" \
    -e "s/^n_pivot.*/n_pivot = $NCH                # one pivot attempt per chain per sweep/" \
    -e "s/^n_flip.*/n_flip = $NCH                 # one flip attempt per chain per sweep/" \
    -e "s/^max_disp.*/max_disp    = $MAXDISP             # dense box: small steps (see box_gen.sh)/" ../sample_data.dat > inputs/$b.dat
cat >> inputs/$b.dat <<EOT
# ---- periodic box: $NCH free chains x $NRES residues, residue density $RHO (monomer density $(python3 -c "print(7 * $RHO)"), 7 monomers per residue), box $BOX a (user 2026-09-24) ----
n_arms = $NCH
box = $BOX
E_helix = 0
db_file = $DB
n_twist = $((NCH * NRES))
twist_step = 30
nl_skin = $NLSKIN
max_rot = $MAXROT
n_seg = $((NCH * NRES))
seg_len_max = $SEGLEN
seg_rot = $SEGROT
n_segflip = $NCH
segflip_len_max = $SEGFLIPLEN
n_torsion = $NCH
torsion_rot = $TORROT
EOT
if [ -n "$RESTART_DIR" ]; then
  prev=$(ls -t $RESTART_DIR/out/J${J}_s*_conf.xyz 2>/dev/null | head -1)
  [ -n "$prev" ] && { printf 'restart_file = %s\nrestart_frame = -1\n' "$(cd "$(dirname "$prev")" && pwd)/$(basename "$prev")" >> inputs/$b.dat; sed -i 's/^n_equil .*/n_equil     = 0                # restarted from the last frame of the previous run/' inputs/$b.dat; }
fi
done; done
echo "inputs: $(ls inputs/J*.dat | wc -l)  (box $BOX a)"
[ -n "$GEN_ONLY" ] && exit 0
ls inputs/J*${SUFFIX:-}.dat | LC_ALL=C sort -t_ -k2,2 -k1,1r | xargs -P ${NPROC:-4} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; mkdir logs/$b.lock 2>/dev/null || exit 0; '"$BIN"' {} > logs/$b.log 2>&1; rmdir logs/$b.lock'
echo "done: $(grep -l summary logs/J*.log 2>/dev/null | wc -l) / $(ls inputs/J*.dat | wc -l)"
