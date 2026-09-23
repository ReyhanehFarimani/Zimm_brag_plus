#!/bin/bash
# Gyration-tensor test of the tabulated helix-helix potential (nb_hh = db, HELIX_PAIR_DB v2,
# theta0 = 45 table with the registry + twist move) against the measured radial fit (nb_hh = fit)
# that the theta0 = 45 campaign used before.  N = 40 chains with FROZEN state patterns:
#   allL     40 L                  one helical run: the registry is one angle, twist moves rotate it
#   allR     40 R                  mirror of allL -- must agree statistically (table mirror check)
#   halfL_C  20 L + 20 C           helix-coil
#   halfL_R  20 L + 20 R           two runs of opposite hand: the R.L block of the table
# Bonds, bends and the coil core are the theta0 = 45 measured fits (runs/t45/sample_data.dat).
#   usage: [NPROC=4] [POT="db fit"] [SEEDS="1 2 3"] [DB=...] ./gen_gyration_db.sh
cd "$(dirname "$0")"
BIN=${BIN:-../../zimm}
DB=${DB:-/home/reyhaneh/Documents/Zimm_brag_plus/helix_pair_db_t45_win.bin}
POT=${POT:-"db fit"}
SEEDS=${SEEDS:-"1 2 3"}
NSW=${NSW:-500000}
NEQ=${NEQ:-50000}
mkdir -p inputs out logs
mk() {  # tag pattern potential seed
  local name=$1_$3_s$4
  local ntw=0; [ "$3" = db ] && ntw=40          # twist moves only change the energy of the table
  sed -e "s/^N .*/N           = 40/" -e "s/^init .*/init        = walk/" \
      -e "s/^n_equil .*/n_equil     = $NEQ/" -e "s/^n_sweeps .*/n_sweeps    = $NSW/" \
      -e "s/^log_every .*/log_every   = 100/" -e "s/^dump_every .*/dump_every  = 100/" \
      -e "s/^seed .*/seed        = 45$4$4$4/" -e "s#^out_prefix .*#out_prefix  = out/$name#" \
      -e "s/^nb_hh .*/nb_hh = $3/" ../t45/sample_data.dat > inputs/$name.dat
  cat >> inputs/$name.dat <<INP
# ---- gyration test: frozen pattern, registry twist moves (used by nb_hh = db only) ----
state_pattern = $2
freeze_states = 1
db_file = $DB
n_twist = $ntw
twist_step = 30
n_pivot = 4
nl_skin = 1.0
INP
}
L20=LLLLLLLLLLLLLLLLLLLL; R20=RRRRRRRRRRRRRRRRRRRR; C20=CCCCCCCCCCCCCCCCCCCC
for p in $POT; do for s in $SEEDS; do
  mk allL    $L20$L20 $p $s
  mk allR    $R20$R20 $p $s
  mk halfL_C $L20$C20 $p $s
  mk halfL_R $L20$R20 $p $s
done; done
ls inputs/*.dat | xargs -P ${NPROC:-4} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; nice -n 5 '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
