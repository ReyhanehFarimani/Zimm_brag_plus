#!/bin/bash
# Gyration tensor of N = 40 chains with frozen state patterns: all-L helix, half L helix + half coil,
# half L + half R. User's bonds/bends (theta0 = 180 deg), Gaussian coil core + Gay-Berne helix rods.
cd "$(dirname "$0")"
BIN=${BIN:-../../zimm_dev}
mkdir -p inputs out logs
mk() {  # tag pattern seed
sed -e "s/^N .*/N           = 40/" -e "s/^init .*/init        = walk/" -e "s/^max_disp .*/max_disp    = 1.0/" \
    -e "s/^n_equil .*/n_equil     = 50000/" -e "s/^n_sweeps .*/n_sweeps    = 300000/" -e "s/^log_every .*/log_every   = 100/" \
    -e "s/^dump_every .*/dump_every  = 100/" -e "s/^seed .*/seed        = $3/" -e "s/^\(theta0_[A-Z]*\) .*/\1 = 180 deg/" \
    -e "s#^out_prefix .*#out_prefix  = out/$1_s$3#" ../no_steric/data_sample.dat > inputs/$1_s$3.dat
cat >> inputs/$1_s$3.dat <<INP
state_pattern = $2
freeze_states = 1
nb_type = gauss
nb_A = 9.01
nb_sigma = 1.5
nb_hh = gb
gb_eps0 = 1.0
rod_L = 2.5
rod_r = 1.0
n_pivot = 4
INP
}
L20=LLLLLLLLLLLLLLLLLLLL; R20=RRRRRRRRRRRRRRRRRRRR; C20=CCCCCCCCCCCCCCCCCCCC
for s in 1 2 3; do
  mk allL      $L20$L20 $s
  mk halfL_C   $L20$C20 $s
  mk halfL_R   $L20$R20 $s
done
ls inputs/*.dat | xargs -P ${NPROC:-6} -I{} sh -c 'b=$(basename {} .dat); grep -q summary logs/$b.log 2>/dev/null && exit 0; '"$BIN"' {} > logs/$b.log 2>&1'
echo "done: $(grep -l summary logs/*.log | wc -l) / $(ls inputs | wc -l)"
