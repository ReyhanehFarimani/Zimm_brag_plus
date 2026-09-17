#!/bin/bash
# Master chain (2026-09-17): each stage starts when the previous one is complete and no zimm is running. 20 cores.
#   1 t45/jscan_sigc (running; its own helpers finish the reference, control and extra-J runs)
#   2 t45/jscan_N400     finite-size test of the suspicious crossover cells
#   3 t100/jscan_eps     corrected re-run, hf_len = 1
#   4 t100/jscan_sigc    corrected re-run, hf_len = 1
#   + t45/jscan_attr after 2 and t100/jscan_attr after 4: weak same-handed attraction
#   5 t45/jscan_scale, 6 t100/jscan_scale   whole-potential amplification
cd "$(dirname "$0")"
complete () { [ "$(grep -l summary $1/logs/J*.log 2>/dev/null | wc -l)" -ge "$(ls $1/inputs/J*.dat 2>/dev/null | wc -l)" ]; }
wait_for () { until complete $1 && ! pgrep -x zimm > /dev/null; do sleep 30; done; echo "$(date '+%F %T') $1 complete"; }
stage () { echo "$(date '+%F %T') start $1"; (cd $1 && NPROC=20 ./gen_run.sh > batch.log 2>&1); wait_for $1; }
wait_for t45/jscan_sigc
for d in t45/jscan_N400 t45/jscan_attr t100/jscan_eps t100/jscan_sigc t100/jscan_attr t45/jscan_scale t100/jscan_scale; do stage $d; done
echo "$(date '+%F %T') chain finished"
