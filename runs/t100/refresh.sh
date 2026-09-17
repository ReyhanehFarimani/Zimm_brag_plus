#!/bin/bash
# Refresh the jscan_eps and jscan_sigc analyses (finished runs only) and print the lines worth reading.
# Figures: jscan_eps/transitions.png, jscan_sigc/range_scan.png.  Full tables: the two .out files.
cd "$(dirname "$0")"
PY=${PY:-/home/reyhaneh/.conda/envs/sim_analysis/bin/python}
date "+%Y-%m-%d %H:%M:%S"
echo "running zimm: $(pgrep -c -x zimm)   jscan_eps: $(grep -l summary jscan_eps/logs/J*.log 2>/dev/null | wc -l)/$(ls jscan_eps/inputs/J*.dat | wc -l)   jscan_sigc: $(grep -l summary jscan_sigc/logs/J*.log 2>/dev/null | wc -l)/$(ls jscan_sigc/inputs/J*.dat | wc -l)   errors: $(grep -l -i error jscan_eps/logs/J*.log jscan_sigc/logs/J*.log 2>/dev/null | wc -l)"
(cd jscan_sigc && $PY analyze_plot.py > analyze_plot.out 2>&1; grep -E "^range scan|^== |sigma_c = |RMS z|at the largest|Traceback|Error" analyze_plot.out)
echo "--- cells beyond 3 sigma in the range scan:"
grep -oE "J = [0-9]: .*" jscan_sigc/analyze_plot.out | grep -E "\[[+-]([3-9]|[1-9][0-9])\.[0-9]\]" | head -10 || true
[ -n "$WITH_EPS" ] && (cd jscan_eps && $PY plot_transitions.py > plot_transitions.out 2>&1; head -1 plot_transitions.out)
exit 0
