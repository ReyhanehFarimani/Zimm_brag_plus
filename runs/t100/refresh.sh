#!/bin/bash
# Refresh the jscan_eps and jscan_sigc analyses (finished runs only) and print the lines worth reading.
# Figures: jscan_eps/transitions.png, jscan_sigc/range_scan.png.  Full tables: the two .out files.
cd "$(dirname "$0")"
PY=${PY:-/home/reyhaneh/.conda/envs/sim_analysis/bin/python}
date "+%Y-%m-%d %H:%M:%S"
echo "running zimm: $(pgrep -c -x zimm)   jscan_eps: $(grep -l summary jscan_eps/logs/J*.log 2>/dev/null | wc -l)/$(ls jscan_eps/inputs/J*.dat | wc -l)   jscan_sigc: $(grep -l summary jscan_sigc/logs/J*.log 2>/dev/null | wc -l)/$(ls jscan_sigc/inputs/J*.dat | wc -l)   errors: $(grep -l -i error jscan_eps/logs/J*.log jscan_sigc/logs/J*.log 2>/dev/null | wc -l)"
NEPS=$(grep -l summary jscan_eps/logs/J*_es*.log 2>/dev/null | wc -l); NSIG=$(grep -l summary jscan_sigc/logs/J*.log 2>/dev/null | wc -l)
echo "t100 (hf_len = 1 re-run; the hf_len = 1.2311 data is in prev_hflen1.2311/): eps scan $NEPS/108, range scan $NSIG/180"
[ "$NEPS" -gt 0 ] && (cd jscan_eps && $PY plot_transitions.py > plot_transitions.out 2>&1; sed -n '1p;/WHERE DOES/,/RMS z:  control.*<|m|>/p' plot_transitions.out | head -24; grep -E "all J combined" plot_transitions.out)
[ "$NSIG" -gt 0 ] && (cd jscan_sigc && $PY analyze_plot.py > analyze_plot.out 2>&1; grep -E "^range scan|^== |sigma_c = |RMS z|at the largest|Traceback|Error" analyze_plot.out)
echo "--- t100 cells beyond 3 sigma in the range scan:"
[ "$NSIG" -gt 0 ] && grep -oE "J = [0-9]: .*" jscan_sigc/analyze_plot.out | grep -E "\[[+-]([3-9]|[1-9][0-9])\.[0-9]\]" | head -10 || true
[ "$NSIG" -gt 0 ] && (cd jscan_sigc && $PY plot_transitions.py > plot_transitions.out 2>&1; head -1 plot_transitions.out)
echo "=================== theta0 = 45 scan (../t45/jscan_sigc): finished $(grep -l summary ../t45/jscan_sigc/logs/J*.log 2>/dev/null | wc -l)/$(ls ../t45/jscan_sigc/inputs/J*.dat 2>/dev/null | wc -l), running t45 sigc/attr: $(find ../t45/jscan_sigc/logs -maxdepth 1 -name "*.lock" ! -exec test -e {}/DROPPED ; -print 2>/dev/null | wc -l)/$(ls -d ../t45/jscan_attr/logs/*.lock 2>/dev/null | wc -l), errors: $(grep -l -i error ../t45/jscan_sigc/logs/J*.log 2>/dev/null | wc -l)"
(cd ../t45/jscan_sigc && $PY analyze.py > analyze.out 2>&1; grep -E "^t45 scan|^== |RMS z|Traceback|Error" analyze.out; sed -n '/MC vs exact/,$p' analyze.out
 echo "--- t45 cells beyond 3 sigma:"; grep -oE "sigma_c [0-9.]+  J = [0-9]: .*" analyze.out | grep -E "\[[+-]([3-9]|[1-9][0-9])\.[0-9]\]" | head -10 || true)
(cd ../t45/jscan_attr && $PY zoom.py > zoom.out 2>&1; echo "=== ZOOM (strongest attraction, J = 2.5 .. 5):"; cat zoom.out | cut -c1-400)
exit 0
