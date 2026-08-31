#!/bin/bash
# Progress of the batch: finished runs, running processes, and per-run % of production sweeps done.
cd "$(dirname "$0")"
echo "finished: $(grep -l summary logs/*.log 2>/dev/null | wc -l | tr -d ' ') / $(ls inputs | wc -l | tr -d ' ')   running zimm processes: $(pgrep -x zimm | wc -l | tr -d ' ')"
for f in inputs/*.dat; do
  b=$(basename $f .dat); o=out/${b}_obs.dat; [ -f "$o" ] || continue
  nsw=$(awk '/^n_sweeps/{print $3}' $f); last=$(tail -1 $o | awk '!/^#/{print $1}')
  [ -n "$last" ] && printf "  %-10s %5.1f%%\n" $b $(echo "100*$last/$nsw" | bc -l)
done
