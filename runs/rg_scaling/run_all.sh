#!/bin/bash
# Run every input in parallel (one process per core). Logs go to logs/, data to out/.
cd "$(dirname "$0")"
NPROC=${NPROC:-$(sysctl -n hw.ncpu)}
# smallest N first: the cheap runs finish within the hour and let us check tau_int before the
# multi-day N=4096 runs are far along (they start ~20 min later either way, and are the tail regardless)
ls inputs/*.dat | sort -t N -k2 -n | xargs -P "$NPROC" -I{} sh -c 'b=$(basename {} .dat); ../../zimm {} > logs/$b.log 2>&1'
echo "finished: $(grep -l 'summary' logs/*.log | wc -l) / $(ls inputs | wc -l)"
