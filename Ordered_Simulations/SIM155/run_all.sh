#!/bin/bash
# SIM155: all best fits, D1 first. 4 parallel jobs x 3 OpenMP threads (12 cores).
cd "$(dirname "$0")"
for d in D1 D2 D3; do for m in lcdm w0wa cmstgA; do for s in 1 2 3; do echo "$m $d $s"; done; done; done |
  xargs -P 4 -L 1 bash -c 'OMP_NUM_THREADS=3 /home/aion/cmstg-venv/bin/python run_fit.py $0 $1 $2 > fits/$0_$1_s$2.log 2>&1; echo "done $0 $1 s$2: $(grep -o "\"chi2_total\": [0-9.]*" fits/$0_$1_s$2_result.json 2>/dev/null)"'
echo ALL_DONE
