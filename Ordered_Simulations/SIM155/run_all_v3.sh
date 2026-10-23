#!/bin/bash
# SIM155 v3: all 27 best fits (3 models x D1-D3 x 3 seeds), longest (CMSTG) and D1 first.
# 6 parallel jobs x 2 OpenMP threads (12 cores). Each fit writes fits_v3/<model>_<data>_s<seed>.{json,log}.
cd "$(dirname "$0")"
mkdir -p fits_v3
for ms in "cmstgA D1" "lcdm D1" "w0wa D1" "cmstgA D2" "cmstgA D3" "lcdm D2" "lcdm D3" "w0wa D2" "w0wa D3"; do
  for s in 1 2 3; do echo "$ms $s"; done
done | xargs -P 6 -L 1 bash -c 'OMP_NUM_THREADS=2 /home/aion/cmstg-venv/bin/python fit.py $0 $1 $2 > fits_v3/$0_$1_s$2.log 2>&1; echo "done $0 $1 s$2 exit $?"'
echo ALL_DONE
