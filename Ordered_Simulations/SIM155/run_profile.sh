#!/bin/bash
# SIM155 profile likelihood (exploratory method for locating the CMSTG best fit; see profile.py).
# Fills CPU slots as the main fits finish: at most 6 fit/profile python processes at once (anchored match, so the
# xargs wrappers are not counted), 2 OpenMP threads each.
# Usage: ./run_profile.sh DATA dF1 dF2 ...
cd "$(dirname "$0")"
mkdir -p profile
data=$1; shift
for dF in "$@"; do
  while [ "$(pgrep -c -f '^/home/aion/cmstg-venv/bin/python (fit|profile)\.py')" -ge 6 ]; do sleep 20; done
  echo "start $data dF=$dF $(date +%H:%M)"
  OMP_NUM_THREADS=2 nohup /home/aion/cmstg-venv/bin/python profile.py $data $dF > profile/cmstgA_${data}_dF$(printf %.4f $dF).log 2>&1 &
  sleep 5
done
wait
echo PROFILE_DONE
