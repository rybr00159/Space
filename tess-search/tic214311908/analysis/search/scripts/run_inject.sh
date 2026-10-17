#!/bin/bash
# s09 injection-recovery in chunks of 100 (seeds 1-8), 1 thread
cd tic214311908/analysis/search
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for sd in 1 2 3 4 5 6 7 8; do /home/user/work/venv/bin/python scripts/s09_inject.py 100 $sd > logs/s09_inject_$sd.log 2>&1; done
