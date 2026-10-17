#!/bin/bash
# s02 null calibration runs (sequential, 2 numba threads); results saved after every realisation
cd tic214311908/analysis/search
export OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/home/user/work/venv/bin/python
$PY scripts/s02_bls_null.py ctrl 1 0 clean > logs/s02_ctrl.log 2>&1
for k in pair_A pair_A_inv Mira_B Mira_B_inv; do $PY scripts/s02_bls_null.py $k 1 0 clean > logs/s02_$k.log 2>&1; done
$PY scripts/s02_bls_null.py roll 30 2 clean > logs/s02_roll_2.log 2>&1
$PY scripts/s02_bls_null.py inv_roll 40 3 clean > logs/s02_inv_roll_3.log 2>&1
