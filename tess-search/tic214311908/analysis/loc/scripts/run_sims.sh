#!/bin/sh
# Simulation batch: baseline (PRF) injections at all 93 catalogued stars, both image methods.
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
PY=/home/user/work/venv/bin/python
cd tic214311908/analysis/loc/scripts
$PY r01_sim.py DIFF prf 1200 1000 all diff_prf_all > ../out/log_sim_diff_prf_all.txt 2>&1
$PY r01_sim.py AMP prf 1200 1000 all amp_prf_all > ../out/log_sim_amp_prf_all.txt 2>&1
