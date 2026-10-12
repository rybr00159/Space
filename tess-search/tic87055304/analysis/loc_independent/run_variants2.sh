#!/bin/bash
cd /home/user/work/tic87/loc_independent
PY=/home/user/work/venv/bin/python
export OMP_NUM_THREADS=1
(LOC_PLANE=1 nice -n 10 $PY s02_diffimg.py prf 3 local_r4 prf_plane > log_var_prf_plane.txt 2>&1) &
(LOC_CLEAN=1 nice -n 10 $PY s02_diffimg.py prf 3 local_r4 prf_clean > log_var_prf_clean.txt 2>&1; LOC_CLEAN=1 LOC_PLANE=1 nice -n 10 $PY s02_diffimg.py prf 3 local_r4 prf_clean_plane > log_var_prf_clean_plane.txt 2>&1) &
wait
echo ALLDONE > log_variants2_done.txt
