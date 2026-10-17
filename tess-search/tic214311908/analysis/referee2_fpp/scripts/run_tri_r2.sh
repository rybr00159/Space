#!/bin/bash
# referee-2 TRICERATOPS runs (sequential, 2 numba threads)
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
cd tic214311908/analysis/referee2_fpp
for s in 101 102 103 104 105 106; do
  timeout 1500 /home/user/work/venv_tri/bin/python scripts/r06b_tri_run.py $s 200000 emp > logs/tri_r2_emp_s$s.log 2>&1
done
for s in 101 102 103; do
  timeout 1500 /home/user/work/venv_tri/bin/python scripts/r06b_tri_run.py $s 200000 prop > logs/tri_r2_prop_s$s.log 2>&1
done
echo done > logs/tri_r2_batch.done
