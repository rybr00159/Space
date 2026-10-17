#!/bin/bash
# second referee-2 TRICERATOPS batch: REF light curve with the analyst's conservative contrast curve (isolates the cc effect)
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
cd tic214311908/analysis/referee2_fpp
until [ -f logs/tri_r2_batch.done ]; do sleep 20; done
for s in 101 102 103 104; do
  timeout 1500 /home/user/work/venv_tri/bin/python scripts/r06b_tri_run.py $s 200000 emp cccons > logs/tri_r2_emp_cccons_s$s.log 2>&1
done
echo done > logs/tri_r2b_batch.done
