#!/bin/bash
# Parallel TRICERATOPS stream 3 (1 numba thread): the E variants and the D runs in reverse seed order; runs already done
# or running in stream 1 (run_tri_batch.sh) are skipped by f06 (it checks for out/tri/<tag>.json or .running).
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
PY=/home/user/work/venv_tri/bin/python
S=tic214311908/analysis/fpp/scripts/f06_tri_run.py
L=tic214311908/analysis/fpp/logs
run() { tag=$1; shift; echo "$(date) start $tag" >> $L/tri_batch4.log; $PY $S "$@" $tag > $L/tri4_$tag.log 2>&1; echo "$(date) end $tag rc=$?" >> $L/tri_batch4.log; }
for s in 1 2 3; do run E1_s$s ref200 $s 200000 cons loc 3 emp; done
for s in 1 2 3; do run E2_s$s inv600 $s 200000 cons loc 3 emp; done
for s in 1 2 3; do run E3_s$s inv200 $s 200000 opt loc 3 emp; done
for s in 1 2 3; do run E4_s$s inv200 $s 200000 cons loc 3 prop; done
for s in 1 2; do run E5_s$s inv200 $s 30000 cons keep 1.5 emp; done
for s in 10 9 8 7 6 5 4 3 2 1; do run D_s$s inv200 $s 30000 cons keep 3 emp; done
echo "$(date) ALL DONE" >> $L/tri_batch4.log
