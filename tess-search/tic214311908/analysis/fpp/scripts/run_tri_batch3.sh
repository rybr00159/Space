#!/bin/bash
# Parallel TRICERATOPS stream 2 (network-bound; 1 numba thread): F = A with circular EB orbits (10 seeds), G = C with circular EBs (3)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
PY=/home/user/work/venv_tri/bin/python
S=tic214311908/analysis/fpp/scripts/f06_tri_run.py
L=tic214311908/analysis/fpp/logs
run() { tag=$1; shift; echo "$(date) start $tag" >> $L/tri_batch3.log; TRI_EBCIRC=1 $PY $S "$@" $tag > $L/tri_$tag.log 2>&1; echo "$(date) end $tag rc=$?" >> $L/tri_batch3.log; }
for s in 1 2 3 4 5 6 7 8 9 10; do run F_s$s inv200 $s 200000 cons loc 3 emp; done
for s in 1 2 3; do run G_s$s inv200 $s 200000 none loc 3 emp; done
echo "$(date) ALL DONE" >> $L/tri_batch3.log
