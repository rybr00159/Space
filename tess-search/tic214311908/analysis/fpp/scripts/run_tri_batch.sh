#!/bin/bash
# Sequential TRICERATOPS run set (f06_tri_run.py). Logs in logs/tri_*.log, results in out/tri/*.json.
# A  main     : inv200, conservative IPD/DECaPS2 contrast curve, nearby = loc (3.1" and 5.5" stars only), N=2e5, 12 seeds
# C  no CC    : as A without contrast curve, 10 seeds
# D  raw      : as A but all TIC neighbours kept (TRICERATOPS default treatment), N=3e4, 10 seeds
# E1-E5       : REF light curve; 600-s binning incl. S39; optimistic contrast curve; propagated errors; 1.5-px aperture raw
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/home/user/work/venv_tri/bin/python
S=tic214311908/analysis/fpp/scripts/f06_tri_run.py
L=tic214311908/analysis/fpp/logs
run() { tag=$1; shift; echo "$(date) start $tag" >> $L/tri_batch.log; $PY $S "$@" $tag > $L/tri_$tag.log 2>&1; echo "$(date) end $tag rc=$?" >> $L/tri_batch.log; }
for s in 1 2 3 4 5 6 7 8 9 10 11 12; do run A_s$s inv200 $s 200000 cons loc 3 emp; done
for s in 1 2 3 4 5 6 7 8 9 10; do run C_s$s inv200 $s 200000 none loc 3 emp; done
for s in 1 2 3 4 5 6 7 8 9 10; do run D_s$s inv200 $s 30000 cons keep 3 emp; done
for s in 1 2 3; do run E1_s$s ref200 $s 200000 cons loc 3 emp; done
for s in 1 2 3; do run E2_s$s inv600 $s 200000 cons loc 3 emp; done
for s in 1 2 3; do run E3_s$s inv200 $s 200000 opt loc 3 emp; done
for s in 1 2 3; do run E4_s$s inv200 $s 200000 cons loc 3 prop; done
for s in 1 2; do run E5_s$s inv200 $s 30000 cons keep 1.5 emp; done
echo "$(date) ALL DONE" >> $L/tri_batch.log
