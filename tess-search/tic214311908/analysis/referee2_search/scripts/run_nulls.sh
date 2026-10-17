#!/bin/bash
# referee-2 null runs for r03 (sequential, 2 numba threads)
cd tic214311908/analysis/referee2_search
export OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/home/user/work/venv/bin/python
$PY scripts/r03_bls.py null blk 24 11 > logs/r03_null_blk.log 2>&1
$PY scripts/r03_bls.py null blki 16 12 > logs/r03_null_blki.log 2>&1
$PY scripts/r03_bls.py null shift 10 13 > logs/r03_null_shift.log 2>&1
echo done > logs/run_nulls.done
