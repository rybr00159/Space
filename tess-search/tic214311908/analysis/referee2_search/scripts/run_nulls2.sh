#!/bin/bash
# referee-2 null runs for r03, second batch (the first batch, blk seed 11, was stopped after 4 realisations to
# share the CPUs; its 4 saved realisations are used)
cd tic214311908/analysis/referee2_search
export OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/home/user/work/venv/bin/python
$PY scripts/r03_bls.py null blki 10 12 > logs/r03_null_blki.log 2>&1
$PY scripts/r03_bls.py null shift 6 13 > logs/r03_null_shift.log 2>&1
$PY scripts/r03_bls.py null blk 12 21 > logs/r03_null_blk21.log 2>&1
echo done > logs/run_nulls2.done
