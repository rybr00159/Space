#!/bin/bash
# runs the s03 fits two at a time (machine limit: 4 CPUs shared)
cd tic214311908/analysis/shape
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
PY=/home/user/work/venv/bin/python
run() { $PY scripts/s03_planet_fit.py $1 $2 > out/s03_$1_$2.log 2>&1; }
for m in "$@"; do
  run inv $m & run ref $m; wait
done
