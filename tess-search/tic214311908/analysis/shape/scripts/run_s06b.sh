#!/bin/bash
# s06 remaining scenarios (second pass, after the WD-radius and secondary-dilution fixes)
cd tic214311908/analysis/shape
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/home/user/work/venv/bin/python
for a in "ref twin" "inv twin" "ref lateM" "ref any" "ref single_ecc" "inv lateM"; do
  $PY scripts/s06_compact.py $a >> out/s06.log 2>&1
done
