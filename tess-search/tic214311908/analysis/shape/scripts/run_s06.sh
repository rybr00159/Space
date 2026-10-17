#!/bin/bash
# s06 compact-star occultation scenarios, sequential (one process)
cd tic214311908/analysis/shape
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for sc in single twin single_ecc lateM any; do for ds in ref inv; do
  /home/user/work/venv/bin/python scripts/s06_compact.py $ds $sc >> out/s06.log 2>&1
done; done
/home/user/work/venv/bin/python scripts/s06_compact.py ref single_Rfree >> out/s06.log 2>&1
/home/user/work/venv/bin/python scripts/s06_compact.py ref twin >> out/s06.log 2>&1; /home/user/work/venv/bin/python scripts/s06_compact.py inv twin >> out/s06.log 2>&1
