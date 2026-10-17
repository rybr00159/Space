#!/bin/bash
# s12 background-EB scenario fits, sequential (one process)
cd tic214311908/analysis/shape
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/home/user/work/venv/bin/python
for src in bg3 bg5; do for cls in cwd csdb celm msms5 msms11 msbd; do
  $PY scripts/s12_bg_eb.py ref $src $cls >> out/s12.log 2>&1
done; done
for src in bg3 bg5; do for cls in cwd celm; do $PY scripts/s12_bg_eb.py ref $src $cls lconly >> out/s12.log 2>&1; done; done
for src in bg3 bg5; do $PY scripts/s12_bg_eb.py inv $src cwd >> out/s12.log 2>&1; $PY scripts/s12_bg_eb.py inv $src celm >> out/s12.log 2>&1; done
