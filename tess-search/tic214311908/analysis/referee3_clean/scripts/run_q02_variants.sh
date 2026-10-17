#!/bin/bash
# Runs q02 (large-aperture raw-pixel flux losses) robustness variants, then q03/q04 summaries for each.
cd tic214311908/analysis/referee3_clean
PY=/home/user/work/venv/bin/python
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run () { env "$@" $PY scripts/q02_bigap.py && T=$(env "$@" $PY -c 'import os;print(os.environ.get("Q2_TAG",""))') && $PY scripts/q03_bigap_summary.py $T > logs/q03$T.log && $PY scripts/q04_paired.py $T > logs/q04$T.log; }
run Q2_TAG=_bkgnone Q2_BKG=none
run Q2_TAG=_bkgwide Q2_BKG=wide
run Q2_TAG=_ord1 Q2_ORDER=1
run Q2_TAG=_half40 Q2_HALF=0.40
run Q2_TAG=_half25o1 Q2_HALF=0.25 Q2_ORDER=1
run Q2_TAG=_q2048 Q2_QUAL=q0_2048
echo done
