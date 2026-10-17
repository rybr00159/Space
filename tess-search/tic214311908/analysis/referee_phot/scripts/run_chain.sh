#!/bin/bash
# Re-runs the local (no-network) referee chain in order, one log per script.
cd tic214311908/analysis/referee_phot
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
PY=/home/user/work/venv/bin/python
for s in r01_wcs r01b_fluxscale r02_time r02b_qlpstamps r03_aux r04_mirab r04b_mirastatic r05_local r06_transit r07_pertransit r08_oddeven_timing r09_xpflux r10_diffimg r11_compare r12_centroid_inject r13_csvcheck r14_shape_crosscheck r15_timing_dist; do
  echo "== $s $(date +%T)"; $PY scripts/$s.py > logs/$s.log 2>&1 || echo "FAILED $s"
done
echo "== done $(date +%T)"
