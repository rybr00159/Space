#!/bin/sh
# Simulation batch 2: uncatalogued-source ring injections around the M3V, PSF-mismatch injections (alt/blur/emp,
# empirical Mira PSF), and injections with each hypothesis's own best-fit amplitude.
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
PY=/home/user/work/venv/bin/python
cd tic214311908/analysis/loc/scripts
$PY r01_sim.py DIFF prf 1200 600 ring diff_ring > ../out/log_sim_diff_ring.txt 2>&1
$PY r01_sim.py AMP prf 1200 600 ring amp_ring > ../out/log_sim_amp_ring.txt 2>&1
for m in alt blur emp; do
  $PY r01_sim.py DIFF $m 1200 1000 main diff_${m}_main > ../out/log_sim_diff_${m}_main.txt 2>&1
  $PY r01_sim.py AMP $m 1200 1000 main amp_${m}_main > ../out/log_sim_amp_${m}_main.txt 2>&1
done
$PY r01_sim.py DIFF miraemp 1200 1000 mira diff_miraemp > ../out/log_sim_diff_miraemp.txt 2>&1
$PY r01_sim.py AMP miraemp 1200 1000 mira amp_miraemp > ../out/log_sim_amp_miraemp.txt 2>&1
AMPMODE=own $PY r01_sim.py DIFF prf 1200 1000 main diff_prf_main_ownamp > ../out/log_sim_diff_prf_main_ownamp.txt 2>&1
echo BATCH2_DONE
