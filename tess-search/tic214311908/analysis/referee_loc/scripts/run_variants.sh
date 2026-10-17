#!/bin/bash
# Referee variants of rq04_fit.py (real data; a few with injections). Run from the referee folder.
cd tic214311908/analysis/referee_loc
PY="/home/user/work/venv/bin/python scripts/rq04_fit.py"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
$PY W0.30 linear 500 prf v_linear OFFKEY=linear_R2.5_F0.03 > out/log_fit_v_linear.txt 2>&1
$PY W0.30 none 500 prf v_nocolour OFFKEY=none_R2.5_F0.03 > out/log_fit_v_nocolour.txt 2>&1
$PY W0.30 curve 300 prf v_cs0.7 CSCALE=0.7 > out/log_fit_v_cs0.7.txt 2>&1
$PY W0.30 curve 300 prf v_cs1.3 CSCALE=1.3 > out/log_fit_v_cs1.3.txt 2>&1
$PY W0.30 curve 0 prf v_locpos OFF=loc MIRA=loc > out/log_fit_v_locpos.txt 2>&1
$PY W0.20 curve 300 prf v_W020 > out/log_fit_v_W020.txt 2>&1
$PY W0.45 curve 300 prf v_W045 > out/log_fit_v_W045.txt 2>&1
$PY W0.30 curve 0 prf v_R3.5 RFIT=3.5 > out/log_fit_v_R3.5.txt 2>&1
$PY W0.30 curve 0 prf v_R7 RFIT=7 > out/log_fit_v_R7.txt 2>&1
$PY W0.30 curve 300 prf v_inflpix INFL=pix > out/log_fit_v_inflpix.txt 2>&1
for s in 12 39 66 93 103 104; do $PY W0.30 curve 300 prf v_drop$s DROP=$s > out/log_fit_v_drop$s.txt 2>&1; done
$PY W0.30 curve 1000 blur v_blur > out/log_fit_v_blur.txt 2>&1
$PY W0.30 curve 1000 prf v_ownamp AMP=own > out/log_fit_v_ownamp.txt 2>&1
$PY W0.30 curve 500 prf v_c175_2.0 C175=2.0 > out/log_fit_v_c175_2.0.txt 2>&1
$PY W0.30 curve 500 prf v_c175_2.76 C175=2.76 > out/log_fit_v_c175_2.76.txt 2>&1
$PY W0.30 curve 0 prf v_cfree1.3 CFREE=1.3 > out/log_fit_v_cfree1.3.txt 2>&1
echo ALLDONE > out/variants_done.txt
