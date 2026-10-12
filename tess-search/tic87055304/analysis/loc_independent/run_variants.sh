#!/bin/bash
# Robustness variants of the real difference-image analysis (two sequential queues = 2 processes max).
cd /home/user/work/tic87/loc_independent
PY=/home/user/work/venv/bin/python
export OMP_NUM_THREADS=1
q1() {
  nice -n 10 $PY s02_diffimg.py prf 2.0 local_r4 prf_r2 > log_var_prf_r2.txt 2>&1
  nice -n 10 $PY s02_diffimg.py prf 2.5 local_r4 prf_r2.5 > log_var_prf_r2.5.txt 2>&1
  nice -n 10 $PY s02_diffimg.py prf 4.0 local_r4 prf_r4 > log_var_prf_r4.txt 2>&1
  nice -n 10 $PY s02_diffimg.py prf 3 external_2.5-8 prf_regext > log_var_prf_regext.txt 2>&1
  nice -n 10 $PY s02_diffimg.py prf 3 local_r4_ratiofix prf_regratio > log_var_prf_regratio.txt 2>&1
}
q2() {
  nice -n 10 $PY s02_diffimg.py prf 3 local_r2.5 prf_regr2.5 > log_var_prf_regr2.5.txt 2>&1
  LOC_BASELINE=const nice -n 10 $PY s02_diffimg.py prf 3 local_r4 prf_constbase > log_var_prf_constbase.txt 2>&1
  LOC_OOT_IN=1.5 LOC_OOT_OUT=3.0 nice -n 10 $PY s02_diffimg.py prf 3 local_r4 prf_ootfar > log_var_prf_ootfar.txt 2>&1
  LOC_IN_HALF=0.25 nice -n 10 $PY s02_diffimg.py prf 3 local_r4 prf_in50 > log_var_prf_in50.txt 2>&1
  nice -n 10 $PY s02_diffimg.py gauss 3 local_r4_ratiofix gauss_pointing pointing > log_var_gauss_pointing.txt 2>&1
}
q1 &
q2 &
wait
echo ALLDONE > log_variants_done.txt
