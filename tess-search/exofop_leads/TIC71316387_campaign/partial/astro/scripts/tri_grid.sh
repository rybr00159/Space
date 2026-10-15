#!/bin/bash
# runs TRICERATOPS variant grid with 3 concurrent jobs
cd /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/tri
export TRI_MASS=1.05 TRI_RAD=1.07 TRI_PLX=3.3622 TRI_DROP_NEARBY=1 TRI_DROP=TP,PTP,STP,DTP,BTP
PY=/home/user/work/venv_tri/bin/python; R=../scripts/run_tri_astro.py
jobs_list=()
for s in 11 12 13 14 15 16 17 18; do jobs_list+=("$s 28,68 0.6518369 pdc_s28s68_fold.npz gA_s28s68_seed$s"); done
for s in 21 22 23 24; do jobs_list+=("$s 68 0.6518369 pdc_s68_fold.npz gB_s68_seed$s"); done
for s in 31 32 33 34; do jobs_list+=("$s 95,105 0.6518369 qlp_s95s105_fold.npz gC_qlp_seed$s"); done
for s in 41 42 43 44; do jobs_list+=("$s 28,68 1.3036712 pdc_s28s68_deeponly.npz gD_deepPEB_seed$s"); done
for j in "${jobs_list[@]}"; do
  set -- $j
  while [ $(jobs -r | wc -l) -ge 3 ]; do sleep 5; done
  TRI_SEED=$1 timeout 3600 $PY $R 71316387 $2 $3 0 0.064 $4 $5.json > $5.log 2>&1 &
done
wait
echo ALLDONE
