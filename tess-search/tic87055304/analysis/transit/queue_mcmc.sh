#!/bin/bash
# Sequential queue of MCMC runs (one heavy process at a time). Log per run in log_mcmc_<name>.txt
cd /home/user/work/tic87/transit
export OMP_NUM_THREADS=1
PY=/home/user/work/venv/bin/python
run() { name=$1; shift; nice -n 10 $PY run_mcmc.py $name "$@" > log_mcmc_$name.txt 2>&1; }
run subgiant --rho uniform,0.04,0.10 --maxsteps 30000
run dil051 --rho uniform,0.003,5.0 --dil 0.51 --maxsteps 30000
for s in 13 39 66 93 104; do run drop$s --drop $s --maxsteps 24000; done
run oddeven --oddeven --maxsteps 30000
echo QUEUE_DONE > queue_done.txt
