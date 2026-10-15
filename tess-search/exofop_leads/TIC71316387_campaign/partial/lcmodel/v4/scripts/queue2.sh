#!/bin/bash
V=/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4
PY=/home/user/work/venv/bin/python
until grep -qE '^time|Traceback' $V/results/inject_A.log; do sleep 20; done
NW=32 $PY -I $V/scripts/s08_mcmc.py geo 6000 geo 2 > $V/results/mcmc_geo.log 2>&1 &
NW=32 $PY -I $V/scripts/s08_mcmc.py ms 6000 ms 2 0.15 1.04 0.05 > $V/results/mcmc_ms.log 2>&1 &
NPROC=2 $PY -I $V/scripts/s03b_ms_ml_lite.py default 0.15 1.04 0.05 base > $V/results/msml_base.log 2>&1
NPROC=2 $PY -I $V/scripts/s14_trapezoid.py A > $V/results/trapezoid_A.log 2>&1
NPROC=2 $PY -I $V/scripts/s02b_variants.py > $V/results/variants.log 2>&1
NPROC=2 $PY -I $V/scripts/s03b_ms_ml_lite.py default 0.30 1.05 0.10 loose > $V/results/msml_loose.log 2>&1
wait
echo QUEUE_DONE > $V/results/queue2.done
