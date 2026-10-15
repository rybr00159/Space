#!/bin/bash
# after the injection batch: MS-model ML profile (4 workers), then geo and MS MCMC chains (2 workers each), then variants
V=/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4
PY=/home/user/work/venv/bin/python
until grep -qE '^time|Traceback' $V/results/inject_A.log; do sleep 20; done
NPROC=4 $PY -I $V/scripts/s03_ms_ml.py default 0.15 1.04 0.05 base > $V/results/msml_base.log 2>&1
NW=40 $PY -I $V/scripts/s08_mcmc.py geo 8000 geo 2 > $V/results/mcmc_geo.log 2>&1 &
NW=40 $PY -I $V/scripts/s08_mcmc.py ms 8000 ms 2 0.15 1.04 0.05 > $V/results/mcmc_ms.log 2>&1 &
wait
NPROC=4 $PY -I $V/scripts/s02b_variants.py > $V/results/variants.log 2>&1
NPROC=4 $PY -I $V/scripts/s03_ms_ml.py default 0.30 1.05 0.10 loose > $V/results/msml_loose.log 2>&1
echo QUEUE_DONE > $V/results/queue_main.done
