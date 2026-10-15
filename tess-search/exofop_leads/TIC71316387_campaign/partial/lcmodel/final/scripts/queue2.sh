#!/bin/bash
# MS-model MCMC (2 workers) after the first injection batch
F=/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final
until grep -q '^time' $F/results/inject_A.log; do sleep 20; done
NW=48 /home/user/work/venv/bin/python -I $F/scripts/s08_mcmc.py ms 6000 ms 2 0.15 1.04 0.05 > $F/results/mcmc_ms.log 2>&1
