#!/bin/bash
# sequential queue (2 workers) after the first injection batch: variants, loose-prior MS profile
F=/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final
until grep -q '^time' $F/results/inject_A.log; do sleep 20; done
NPROC=2 /home/user/work/venv/bin/python -I $F/scripts/s02b_variants.py > $F/results/variants.log 2>&1
NPROC=2 /home/user/work/venv/bin/python -I $F/scripts/s03_ms_ml.py default 0.30 1.05 0.10 loose > $F/results/msml_loose.log 2>&1
