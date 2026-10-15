#!/bin/bash
# combination variants (baseline injection noise unless stated)
cd /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/scripts
PY="/home/user/work/venv/bin/python -I"
c() { lab=$1; shift; env "$@" V4C_LABEL=$lab V4C_NS=${NS:-2000} OMP_NUM_THREADS=1 $PY s05_combine.py > ../results/log_s05_$lab.txt 2>&1; }
c deep_only V4C_TYPES=deep
c shallow_only V4C_TYPES=shallow
c kappa_equal V4C_KAPPA=equal
c kappa_theta V4C_KAPPA=theta
c kappa_red V4C_KAPPA=red
c no_bias V4C_NOBIAS=1
c err_x1.5 V4C_ERRX=1.5
for s in 1 28 68 95 105 106; do c drop_s$s V4C_DROP=$s; done
for v in box2 box4 plane quad ushape w2; do [ -f ../results/real_s106_$v.json ] && c var_$v V4C_REALTAG=_$v; done
