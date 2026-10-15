#!/bin/bash
# analysis variants of the real per-sector fits (injection noise of the baseline is reused in the combination)
cd /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/scripts
PY="/home/user/work/venv/bin/python -I"
run() { tag=$1; shift; env "$@" V4_TAG=_$tag OMP_NUM_THREADS=1 $PY s02_real.py > ../results/log_s02_$tag.txt 2>&1; }
if [ "$1" = "a" ]; then
  run box2 V4_BOX=2; run box4 V4_BOX=4; run plane V4_PLANE=1
else
  run quad V4_BASEDEG=2; run ushape V4_SHAPE=U; run w2 V4_W=2.0
fi
