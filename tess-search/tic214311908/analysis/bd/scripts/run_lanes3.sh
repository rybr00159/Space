#!/bin/bash
# Lane A: wait for the grid (PID $1) -> preliminary posterior (no injection calibration) -> compBD40 injections.
# Lane B: wait for the INV planet injections (PID $2) -> grazBD injections (REF, INV).
cd tic214311908/analysis/bd
PY=/home/user/work/venv/bin/python
( while [ -e /proc/$1 ]; do sleep 20; done
  $PY scripts/b09_posterior.py > logs/b09_prelim.log 2>&1
  $PY scripts/b06_inject.py ref compBD40 30 303 > logs/b06_ref_compBD40.log 2>&1
  $PY scripts/b06_inject.py inv compBD40 30 313 > logs/b06_inv_compBD40.log 2>&1
  echo done > logs/laneA.done ) &
( while [ -e /proc/$2 ]; do sleep 20; done
  $PY scripts/b06_inject.py ref grazBD 40 202 > logs/b06_ref_grazBD.log 2>&1
  $PY scripts/b06_inject.py inv grazBD 48 212 > logs/b06_inv_grazBD.log 2>&1
  echo done > logs/laneB.done ) &
wait
