#!/bin/bash
# Injection runs in two lanes (at most 2 heavy processes). Lane A waits for PID $1 (ref/planet job), lane B for PID $2 (b05b).
cd tic214311908/analysis/bd
PY=/home/user/work/venv/bin/python
( while [ -e /proc/$1 ]; do sleep 20; done
  $PY scripts/b06_inject.py ref grazBD 40 202 > logs/b06_ref_grazBD.log 2>&1
  $PY scripts/b06_inject.py ref compBD40 30 303 > logs/b06_ref_compBD40.log 2>&1
  echo done > logs/laneA.done ) &
( while [ -e /proc/$2 ]; do sleep 20; done
  $PY scripts/b06_inject.py inv planet 48 111 > logs/b06_inv_planet.log 2>&1
  $PY scripts/b06_inject.py inv grazBD 48 212 > logs/b06_inv_grazBD.log 2>&1
  $PY scripts/b06_inject.py inv compBD40 30 313 > logs/b06_inv_compBD40.log 2>&1
  echo done > logs/laneB.done ) &
wait
