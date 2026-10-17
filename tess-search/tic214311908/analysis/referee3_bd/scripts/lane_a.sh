#!/bin/bash
# lane A: wait for the running REF grazBD injection, then the other injection truths
cd tic214311908/analysis/referee3_bd
PY=/home/user/work/venv/bin/python
while pgrep -f "s04_inject.py ref grazBD" > /dev/null; do sleep 10; done
$PY scripts/s04_inject.py inv grazBD 28 202 > logs/s04_inv_grazBD.log 2>&1
$PY scripts/s04_inject.py ref planet 50 303 > logs/s04_ref_planet.log 2>&1
$PY scripts/s04_inject.py inv planet 28 404 > logs/s04_inv_planet.log 2>&1
