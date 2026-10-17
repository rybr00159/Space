#!/bin/bash
# lane B: wait for s02, then the evidence grid
cd tic214311908/analysis/referee3_bd
PY=/home/user/work/venv/bin/python
while pgrep -f "s02_fpp_fixed.py" > /dev/null; do sleep 10; done
$PY scripts/s06_evidence.py > logs/s06.log 2>&1
$PY scripts/s07_ptp_followup.py > logs/s07.log 2>&1
