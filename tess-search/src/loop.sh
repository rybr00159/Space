#!/bin/bash
# Re-run downstream steps every 20 min while the search orchestrator runs, then once more at the end.
cd /home/user/work/src
while true; do
  for s in 105 104; do
    if ls ../results/s0$s/chunk_*.parquet >/dev/null 2>&1; then
      ../venv/bin/python driver.py $s >> ../results/driver_s0$s.log 2>&1
      ../venv/bin/python eb_triage.py $s >> ../vet/eb_triage_s0$s.log 2>&1
    fi
  done
  pgrep -f "orchestrate.sh" > /dev/null || break
  sleep 1200
done
echo LOOP_DONE >> ../results/driver_s0105.log
