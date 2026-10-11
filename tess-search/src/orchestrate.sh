#!/bin/bash
# Run the remaining searches in sequence: wait for S105 -> S104 single-sector (overlap excluded) -> S104+S105 combined.
cd /home/user/work/src
while pgrep -f "^../venv/bin/python run_search.py ../data/s0105" > /dev/null; do sleep 30; done
../venv/bin/python run_search.py ../data/s0104 ../results/s0104 4 ../meta/badwin_s0104.json > ../results_s0104.log 2>&1
../venv/bin/python run_search_multi.py ../data/s0104 ../data/s0105 ../meta/badwin_s0104.json ../meta/badwin_s0105.json ../results/multi > ../results_multi.log 2>&1
echo ORCH_DONE >> ../results_multi.log
