#!/bin/bash
# Rebuild all write-up figures (../figures/*.png). Inputs are read-only outputs of the analysis folders under
# ../../ (clean, referee3_clean, referee_loc, referee2_search, referee2_fpp, referee_star, referee2_shape, referee3_bd,
# referee3_followup) and the TESScut cubes in /home/user/work/tesscut5/214311906/. f0 must run before f1.
set -e
cd "$(dirname "$0")"
PY=/home/user/work/venv/bin/python
export OMP_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
$PY f0_rawap.py            # data/rawap_box3.npz (raw-pixel 3x3-box curves; checks against referee3_clean q02)
$PY f1_fold.py             # fold.png
$PY f2_per_transit.py      # per_transit.png
$PY f3_false_positives.py  # false_positives.png
$PY f4_host.py             # host.png
$PY f5_scenarios.py        # scenarios.png
$PY f6_localisation.py     # localisation.png
$PY f7_followup.py         # followup_2027.png
