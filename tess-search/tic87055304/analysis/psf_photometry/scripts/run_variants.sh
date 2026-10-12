#!/bin/bash
# Robustness variants of the PSF photometry (run sequentially, niced, single-threaded)
export OMP_NUM_THREADS=1
cd /home/user/work/tic87/psf_photometry/scripts
PY=/home/user/work/venv/bin/python
run() { tag=$1; shift; mkdir -p ../lc/$tag; nice -n 10 $PY psfphot.py --tag $tag "$@" > ../lc/$tag/log.txt 2>&1; nice -n 10 $PY analyze_lc.py --tag $tag >> ../lc/$tag/log.txt 2>&1; tail -1 ../lc/$tag/log.txt; }
run noreg --mode none
run pairreg --pairreg
run off_p05 --offset_nb 0.05
run off_m05 --offset_nb -0.05
run off_p10 --offset_nb 0.10
run off_m10 --offset_nb -0.10
run rfit2p5 --rfit 2.5
run rfit5 --rfit 5.0
run tfree14 --tfree 14.0
run tfree16p5 --tfree 16.5
run smear025 --sig 0.25
run smear050 --sig 0.50
run gauss --prf gauss
run prf_corner --ccdpos 300,300
run freereg --mode free
