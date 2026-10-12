#!/bin/bash
# Waits for run_tri_queue.sh to finish (one TRICERATOPS process at a time), then runs the on-target variant
# (scripts/run_tri_ontarget.py = copy of the fixed runner + remove_nearby_stars()).
cd /home/user/work/tic87/neighbour/tri
until grep -q QUEUE_DONE queue.log; do sleep 30; done
LC=/home/user/work/tic87/neighbour/data/tri_lc_nb_200s.npz
P=4.1362006; T0=3846.46449; DUR=0.17593
name=D_ontarget_adopted
echo "START $name $(date -u +%FT%T)" >> queue.log
env TRI_ONTARGET=1 TRI_DEPTH=0.0212 TRI_MASS=1.17 TRI_RAD=1.93 TRI_PLX=0.7785 TRI_RADIUS=10 OMP_NUM_THREADS=1 nice -n 10 timeout 5400 \
  /home/user/work/venv_tri/bin/python /home/user/work/tic87/neighbour/scripts/run_tri_ontarget.py \
  87055315 66,93,104 $P $T0 $DUR $LC /home/user/work/tic87/neighbour/tri/tri_$name.json > tri_$name.log 2>&1
echo "EXIT $? $name $(date -u +%FT%T)" >> queue.log
echo QUEUE2_DONE >> queue.log
