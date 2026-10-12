#!/bin/bash
# TRICERATOPS runs with TIC 87055315 as the target, using the fixed runner /home/user/work/src/run_tri.py (unmodified).
# One TRICERATOPS process at a time (house rule). LC: ../data/tri_lc_nb_200s.npz (undiluted neighbour LC, S66/S93/S104).
cd /home/user/work/tic87/neighbour/tri
LC=/home/user/work/tic87/neighbour/data/tri_lc_nb_200s.npz
P=4.1362006; T0=3846.46449; DUR=0.17593
run() {  # name envs...
  name=$1; shift
  echo "START $name $(date -u +%FT%T)" >> queue.log
  env "$@" TRI_RADIUS=10 OMP_NUM_THREADS=1 nice -n 10 timeout 5400 /home/user/work/venv_tri/bin/python /home/user/work/src/run_tri.py \
      87055315 66,93,104 $P $T0 $DUR $LC /home/user/work/tic87/neighbour/tri/tri_$name.json > tri_$name.log 2>&1
  echo "EXIT $? $name $(date -u +%FT%T)" >> queue.log
}
run A_ticparams TRI_DEPTH=0.0212
run B_adopted TRI_DEPTH=0.0212 TRI_MASS=1.17 TRI_RAD=1.93 TRI_PLX=0.7785
run C_adopted_repeat TRI_DEPTH=0.0212 TRI_MASS=1.17 TRI_RAD=1.93 TRI_PLX=0.7785
echo QUEUE_DONE >> queue.log
