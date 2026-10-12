#!/bin/bash
cd /home/user/work/src
run() { name=$1; tic=$2; secs=$3; P=$4; t0=$5; dur=$6; lc=$7; shift 7
  env "$@" TRI_RADIUS=10 nice -n 10 timeout 7200 ../venv_tri/bin/python run_tri.py $tic $secs $P $t0 $dur $lc /home/user/work/tri/fixed/tri_$name.json > /home/user/work/tri/fixed/tri_$name.log 2>&1
  echo "TRI_EXIT $? $name" >> /home/user/work/tri/fixed/queue.log; }
run 245976814 245976814 39,66,93,104 0.5070757303 2362.4568297 0.06167 /home/user/work/tri/lc_245976814.npz TRI_MASS=1.15 TRI_PLX=2.87 &
run 155991187 155991187 2,29,69,96,103,104 0.5012109557 1354.4850784 0.0206 /home/user/work/tri/lc_155991187.npz TRI_MASS=0.68 TRI_RAD=0.69 TRI_PLX=6.58 &
wait; echo QUEUE2_DONE >> /home/user/work/tri/fixed/queue.log
