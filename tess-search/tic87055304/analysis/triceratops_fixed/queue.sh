#!/bin/bash
# Re-run TRICERATOPS for every candidate quoted in the README with the corrected depth input (run_tri.py fix).
cd /home/user/work/src
run() {  # name tic sectors P t0 dur lc envs...
  name=$1; tic=$2; secs=$3; P=$4; t0=$5; dur=$6; lc=$7; shift 7
  env "$@" TRI_RADIUS=10 nice -n 10 timeout 7200 ../venv_tri/bin/python run_tri.py $tic $secs $P $t0 $dur $lc /home/user/work/tri/fixed/tri_$name.json > /home/user/work/tri/fixed/tri_$name.log 2>&1
  echo "TRI_EXIT $? $name" >> /home/user/work/tri/fixed/queue.log
}
( run 87055304_tic 87055304 13,39,66,93,104 4.136190054 1658.411589 0.12492 /home/user/work/tri/lc_87055304_sap.npz TRI_MASS=1.3 TRI_PLX=2.807
  run 155991187 155991187 2,29,69,96,103,104 0.5012109557 1354.4850784 0.0206 /home/user/work/tri/lc_155991187.npz TRI_MASS=0.72
  run 158159741 158159741 13,66,93,94,104 1.8068919348 1658.3739087 0.03453 /home/user/work/tri/lc_158159741.npz ) &
( run 245976814 245976814 39,66,93,104 0.5070757303 2362.4568297 0.06167 /home/user/work/tri/lc_245976814.npz TRI_MASS=1.15
  run 72651793 72651793 3,30,97 0.6345290913 1386.2898695 0.0343 /home/user/work/tri/lc_72651793.npz
  run 320686138 320686138 13,27,67,94,104 11.914164772 1668.1597044 0.09026 /home/user/work/tri/lc_320686138_sap.npz ) &
wait
echo QUEUE_DONE >> /home/user/work/tri/fixed/queue.log
