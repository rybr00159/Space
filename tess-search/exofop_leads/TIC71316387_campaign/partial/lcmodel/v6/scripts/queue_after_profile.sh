#!/bin/bash
# wait for the real-data L3 profile to finish, then run the geometric MCMC and the robustness variants
B=/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6
cd $B
while pgrep -f "[c]04_profile.py default results/seeds_default.json default SPOC1,SPOC28,TC95,TC105,TC106 L3" >/dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 30000 /home/user/work/venv/bin/python -I scripts/c03_mcmc.py default results/best_default.json geo 10000 56 2 3.0 SPOC1,SPOC28,TC95,TC105,TC106 > logs/c03_geo.log 2>&1 &
VFILE=results/variants_a.json timeout 30000 /home/user/work/venv/bin/python -I scripts/c08_variants.py e0 ldclaret noscale tcall qlp > logs/c08_a.log 2>&1
VFILE=results/variants_b.json timeout 30000 /home/user/work/venv/bin/python -I scripts/c08_variants.py lin early late esin0 > logs/c08_b.log 2>&1
wait
