cd /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v3
P=/home/user/work/venv/bin/python
export V3_NS=1500
V3_TYPES=deep V3_LABEL=deep_only $P scripts/s05_combine.py
V3_TYPES=shallow V3_LABEL=shallow_only $P scripts/s05_combine.py
V3_KMODE=equal V3_LABEL=k_equal $P scripts/s05_combine.py
V3_KMODE=theta V3_LABEL=k_theta $P scripts/s05_combine.py
V3_NOBIAS=1 V3_LABEL=no_bias $P scripts/s05_combine.py
V3_ERRSCALE=1.5 V3_LABEL=err_x1.5 $P scripts/s05_combine.py
for s in 1 28 68 95 105 106; do V3_DROP=$s V3_LABEL=drop_s$s $P scripts/s05_combine.py; done
