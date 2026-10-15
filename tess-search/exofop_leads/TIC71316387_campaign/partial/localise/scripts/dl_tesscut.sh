#!/bin/bash
D=/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/data/tesscut
for s in 1 28 68 95 105 106; do
  f=$D/tess_s${s}_21x21.zip
  timeout 1500 curl -sS -L -o $f "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=330.887948&dec=-42.835511&y=21&x=21&sector=${s}"
  echo "done $s $(ls -la $f)"
  mkdir -p $D/s$s; unzip -o -q $f -d $D/s$s && rm -f $f
done
ls -laR $D
