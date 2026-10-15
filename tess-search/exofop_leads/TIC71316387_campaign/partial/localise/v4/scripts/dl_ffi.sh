#!/bin/bash
cd /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/ffi
for b in tess2018216022942-s0001-1-2-0120-s tess2018228192941-s0001-1-2-0120-s tess2020218055915-s0028-1-1-0190-s tess2020232100914-s0028-1-1-0190-s tess2023216082832-s0068-1-1-0262-s tess2023230060152-s0068-1-1-0262-s tess2025213031924-s0095-1-2-0292-s tess2025226172923-s0095-1-2-0292-s tess2026172020230-s0105-2-3-0307-s tess2026185222910-s0105-2-3-0307-s tess2026199201909-s0106-1-4-0308-s tess2026216060909-s0106-1-4-0308-s; do
  f=${b}_ffic.fits
  [ -s $f ] || timeout 600 curl -sS -L -o $f "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/$f"
  echo "$f $(stat -c %s $f)"
done
