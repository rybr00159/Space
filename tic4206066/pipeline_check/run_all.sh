#!/bin/sh
# Re-run the TIC 4206066 pipeline check. Needs numpy, numba, wotan, astropy, astroquery, scipy, pandas.
# Light curves download to $LC_CACHE (default ~/.cache/tess_lc) and from MAST directly.
cd "$(dirname "$0")"
PY=${PY:-python}
for s in would_detect.py would_detect2.py null_hires.py s6_b_check.py; do
  echo "### $s"
  $PY -W ignore "$s"
done
