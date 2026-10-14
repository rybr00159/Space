import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import x20_skymapper as X
from loclib import to_json
files = [f for f in sorted(glob.glob(f'{X.IMG}/sm_*.fits')) if '20160926' not in f]
with ProcessPoolExecutor(4) as ex:
    res = list(ex.map(X.run, files))
to_json([r for r in res if r], f'{X.R}/results/x20_skymapper.json')
