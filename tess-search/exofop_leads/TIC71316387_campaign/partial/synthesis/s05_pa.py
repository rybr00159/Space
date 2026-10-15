import numpy as np, json
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
z=np.load(B+'localise/v3/results/s05_baseline_map.npz'); g=z['g']; M=z['M']
from scipy.interpolate import RegularGridInterpolator
I=RegularGridInterpolator((g,g),M)
out={}
pa=np.arange(0,360,1.0)
for s in (0.4,0.45,0.5,0.6,0.7,0.8):
    dE=s*np.sin(np.radians(pa)); dN=s*np.cos(np.radians(pa))
    v=I(np.c_[dN,dE])
    def rng_(thr):
        ok=pa[v<thr]
        if len(ok)==0: return None
        # express as circular range
        okc=np.where(ok>180,ok-360,ok)
        return [float(okc.min()%360), float(okc.max()%360), int(len(ok))]
    out[str(s)]=dict(min_dchi2=float(v.min()), pa_at_min=float(pa[v.argmin()]), dchi2_at_opposite=float(v[(int(pa[v.argmin()])+180)%360]), pa_range_68=rng_(2.28), pa_range_95=rng_(6.54))
dd=json.load(open(B+'localise/v3/results/s05_drop_s105.json')) if __import__('os').path.exists(B+'localise/v3/results/s05_drop_s105.json') else None
if dd: out['drop_s105']=dict(delta=dd.get('M2_delta'), sd=dd.get('M2_delta_emp_sd'))
print(json.dumps(out,indent=1)); json.dump(out,open(B+'synthesis/s05_pa.json','w'),indent=1)
