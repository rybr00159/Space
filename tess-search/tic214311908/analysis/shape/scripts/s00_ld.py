"""s00: Claret (2017, A&A 600, A30) TESS quadratic limb-darkening coefficients from PHOENIX models (VizieR
J/A+A/600/A30/tableab), interpolated to the M3V (Teff 3350 +- 55 K, logg 4.84) and to the twin components
(3237 K, logg 4.94). Also gives the spread over Teff +-2 sigma, logg 4.5-5.0, both fitting methods (LSM, PCM) and
the PHOENIX 'r'-quasi-spherical vs plane-parallel options to set the prior width.
Output: out/ld_claret17.json, data/claret17_tableab.csv"""
import json, numpy as np
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = -1
t = Vizier.get_catalogs('J/A+A/600/A30/tableab')[0].to_pandas()
t.to_csv('data/claret17_tableab.csv', index=False)
print(t.columns.tolist()); print(t[['Teff','logg','Z','L/HP','Type','Mod']].drop_duplicates(['Type','Mod']).head(10))
print(t.Teff.min(), t.Teff.max(), sorted(t.logg.unique()), t.Z.unique(), t['Mod'].unique(), t['Type'].unique() if 'Type' in t else '')
from scipy.interpolate import griddata
res = {}
def interp(sub, teff, logg, col):
    pts = sub[['Teff','logg']].values.astype(float); v = sub[col].values.astype(float)
    return float(griddata(pts, v, (teff, logg), method='linear'))
for mod in ['PC','PD']:
    for typ in ['q','r']:
        sub = t[(t.Mod==mod)&(t.Type==typ)]
        if len(sub)==0: continue
        for meth in ['LSM','PCM']:
            for name,(te,lg) in {'single':(3350,4.84),'twin':(3237,4.94),'single_m2s':(3240,4.84),'single_p2s':(3460,4.84),'single_lg45':(3350,4.5),'single_lg50':(3350,5.0)}.items():
                a = interp(sub,te,lg,'a'+meth); b = interp(sub,te,lg,'b'+meth)
                if np.isfinite(a): res[f"{mod}_{typ}_{meth}_{name}"] = [a,b]
for k,v in res.items(): print(k, np.round(v,3))
A = np.array([v for k,v in res.items() if 'single' in k])
out = dict(source='Claret 2017 A&A 600 A30, VizieR J/A+A/600/A30/tableab (PHOENIX, TESS band, Z=0 only)', table=res,
           adopted_single=dict(u1=float(np.median([v[0] for k,v in res.items() if k.endswith('_single')])), u2=float(np.median([v[1] for k,v in res.items() if k.endswith('_single')]))),
           adopted_twin=dict(u1=float(np.median([v[0] for k,v in res.items() if k.endswith('_twin')])), u2=float(np.median([v[1] for k,v in res.items() if k.endswith('_twin')]))),
           spread_single=dict(u1_min=float(A[:,0].min()), u1_max=float(A[:,0].max()), u2_min=float(A[:,1].min()), u2_max=float(A[:,1].max())),
           prior_note='Gaussian prior with sigma 0.1 on each of u1, u2 around the adopted values (covers model/method spread and the absent [Fe/H]=+0.2 grid; PHOENIX vs observed M-dwarf LD differences are ~0.1, e.g. Patel & Espinoza 2022)')
print(out['adopted_single'], out['adopted_twin'], out['spread_single'])
json.dump(out, open('out/ld_claret17.json','w'), indent=1)
