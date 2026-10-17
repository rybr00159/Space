import pandas as pd, numpy as np, json
d=pd.read_csv('out/p05_draws.csv'); v=pd.read_csv('out/p05_variants.csv')
d['FPP_half']=d.FPP+0.5*d.TP_twin
out={}
for name,s in [('all',d)]+[(k,g) for k,g in d.groupby('ds')]:
    out[name]=dict(n=len(s),FPP_5_50_95=[float(np.percentile(s.FPP,p)) for p in (5,50,95)],
                   frac_ge_0p015=float(np.mean(s.FPP>=0.015)),
                   FPPhalf_50_95=[float(np.percentile(s.FPP_half,p)) for p in (50,95)], frac_half_lt_0p015=float(np.mean(s.FPP_half<0.015)),
                   strict_50_95=[float(np.percentile(s.FPP_strict,p)) for p in (50,95)], frac_strict_lt_0p015=float(np.mean(s.FPP_strict<0.015)),
                   NFPP_95=float(np.percentile(s.NFPP,95)), frac_NFPP_ge_1em3=float(np.mean(s.NFPP>=1e-3)))
# strict by cmd
out['strict_median_by_cmd']=d.groupby('cmd').FPP_strict.median().to_dict()
out['FPP_median_by_cmd']=d.groupby('cmd').FPP.median().to_dict()
out['frac_strict_lt_0p015_colour_cmd']=float(np.mean(d[d.cmd=='colour'].FPP_strict<0.015))
out['FPP_ge_0p015_rows']=d[d.FPP>=0.015][['ds','ecc','cmd','ipd','f_mult','S_close','beta_N','FPP']].round(4).to_dict('records')
json.dump(out,open('out/p07_draw_stats.json','w'),indent=1)
print(json.dumps({k:v for k,v in out.items() if k!='FPP_ge_0p015_rows'},indent=1)); print(pd.DataFrame(out['FPP_ge_0p015_rows']).to_string())
print(v[['ds','variant','FPP','FPP_strict','NFPP']].sort_values('FPP').tail(14).to_string())
