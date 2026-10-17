import pandas as pd, sys
B='tic214311908/analysis/'
d=pd.read_csv(B+'referee3_bd/out/s02_fpp_fixed.csv')
pd.set_option('display.width',300); pd.set_option('display.max_rows',300)
d['f_BD']=d['f_BD'].astype(str)
d['ecc']=d['ecc'].astype(str)
cols=['ds','ecc','f_BD','sign','rvf','FPP','FPP_system','PTP','BD_total','NTP_pairA','WD','HEB','BG','TP_single','TP_diluted','TP_twin','P_planet_any_star']
f=d[d.sign=='fixed']
print(f[cols].to_string(float_format=lambda x:'%.3g'%x))
