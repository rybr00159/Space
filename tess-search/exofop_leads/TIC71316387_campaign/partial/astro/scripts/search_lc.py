import lightkurve as lk, warnings; warnings.filterwarnings('ignore')
sr=lk.search_lightcurve('TIC 71316387')
print(sr)
