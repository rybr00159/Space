"""v4: what bound pair (Mamajek dwarf sequence, at the G star's distance 297.4 pc) is consistent with
 (i) the TESS colour levers (field-star curve f03), (ii) the eclipse depths (deep 0.0059, shallow 0.0027 of total TESS light;
     the light lost from star X cannot exceed its TESS flux times min(1, (R_other/R_X)^2)),
 (iii) the VHS feature contrast (dKs 2.9-3.3, dJ 3.7-4.3, if it is the pair), (iv) Gaia non-detection.
Prints, for every primary/secondary combination K5..M5, contrasts in G, T, J, Ks, Ic and the minimum geometric depth factor."""
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4lib import tmag, V4
rows = {}
for line in open(f'{V4}/ext/mamajek.txt'):
    if line.startswith('#') or len(line.split()) < 30:
        continue
    p = line.split()
    if p[0] in ('K5V', 'K7V', 'M0V', 'M0.5V', 'M1V', 'M1.5V', 'M2V', 'M2.5V', 'M3V', 'M3.5V', 'M4V', 'M4.5V', 'M5V'):
        def f(i):
            try: return float(p[i])
            except: return np.nan
        rows[p[0]] = dict(Teff=f(1), R=f(6), Mv=f(7), BpRp=f(11), MG=f(13), VIc=f(17), MJ=f(21), MKs=f(22), M=f(29))
DM = 5 * np.log10(1000 / 3.3622) - 5
G0 = dict(G=11.775, T=11.358, J=10.728, Ks=10.341, Ic=None)
# G star Ic from G and BP-RP: use V ~ 11.86 (APASS) and V-Ic ~0.69 for G2V
G0['Ic'] = 11.857 - 0.70
F3 = json.load(open(f'{V4}/results/f03_theta_red.json'))
cv = F3['curve']
cC = np.array([0.786] + [c['C'] for c in cv if c['C'] > 0.8]); cR = np.array([0.0] + [-c['rad'] for c in cv if c['C'] > 0.8])
def lever(c):
    if c <= cC[-1]:
        return float(np.interp(c, cC, cR))
    slope = (cR[-1] - cR[-2]) / (cC[-1] - cC[-2])
    return float(cR[-1] + slope * (c - cC[-1]))
def mags(s):
    r = rows[s]
    G = r['MG'] + DM; T = float(tmag(G, r['BpRp'])); J = r['MJ'] + DM; K = r['MKs'] + DM; I = r['Mv'] - r['VIc'] + DM
    return dict(G=G, T=T, J=J, Ks=K, Ic=I)
def comb(a, b):
    return -2.5 * np.log10(10 ** (-0.4 * a) + 10 ** (-0.4 * b))
print(f'DM = {DM:.3f}')
print(f"{'pair':12s} {'BpRp':>9s} {'lever_P':>7s} {'lever_S':>7s} {'dG':>5s} {'dT':>5s} {'dIc':>5s} {'dJ':>5s} {'dKs':>5s}  fT_P   fT_S  lostP/fT_P(<=(Rs/Rp)^2?)  lostS/fT_S")
out = []
for P in rows:
    for S in rows:
        if rows[S]['Teff'] > rows[P]['Teff']:
            continue
        a = mags(P); b = mags(S)
        d = {k: comb(a[k], b[k]) - G0[k] for k in a}
        fTP = 10 ** (-0.4 * (a['T'] - G0['T'])); fTS = 10 ** (-0.4 * (b['T'] - G0['T']))
        # total light fraction: relative to total system light (G star + pair)
        tot = 1 + fTP + fTS
        lostP = 0.0059 * tot / fTP; lostS = 0.0027 * tot / fTS
        geomP = (rows[S]['R'] / rows[P]['R']) ** 2
        ok = (lostP <= geomP) and (lostS <= 1.0)
        out.append(dict(P=P, S=S, dG=d['G'], dT=d['T'], dIc=d['Ic'], dJ=d['J'], dKs=d['Ks'], fT_pair=(fTP + fTS) / tot, lostP=lostP, lostS=lostS,
                        geomP=geomP, depth_ok=bool(ok), lever_P=lever(rows[P]['BpRp']), lever_S=lever(rows[S]['BpRp'])))
        print(f"{P + '+' + S:12s} {rows[P]['BpRp']:4.2f}/{rows[S]['BpRp']:4.2f} {lever(rows[P]['BpRp']):7.2f} {lever(rows[S]['BpRp']):7.2f} {d['G']:5.2f} {d['T']:5.2f} {d['Ic']:5.2f} {d['J']:5.2f} {d['Ks']:5.2f} "
              f"{fTP / tot:6.3f} {fTS / tot:6.3f}  {lostP:5.2f} vs {geomP:4.2f} {'ok' if lostP <= geomP else '--'}   {lostS:5.2f} {'ok' if lostS <= 1 else '--'}")
json.dump(out, open(f'{V4}/results/p01_pairs.json', 'w'), indent=1)
