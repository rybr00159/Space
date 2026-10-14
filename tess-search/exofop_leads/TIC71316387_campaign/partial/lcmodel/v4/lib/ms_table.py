"""Main-sequence relations from Pecaut & Mamajek (2013, version in ext/mamajek.txt): Teff, R, M, M_G, Bp-Rp, M_Ks.
TESS magnitude from Gaia G and Bp-Rp via Stassun et al. (2019) eq. 1.  Surface brightness in TESS band
S_T = 10^(-0.4 M_T)/R^2 (disc-averaged, since M_T is disc-integrated).  Interpolated in log mass."""
import numpy as np
HERE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel'
def _load():
    rows = []
    with open(f'{HERE}/ext/mamajek.txt') as fh:
        lines = fh.readlines()
    i0 = [i for i, l in enumerate(lines) if l.startswith('#SpT')][0]
    hdr = lines[i0][1:].split()
    for l in lines[i0 + 1:]:
        if l.startswith('#'): break
        v = l.split()
        if len(v) < len(hdr) - 1: continue
        rows.append(v)
    def col(name):
        j = hdr.index(name)
        out = []
        for r in rows:
            try: out.append(float(r[j]))
            except: out.append(np.nan)
        return np.array(out)
    spt = np.array([r[0] for r in rows])
    d = dict(spt=spt, teff=col('Teff'), R=col('R_Rsun'), M=col('Msun'), MG=col('M_G'), bprp=col('Bp-Rp'), MKs=col('M_Ks'), MJ=col('M_J'), Mv=col('Mv'))
    return d
T = _load()
def tmag(G, c):
    return G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
T['MT'] = tmag(T['MG'], T['bprp'])
ok = np.isfinite(T['M']) & np.isfinite(T['R']) & np.isfinite(T['MT']) & np.isfinite(T['teff']) & (T['M'] <= 1.6) & (T['M'] >= 0.09)
for k in list(T): T[k] = T[k][ok]
o = np.argsort(T['M'])
for k in list(T): T[k] = T[k][o]
# remove duplicate masses (keep first)
_, ui = np.unique(T['M'], return_index=True)
for k in list(T): T[k] = T[k][ui]
lm = np.log10(T['M'])
def interp(name, M):
    return np.interp(np.log10(M), lm, T[name])
def radius(M): return 10 ** np.interp(np.log10(M), lm, np.log10(T['R']))
def teff(M): return 10 ** np.interp(np.log10(M), lm, np.log10(T['teff']))
def MT(M): return np.interp(np.log10(M), lm, T['MT'])
def MKs(M): return np.interp(np.log10(M), lm, T['MKs'])
def SB_T(M):
    """TESS-band disc-averaged surface brightness (arbitrary units)"""
    return 10 ** (-0.4 * MT(M)) / radius(M) ** 2
def spt_of_teff(t):
    i = np.argmin(np.abs(T['teff'] - t)); return T['spt'][i]
def spt_of_mass(m):
    i = np.argmin(np.abs(T['M'] - m)); return T['spt'][i]
if __name__ == '__main__':
    for i in range(len(T['M'])):
        print(f"{T['spt'][i]:6s} M={T['M'][i]:.3f} R={T['R'][i]:.3f} Teff={T['teff'][i]:.0f} MG={T['MG'][i]:.2f} BpRp={T['bprp'][i]:.2f} MT={T['MT'][i]:.2f} MKs={T['MKs'][i]:.2f} rho={T['M'][i]/T['R'][i]**3:.2f}")
