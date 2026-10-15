"""Main-sequence relations from the Pecaut & Mamajek dwarf sequence (EEM_dwarf_UBVIJHK_colors_Teff, 2022.04.16):
mass -> Teff, R, M_G, Bp-Rp, M_Ks; TESS absolute magnitude from G and Bp-Rp with the TIC-8 relation
(Stassun et al. 2019); TESS-band surface brightness SB_T = 10^(-0.4 M_T) / R^2 (arbitrary units).
Interpolation is linear in log(mass) over K0V-M6V (0.10-0.88 Msun)."""
import numpy as np
F = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/ext/mamajek.txt'

def _load():
    rows = []
    started = False
    for line in open(F):
        if line.startswith('#SpT'):
            if started: break
            started = True; hdr = line[1:].split(); continue
        if not started or line.startswith('#') or not line.strip(): continue
        v = line.split()
        if len(v) < len(hdr): continue
        rows.append(v[:len(hdr)])
    out = {}
    for j, h in enumerate(hdr):
        col = []
        for r in rows:
            try: col.append(float(r[j]))
            except ValueError: col.append(np.nan)
        out[h] = np.array(col)
    out['SpT'] = np.array([r[0] for r in rows])
    return out

T = _load()
sel = np.isfinite(T['Msun']) & np.isfinite(T['M_G']) & np.isfinite(T['Bp-Rp']) & np.isfinite(T['R_Rsun']) & (T['Msun'] <= 0.9) & (T['Msun'] >= 0.09)
SPT = T['SpT'][sel]
_m = T['Msun'][sel]; o = np.argsort(_m)
SPT = SPT[o]; _m = _m[o]
_lm = np.log(_m)
_teff = T['Teff'][sel][o]; _R = T['R_Rsun'][sel][o]; _MG = T['M_G'][sel][o]; _bprp = T['Bp-Rp'][sel][o]; _MK = T['M_Ks'][sel][o]
_c = _bprp
_MT = _MG - 0.00522555 * _c ** 3 + 0.0891337 * _c ** 2 - 0.633923 * _c + 0.0324473
MMIN, MMAX = float(_m.min()), float(_m.max())

def interp(m, arr):
    return np.interp(np.log(np.asarray(m, float)), _lm, arr)
def teff(m): return interp(m, _teff)
def radius(m): return np.exp(interp(m, np.log(_R)))
def MT(m): return interp(m, _MT)
def MKs(m): return interp(m, _MK)
def MGmag(m): return interp(m, _MG)
def SB_T(m): return 10 ** (-0.4 * MT(m)) / radius(m) ** 2
def spt(m):
    return SPT[np.argmin(np.abs(_lm - np.log(m)))]
def density(m):   # solar units
    return m / radius(m) ** 3
_mg = np.exp(np.linspace(np.log(MMIN), np.log(MMAX), 3000)); _rho = density(_mg)
def mass_from_density(rho):
    # density decreases monotonically with mass on the lower MS
    o = np.argsort(_rho)
    return np.interp(rho, _rho[o], _mg[o], left=np.nan, right=np.nan)
_sb = SB_T(_mg)
def m2_from_J(J, m1):
    """secondary mass whose TESS surface brightness ratio to a primary of mass m1 equals J"""
    o = np.argsort(_sb)
    return np.interp(np.asarray(J) * SB_T(m1), _sb[o], _mg[o], left=np.nan, right=np.nan)

if __name__ == '__main__':
    for m in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
        print(m, spt(m), round(float(teff(m))), round(float(radius(m)), 3), round(float(MT(m)), 2), round(float(density(m)), 2))
