"""Smooth main-sequence relations (0.09-1.25 Msun) from least-squares polynomial fits in x = log10(M/Msun) to the
Pecaut & Mamajek (2013; table version 2022.04.16) dwarf sequence.  The raw table is lumpy at the M2.5V/M3V step
(R jumps 0.36 -> 0.42 Rsun between 0.37 and 0.40 Msun), which would imprint artificial edges on a mass posterior; the
smooth fits keep the trend and the scatter about them is reported (and absorbed by the fR / dlnJ nuisance priors).
TESS magnitude from Gaia G and Bp-Rp (Stassun et al. 2019 eq. 1).  TESS-band disc-averaged surface brightness
SB_T = 10^(-0.4 M_T) / R^2."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import ms_table as MT
T = MT.T
sel = (T['M'] >= 0.09) & (T['M'] <= 1.25)
x = np.log10(T['M'][sel])
DEG = 4
FITS = {}
SCAT = {}
for name, y in [('logR', np.log10(T['R'][sel])), ('logT', np.log10(T['teff'][sel])), ('MT', T['MT'][sel]), ('MKs', T['MKs'][sel]),
                ('MG', T['MG'][sel])]:
    ok = np.isfinite(y)
    c = np.polyfit(x[ok], y[ok], DEG)
    FITS[name] = c
    SCAT[name] = float(np.std(y[ok] - np.polyval(c, x[ok])))
def radius(M): return 10 ** np.polyval(FITS['logR'], np.log10(M))
def teff(M): return 10 ** np.polyval(FITS['logT'], np.log10(M))
def MTess(M): return np.polyval(FITS['MT'], np.log10(M))
def MKs(M): return np.polyval(FITS['MKs'], np.log10(M))
def MG(M): return np.polyval(FITS['MG'], np.log10(M))
def SB_T(M): return 10 ** (-0.4 * MTess(M)) / radius(M) ** 2
def spt_of_mass(m):
    i = np.argmin(np.abs(T['M'] - m)); return str(T['spt'][i])
def mass_of_rho(rho_sun, fR=1.0):
    """invert rho = M / (fR R(M))^3 on 0.09-1.25 Msun (monotonic)"""
    mg = np.exp(np.linspace(np.log(0.09), np.log(1.25), 2000)); rg = mg / (fR * radius(mg)) ** 3
    return float(np.interp(rho_sun, rg[::-1], mg[::-1]))
if __name__ == '__main__':
    print('scatter about smooth fits:', {k: round(v, 4) for k, v in SCAT.items()})
    for m in [0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.03]:
        print(f'M={m:.2f} {spt_of_mass(m):6s} R={radius(m):.3f} Teff={teff(m):.0f} MT={MTess(m):.2f} MKs={MKs(m):.2f} MG={MG(m):.2f} rho={m/radius(m)**3:.2f}')
