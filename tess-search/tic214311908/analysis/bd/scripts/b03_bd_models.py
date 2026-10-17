"""b03: brown-dwarf radius / temperature prior from the ATMO 2020 evolutionary tracks (Phillips et al. 2020,
2020A&A...637A..38P; chemical-equilibrium set, downloaded 2026-10-17 from perso.ens-lyon.fr/isabelle.baraffe/ATMO2020,
data/models/evolutionary_tracks/ATMO_CEQ/MKO_WISE_IRAC/), plus the empirical radius inflation of transiting BDs.
Draws: M_BD from the observed M-dwarf short-period BD mass distribution (Larsen+25, 2025AJ....169..246L: 2/13 systems
with P < 13 d at 13-40 M_J): P(13-40) = 0.18, P(40-78.5) = 0.82, uniform inside; age uniform 1-10 Gyr (thin-disk
field star; star/ kinematics favour several Gyr, the X-ray activity a younger age - both inside this range); radius and
Teff interpolated in (mass, log age); radius inflation (1 + delta), delta ~ N(0.09, 0.10) (mean 5-14 %, intrinsic
dispersion 8-12 %: Beatty 2026, arXiv:2607.27352; 8.7 +- 1.9 %: Mehrotra+26, RNAAS 10, 128); R* ~ N(0.383, 0.012).
TESS-band flux ratio BD/M3V = (R_bd/R*)^2 x SB ratio, SB ratio = blackbody at Teff_BD over blackbody at 3350 K,
integrated over the TESS response (../shape/data/filt_TESS.Red.dat, photon-weighted). Checks: TOI-2119 b (T_b,TESS
2100 K for Teff 2030 K: 2022AJ....163...89C / 2022MNRAS.514.4944C) and the Mamajek dwarf table M_RP (L0V, L5V).
Outputs: data/bd_draws.npz, out/b03_bd_models.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import glob, pandas as pd
from scipy.interpolate import RegularGridInterpolator
rng = np.random.default_rng(3)
files = sorted(glob.glob(f'{BD}/data/models/evolutionary_tracks/ATMO_CEQ/MKO_WISE_IRAC/*_ATMO_CEQ_vega.txt'))
tr = {}
for fn in files:
    m = float(fn.split('/')[-1].split('_')[0])
    a = np.loadtxt(fn)
    tr[m] = a
masses = np.array(sorted(m for m in tr if 0.010 <= m <= 0.0751))
lages = np.linspace(np.log10(0.5), 1.0, 40)
RG = np.zeros((len(masses), len(lages))); TG = np.zeros_like(RG)
for i, m in enumerate(masses):
    a = tr[m]; la = np.log10(a[:, 1])
    RG[i] = np.interp(lages, la, a[:, 4]); TG[i] = np.interp(lages, la, a[:, 2])
fR = RegularGridInterpolator((masses, lages), RG); fT = RegularGridInterpolator((masses, lages), TG)
# TESS response
flt = np.loadtxt(f'{B}/shape/data/filt_TESS.Red.dat'); lam = flt[:, 0]*1e-8; resp = flt[:, 1]   # Angstrom -> cm
h, c, kB = 6.626e-27, 2.998e10, 1.381e-16
def band_sb(T):
    T = np.atleast_1d(T)[:, None]
    Bl = 2*h*c**2/lam**5/(np.exp(h*c/(lam*kB*T)) - 1)
    return np.trapezoid(Bl*resp*lam, lam, axis=1)     # photon-weighted
SB3350 = band_sb(3350.)[0]
N = 200000
MJ_per_Msun = 1047.57
u = rng.random(N)
Mbd = np.where(u < 0.18, rng.uniform(13, 40, N), rng.uniform(40, 78.5, N))
age = rng.uniform(1, 10, N)
mm = np.clip(Mbd/MJ_per_Msun, masses[0], masses[-1])
R0 = fR(np.column_stack([mm, np.log10(age)])); Teff = fT(np.column_stack([mm, np.log10(age)]))
delta = rng.normal(0.09, 0.10, N)
Rbd = R0*(1 + delta)                                   # Rsun
Rst = rng.normal(0.383, 0.012, N)
k = Rbd/Rst
fr = (Rbd/Rst)**2*band_sb(Teff)/SB3350
np.savez_compressed(f'{BD}/data/bd_draws.npz', M_MJ=Mbd.astype(np.float32), age=age.astype(np.float32), R_rsun=Rbd.astype(np.float32),
                    R_model=R0.astype(np.float32), Teff=Teff.astype(np.float32), k=k.astype(np.float32), Rstar=Rst.astype(np.float32),
                    flux_ratio=fr.astype(np.float32))
RJ = RJUP/RSUN
q = lambda x, p: float(np.percentile(x, p))
out = dict(R_RJ=dict(p2p5=q(Rbd/RJ, 2.5), p16=q(Rbd/RJ, 16), p50=q(Rbd/RJ, 50), p84=q(Rbd/RJ, 84), p97p5=q(Rbd/RJ, 97.5)),
           R_model_RJ=dict(p2p5=q(R0/RJ, 2.5), p50=q(R0/RJ, 50), p97p5=q(R0/RJ, 97.5)),
           k=dict(p2p5=q(k, 2.5), p16=q(k, 16), p50=q(k, 50), p84=q(k, 84), p97p5=q(k, 97.5), frac_below_0p19=float(np.mean(k < 0.19))),
           Teff=dict(p16=q(Teff, 16), p50=q(Teff, 50), p84=q(Teff, 84), p95=q(Teff, 95)),
           flux_ratio_ppt=dict(p50=q(fr*1e3, 50), p84=q(fr*1e3, 84), p95=q(fr*1e3, 95), p99=q(fr*1e3, 99), frac_gt_1ppt=float(np.mean(fr > 1e-3)), frac_gt_2ppt=float(np.mean(fr > 2e-3))))
grid = {}
for Mj in (15, 30, 50, 65, 75, 78.5):
    for ag in (1, 3, 5, 10):
        x = np.array([[min(Mj/MJ_per_Msun, masses[-1]), np.log10(ag)]])
        T_ = float(fT(x)[0]); R_ = float(fR(x)[0])
        grid[f'M{Mj}_age{ag}'] = dict(R_RJ=R_/RJ, Teff=T_, flux_ratio_ppt_noinfl=float((R_/0.383)**2*band_sb(T_)[0]/SB3350*1e3))
out['grid'] = grid
# checks
out['check_TOI2119'] = dict(note='TOI-2119 b: R 1.08 R_J around R* 0.50; T_b 2100 K -> predicted TESS secondary (BB, T_host 3621 K)',
                            pred_ppt=float((1.08*RJ/0.50)**2*band_sb(2100.)[0]/band_sb(3621.)[0]*1e3))
eem = open(f'{B}/referee_star/data/EEM_dwarf.txt').read().splitlines()
def row(spt):
    for l in eem:
        if l.startswith(spt + ' '): return l.split()
hdr = [l for l in eem if l.startswith('#SpT')][0].lstrip('#').split()
def val(r, col): return float(r[hdr.index(col)])
m3 = row('M3V'); chk = {}
for sp in ('L0V', 'L5V'):
    r = row(sp)
    MRP_bd = val(r, 'M_G') - val(r, 'G-Rp'); MRP_m3 = val(m3, 'M_G') - val(m3, 'G-Rp')
    ratio_tab = 10**(-0.4*(MRP_bd - MRP_m3))
    ratio_bb = (val(r, 'R_Rsun')/val(m3, 'R_Rsun'))**2*band_sb(val(r, 'Teff'))[0]/band_sb(val(m3, 'Teff'))[0]
    chk[sp] = dict(Teff=val(r, 'Teff'), ratio_RP_table=float(ratio_tab), ratio_BB_TESS=float(ratio_bb), BB_over_table=float(ratio_bb/ratio_tab))
out['check_mamajek'] = chk
print(json.dumps(out, indent=1))
jdump(dict(method=__doc__, **out), f'{BD}/out/b03_bd_models.json')
