"""b10: quantify the photo-eccentric sign error shared by ../fpp/scripts/f04_bayes.py and ../referee2_fpp/scripts/r05_fpp.py.
Both evaluate the light-curve grid at rho_circ = rho_host * g^-3 with g = (1 + e sin w)/sqrt(1 - e^2) while weighting by the
transit probability factor (1 + e sin w)/(1 - e^2). The correct relation (Dawson & Johnson 2012; verified numerically with
batman in this script) is rho_circ = rho_host * g^+3: a transit near periastron (sin w > 0, boosted transit probability)
is SHORTER, i.e. it looks like a DENSER star. The error therefore attaches the periastron probability boost to long
(low-rho_circ) durations - exactly what a grazing eccentric BD on the M3V needs.
This script (i) verifies the relation with batman; (ii) imports referee2_fpp/r05_fpp.py read-only and re-evaluates its
baseline and the ecc/noise variants with the eccentricity function monkey-patched to return 1/g (so that its
rho_c = rho_host*(1/g)^-3 = rho_host*g^3), everything else unchanged.
Output: out/b10_r05_eccsign.json"""
import os
for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'): os.environ[k] = '1'
import sys, json, numpy as np, batman
B = 'tic214311908/analysis'
OUT = f'{B}/bd/out'
# ---- (i) batman check
t = np.linspace(-0.25, 0.25, 50001)
def T14(e, w):
    p = batman.TransitParams(); p.t0 = 0; p.per = 5.6915; p.rp = 0.12; p.a = 25.4; p.inc = 90.; p.ecc = e; p.w = w
    p.u = [0.2, 0.4]; p.limb_dark = 'quadratic'
    f = batman.TransitModel(p, t).light_curve(p); ii = np.where(f < 1)[0]; return (t[ii[-1]] - t[ii[0]])*24
Tc = T14(0., 90.)
chk = []
for e, w in [(0.5, 90.), (0.5, 270.), (0.3, 30.), (0.6, 250.)]:
    g = (1 + e*np.sin(np.radians(w)))/np.sqrt(1 - e**2); T = T14(e, w)
    chk.append(dict(e=e, w_deg=w, T14_h=T, T14_circ_h=Tc, rho_circ_over_rho=(Tc/T)**3, g3=g**3, g_minus3=g**-3,
                    transit_prob_factor=(1 + e*np.sin(np.radians(w)))/(1 - e**2)))
    print(chk[-1])
# ---- (ii) re-evaluate r05 with the corrected sign
sys.path.insert(0, f'{B}/referee2_fpp/scripts')
import r05_fpp as R
orig_ecc = R.ecc
def ecc_fixed(kind, n=2500):
    ge, geo = orig_ecc(kind, n)
    return 1/ge, geo
res = {}
def run(tag, patched, **kw):
    R.ZT.clear(); R.ZCACHE.clear()
    R.ecc = ecc_fixed if patched else orig_ecc
    R.rng = np.random.default_rng(11)
    p = dict(R.BASE); p.update(kw)
    r = R.evaluate(p, detail=True)
    po = r['post']
    bd = po['BD_on_M3V_single'] + po['BD_on_primary_of_binary'] + po['BDHEB_on_companion']
    res[tag] = dict(patched=patched, params=kw, FPP=r['FPP'], FPP_strict=r['FPP_strict'], FPP_system=r['FPP_system'], NFPP=r['NFPP'],
                    PTP=po['PTP_companion'], BD_on_M3V=po['BD_on_M3V_single'], BD_on_primary=po['BD_on_primary_of_binary'],
                    BD_on_companion=po['BDHEB_on_companion'], BD_total=bd, TP_single=po['TP_single'], TP_diluted=po['TP_diluted'],
                    TP_twin=po['TP_on_twin'], dominant=r['dominant_FP'], PTP_mean=r.get('PTP_mean'))
    print(f"{tag:40s} FPP {r['FPP']:.4f} strict {r['FPP_strict']:.4f} PTP {po['PTP_companion']:.2e} BD_M3V {po['BD_on_M3V_single']:.2e} BD_comp {po['BDHEB_on_companion']:.2e} BDtot {bd:.2e} dom {r['dominant_FP']}", flush=True)
for ds in ('ref_beta', 'inv_beta', 'ref_white', 'inv_white'):
    for ecc in ('kipping', 'vaneylen', 'circ'):
        for patched in (False, True):
            run(f'{ds}|{ecc}|{"fixed" if patched else "as_published"}', patched, ds=ds, ecc=ecc)
json.dump(dict(method=__doc__, batman_check=chk, results=res), open(f'{OUT}/b10_r05_eccsign.json', 'w'), indent=1, default=float)
