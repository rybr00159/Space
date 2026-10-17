"""Item 4: is the depth implied on each candidate star physically possible?

For each catalogued star j within 80" (T < 17.6): the signal's total flux loss A (combined amplitude when star j is
assumed to be the source; DIFF and AMP, and the common 14.3 e/s reference) divided by the star's TESS flux
(Gaia -> T, 15000 e/s at T = 10) gives the required depth delta_j.  Star radius/mass: Mira -> 200-500 Rsun; dwarfs
(significant parallax, M_G consistent with the main sequence) -> Mamajek dwarf table interpolated in M_G (M dwarfs) or
BP-RP; otherwise both a dwarf and a giant branch are considered.  Shape constraints from p02_shape.py (matched-filter
light curve): T14 = 1.48 +- 0.08 h, flat fraction T23/T14 >= 0.40 at 2.6 sigma (>= 0.30 at 3.6 sigma).
Tests: (1) delta < 1; (2) dark occultor: k = sqrt(delta) and a flat bottom with T23/T14 >= f_min for some impact
parameter b <= 1-k; (3) duration: T14(b) within 1.32-1.64 h for such b, circular orbit at P = 5.6915 d.
Output: out/physics.json, out/physics_table.csv"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l0lib import *

# Mamajek dwarf table
rows = []
for ln in open('<local>/tic713/localise/v4/ext/mamajek.txt'):
    if ln.startswith('#') or not ln.strip(): continue
    p = ln.split()
    if len(p) < 30: continue
    if not p[0][0] in 'OBAFGKMLTY': continue
    def f(x):
        try: return float(x)
        except ValueError: return np.nan
    rows.append(dict(spt=p[0], teff=f(p[1]), R=f(p[6]), bprp=f(p[11]), MG=f(p[13]), M=f(p[30])))
    if p[0].startswith("M9.5"): break
mt = pd.DataFrame(rows)
mt = mt[np.isfinite(mt.MG) & np.isfinite(mt.R) & np.isfinite(mt.bprp)]
mtm = mt[np.isfinite(mt.M)]

def dwarf_from_MG(MG):
    o = np.argsort(mtm.MG.values)
    return float(np.interp(MG, mtm.MG.values[o], mtm.R.values[o])), float(np.interp(MG, mtm.MG.values[o], mtm.M.values[o]))

def dwarf_from_col(C):
    o = np.argsort(mtm.bprp.values)
    return float(np.interp(C, mtm.bprp.values[o], mtm.R.values[o])), float(np.interp(C, mtm.bprp.values[o], mtm.M.values[o])), float(np.interp(C, mtm.bprp.values[o], mtm.MG.values[o]))

PD = P
T14_OBS, T14_LO, T14_HI = 1.48, 1.32, 1.64
FMIN = 0.40

def feasible(R, M, delta):
    """dark occultor test: any b in [0, 1-k] with T14 in range and flat fraction >= FMIN"""
    if not np.isfinite(delta) or delta >= 1: return False, np.nan, np.nan
    k = np.sqrt(delta)
    a = (M * (PD / 365.25) ** 2) ** (1 / 3) * 215.032          # Rsun
    if a < 1.5 * R: return False, k, a
    ok = False; best = None
    for b in np.linspace(0, max(0, 1 - k), 200):
        s14 = np.sqrt(max((1 + k) ** 2 - b ** 2, 0)); s23 = np.sqrt(max((1 - k) ** 2 - b ** 2, 0))
        T14 = PD / np.pi * np.arcsin(min(1, R / a * s14)) * 24
        T23 = PD / np.pi * np.arcsin(min(1, R / a * s23)) * 24
        f = T23 / T14 if T14 > 0 else 0
        if T14_LO <= T14 <= T14_HI and f >= FMIN:
            ok = True; best = (b, T14, f); break
    return ok, k, a

Zd = np.load(f'{BASE}/out/sim_diff_prf_all.npz', allow_pickle=True) if os.path.exists(f'{BASE}/out/sim_diff_prf_all.npz') else None
Za = np.load(f'{BASE}/out/sim_amp_prf_all.npz', allow_pickle=True) if os.path.exists(f'{BASE}/out/sim_amp_prf_all.npz') else None
g = stars()
sep = np.hypot(g.dE.values, g.dN.values)
H = g[(sep < 80) & (g['T'].values < 17.6)].reset_index(drop=True)
Aref = float(Zd['A_inj'][0]) if Zd is not None else 14.3
out = []
for i, r in H.iterrows():
    nm = next((k for k, v in IDS.items() if v == r.source_id), str(r.source_id))
    Ad = float(Zd['A_own'][i]) if Zd is not None else np.nan
    Aa = float(Za['A_own'][i]) if Za is not None else np.nan
    F = float(r.F)
    plx, eplx = r.parallax, r.parallax_error
    cls = 'unknown'; R = M = np.nan; R2 = M2 = np.nan
    if nm == 'Mira':
        cls = 'Mira (AGB giant)'; R, M = 300.0, 1.5
    elif nm in ('M3V',):
        cls = 'M3V (TIC)'; R, M = 0.38, 0.37
    elif nm in ('pairA', 'pairB'):
        pl = 16.41
        MG = r.G + 5 * np.log10(pl / 100)
        R, M = dwarf_from_MG(MG); cls = f'M dwarf at 61 pc (M_G {MG:.1f})'
    else:
        if np.isfinite(plx) and plx / eplx > 3:
            MG = r.G + 5 * np.log10(plx / 100)
            Rd, Md, MGd = dwarf_from_col(r.C)
            if MG > MGd - 1.5:
                cls = f'dwarf (M_G {MG:.1f})'; R, M = Rd, Md
            else:
                cls = f'evolved? (M_G {MG:.1f} vs MS {MGd:.1f})'; R, M = Rd * 10 ** (0.2 * (MGd - MG)), 1.2
        else:
            Rd, Md, MGd = dwarf_from_col(r.C)
            cls = 'no parallax: dwarf or giant'; R, M = Rd, Md; R2, M2 = 10.0, 1.3
    rec = dict(name=nm, source_id=int(r.source_id), sep_mira=float(sep[(sep < 80) & (g['T'].values < 17.6)][i]),
               dist_M3V=float(np.hypot(r.dE - g[g.source_id == IDS['M3V']].dE.values[0], r.dN - g[g.source_id == IDS['M3V']].dN.values[0])),
               G=float(r.G), C=float(r.C), T=float(r['T']), F=F, plx=float(plx) if np.isfinite(plx) else None,
               plx_snr=float(plx / eplx) if np.isfinite(plx) and eplx > 0 else None, cls=cls, R=R, M=M,
               A_diff=Ad, A_amp=Aa, A_ref=Aref, depth_ref=Aref / F, depth_diff=Ad / F, depth_amp=Aa / F)
    ok, k, a = feasible(R, M, rec['depth_ref'])
    rec.update(dark_occultor_ok=bool(ok), k=float(k) if np.isfinite(k) else None, a_Rsun=float(a) if np.isfinite(a) else None,
               Rcomp_Rearth=float(k * R * 109.1) if np.isfinite(k) else None)
    if np.isfinite(R2):
        ok2, k2, a2 = feasible(R2, M2, rec['depth_ref']); rec['giant_branch_ok'] = bool(ok2)
    # luminous-total-eclipse alternative: depth = light fraction of a small hot star totally eclipsed
    rec['needs_total_eclipse_of_luminous_companion'] = bool(rec['depth_ref'] < 1 and not ok)
    out.append(rec)
df = pd.DataFrame(out)
df.to_csv(f'{BASE}/out/physics_table.csv', index=False)
main = df[df.name.isin(['M3V', 'pairA', 'pairB', 'Mira']) | (df.dist_M3V < 8) | (df.sep_mira < 9)]
pd.set_option('display.width', 250)
print(main[['name', 'sep_mira', 'dist_M3V', 'T', 'F', 'cls', 'R', 'M', 'A_diff', 'A_amp', 'depth_ref', 'k', 'Rcomp_Rearth', 'dark_occultor_ok']].round(3).to_string())
print(f"\nall {len(df)} stars: depth >= 1 for {(df.depth_ref >= 1).sum()}; dark-occultor feasible for {df.dark_occultor_ok.sum()}: "
      + ', '.join(df[df.dark_occultor_ok].name.astype(str)))
to_json(dict(T14_obs=T14_OBS, T14_range=[T14_LO, T14_HI], flat_min=FMIN, A_ref=Aref, stars=out), f'{BASE}/out/physics.json')
