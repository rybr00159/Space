"""f05: seeing-limited contamination of a circular aperture on the M3V by the Mira's PSF wings, pair A/B, the 3.1"/5.5"
stars and all other Gaia DR3 sources within 25" (2027.5 positions, proper motions applied).
PSF: circular Moffat, beta 2.5/3/4, FWHM 1.0/1.5/2.0/2.5"; aperture radii 1.5/2/3/4". Variant 'halo': Moffat + 5% of the
light in an r^-3 aureole beyond 2 FWHM (pessimistic far wings; not a measured PSF).
Bands g r i z V Ic J and the custom I+z, ASTEP blue/red, ExTrA; the Mira at its bright and faint ends (out/mira.json ranges)
and at its predicted phase for the 2027 transits (out/mira_phase_at_transits.csv; G-band offset applied to all bands).
Outputs: out/contamination_grid.csv, out/contamination_summary.json, out/geometry_2027.csv"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, json
g = pd.read_csv(f'{F}/data/gaia_25arcsec_gspc.csv')
A = g[g.source_id == 5925373272850782208].iloc[0]
def pos2027(r, jy=2027.5):
    pmra, pmdec = r.pmra, r.pmdec
    if not np.isfinite(pmra):
        pmra, pmdec = (A.pmra, A.pmdec) if r.source_id == 5925373272827253632 else (0.0, 0.0)
    dt = jy - 2016.0
    return r.ra + pmra*dt/3.6e6/np.cos(np.radians(r.dec)), r.dec + pmdec*dt/3.6e6
m3 = g[g.source_id == 5925373272850784768].iloc[0]
ra0, de0 = pos2027(m3)
geo = []
for r in g.itertuples():
    ra, de = pos2027(r)
    dE = (ra - ra0)*np.cos(np.radians(de0))*3600; dN = (de - de0)*3600
    geo.append(dict(source_id=r.source_id, dE=dE, dN=dN, sep=np.hypot(dE, dN), pa=np.degrees(np.arctan2(dE, dN)) % 360,
                    G=r.phot_g_mean_mag, bp_rp=r.bp_rp))
GEO = pd.DataFrame(geo).sort_values('sep'); GEO.to_csv(f'{OUT}/geometry_2027.csv', index=False)
names = {5925373268528684288: 'Mira', 5925373272850782208: 'pairA', 5925373272827253632: 'pairB', 5925373272827260544: 'star3.1',
         5925373272827262208: 'star5.5', 5925373272827258752: 'star4.3_G18.6', 5925373272827256576: 'star12_G17.0'}
print(GEO[GEO.source_id.isin(list(names))].assign(name=lambda d: d.source_id.map(names))[['name', 'sep', 'pa', 'G']].to_string())
mira_pa_from_m3v = float(GEO[GEO.source_id == 5925373268528684288].pa.iloc[0])
# ---- magnitudes per band (AB for g r i z and custom bands; Vega for V Ic J): adopted values (f02) or colour estimates
ad = json.load(open(f'{OUT}/mags_adopted.json'))['adopted']
mab = json.load(open(f'{OUT}/m3v_band_mags.json'))['AB']
mira = json.load(open(f'{OUT}/mira.json'))['ranges']
BANDS = ['g', 'r', 'i', 'z', 'V', 'Ic', 'J', 'I+z', 'ASTEP_blue', 'ASTEP_red', 'ExTrA']
# custom bands for the M3V: model AB with the offset of the neighbouring observed band
m3v = {b: ad['M3V'][b]['mag'] for b in ('g', 'r', 'i', 'z', 'V', 'Ic', 'J')}
m3v['I+z'] = mab['I+z'] + (m3v['z'] - mab['z'])
m3v['ASTEP_blue'] = mab['ASTEP_blue'] + 0.5*((m3v['g'] - mab['g']) + (m3v['V'] + 0.02 - mab['V']))
m3v['ASTEP_red'] = mab['ASTEP_red'] + 0.5*((m3v['i'] - mab['i']) + (m3v['z'] - mab['z']))
m3v['ExTrA'] = mab['ExTrA'] + (m3v['J'] + 0.91 - mab['J'])
def proxy(star):
    """band magnitudes of a neighbour; colour offsets relative to the M3V-band system where no measurement exists"""
    d = {b: (ad[star][b]['mag'] if b in ad[star] and ad[star][b]['mag'] is not None else np.nan) for b in ('g', 'r', 'i', 'z', 'V', 'Ic', 'J')}
    if star == 'pairA' and not np.isfinite(d['J']):
        d['J'] = 11.36 + (14.654 - 14.281) - 0.12   # 2MASS blended (Qflg FDC), VHS saturated: M3V J + G offset - small colour term; +-0.2
    if star in ('star3.1', 'star5.5'):
        d['J'] = d['z'] - 0.9 + 0.0         # z_AB - J_Vega ~ 0.9 for a reddened K dwarf (VHS J flagged unreliable); +-0.3
    if star == 'pairB':                     # B: G 15.74; colours of an M4V (= M3V colours + 0.1 redder); measured A+B blends are not used
        dG = 15.74 - 14.281
        d = {b: m3v[b] + dG - {'g': 0.25, 'r': 0.2, 'i': 0.05, 'z': -0.05, 'V': 0.2, 'Ic': 0.0, 'J': -0.15}[b] for b in d}
    d['I+z'] = d['z']; d['ASTEP_blue'] = 0.5*(d['g'] + d['V']); d['ASTEP_red'] = 0.5*(d['i'] + d['z']); d['ExTrA'] = d['J'] + 0.91
    return d
NB = {s: proxy(s) for s in ('pairA', 'pairB', 'star3.1', 'star5.5')}
# Mira: bright and faint ends per band
mr = lambda key: (mira[key]['bright'], mira[key]['faint'])
MIRA = {'g': mr('g (AB)'), 'r': mr('r (AB)'), 'i': mr('i (AB)'), 'z': mr('z (AB)'), 'V': mr('V (Vega)'), 'Ic': mr('Ic (Vega)'), 'J': mr('J (Vega)')}
MIRA['I+z'] = MIRA['z']; MIRA['ASTEP_blue'] = (0.5*(MIRA['g'][0] + MIRA['V'][0]), 0.5*(MIRA['g'][1] + MIRA['V'][1]))
MIRA['ASTEP_red'] = (0.5*(MIRA['i'][0] + MIRA['z'][0]), 0.5*(MIRA['i'][1] + MIRA['z'][1])); MIRA['ExTrA'] = (MIRA['J'][0] + 0.91, MIRA['J'][1] + 0.91)
# other faint Gaia sources: G and a flat colour relative to G (they contribute < 1e-3 in any case)
others = GEO[~GEO.source_id.isin(list(names) + [5925373272850784768]) & (GEO.sep < 25)]
def frac_halo(d, r_ap, fwhm, beta, fh=0.05):
    core = (1 - fh)*moffat_frac_in_circle(d, r_ap, fwhm, beta)
    r0 = 2*fwhm; H = fh*r0/(2*np.pi)
    # halo surface brightness H r^-3 for r >= r0, integrated over the aperture (numerically)
    n = 200; rr = (np.arange(n) + 0.5)/n*r_ap; th = (np.arange(360) + 0.5)/360*2*np.pi
    R, TH = np.meshgrid(rr, th, indexing='ij'); x = R*np.cos(TH) - d; y = R*np.sin(TH); s = np.hypot(x, y)
    S = np.where(s >= r0, H*np.clip(s, 1e-3, None)**-3, 0.0)
    return core + float(np.sum(S*R)*(r_ap/n)*(2*np.pi/360))
sepof = {nm: float(GEO[GEO.source_id == sid].sep.iloc[0]) for sid, nm in names.items()}
rows = []
for fwhm in (1.0, 1.5, 2.0, 2.5):
    for beta in (2.5, 3.0, 4.0):
        for rap in (1.5, 2.0, 3.0, 4.0):
            enc = moffat_enclosed(rap, fwhm, beta)
            fr = {nm: moffat_frac_in_circle(sepof[nm], rap, fwhm, beta) for nm in sepof}
            frh = {nm: frac_halo(sepof[nm], rap, fwhm, beta) for nm in ('Mira', 'pairA', 'pairB')}
            fo = np.array([moffat_frac_in_circle(s, rap, fwhm, beta) for s in others.sep.values])
            for b in BANDS:
                F0 = 10**(-0.4*m3v[b])*enc
                c = {}
                for nm in ('pairA', 'pairB', 'star3.1', 'star5.5'):
                    c[nm] = 10**(-0.4*NB[nm][b])*fr[nm]/F0
                c['Mira_bright'] = 10**(-0.4*MIRA[b][0])*fr['Mira']/F0
                c['Mira_faint'] = 10**(-0.4*MIRA[b][1])*fr['Mira']/F0
                c['Mira_bright_halo'] = 10**(-0.4*MIRA[b][0])*frh['Mira']/F0
                c['pairAB_halo'] = (10**(-0.4*NB['pairA'][b])*frh['pairA'] + 10**(-0.4*NB['pairB'][b])*frh['pairB'])/F0
                # faint field stars: G magnitude with the M3V's G - band offset removed (neutral colour, upper-ish bound)
                off = m3v[b] - 14.281 if b not in ('g', 'V', 'ASTEP_blue') else 0.0
                c['others'] = float(np.sum(10**(-0.4*(others.G.values + off))*fo)/F0)
                tot = c['pairA'] + c['pairB'] + c['star3.1'] + c['star5.5'] + c['Mira_bright'] + c['others']
                rows.append(dict(fwhm=fwhm, beta=beta, r_ap=rap, band=b, enclosed_M3V=enc, **{f'c_{k}': v for k, v in c.items()},
                                 c_total_mira_bright=tot, c_total_mira_bright_halo=tot - c['Mira_bright'] - c['pairA'] - c['pairB'] + c['Mira_bright_halo'] + c['pairAB_halo']))
C = pd.DataFrame(rows); C.to_csv(f'{OUT}/contamination_grid.csv', index=False)
# summary
def pick(fwhm, rap, beta=2.5):
    s = C[(C.fwhm == fwhm) & (C.r_ap == rap) & (C.beta == beta)].set_index('band')
    return {b: dict(total=round(float(s.loc[b, 'c_total_mira_bright']), 5), total_halo=round(float(s.loc[b, 'c_total_mira_bright_halo']), 5),
                    mira_bright=float(s.loc[b, 'c_Mira_bright']), star3p1=round(float(s.loc[b, 'c_star3.1']), 5), star5p5=round(float(s.loc[b, 'c_star5.5']), 5),
                    pair=float(s.loc[b, 'c_pairA'] + s.loc[b, 'c_pairB'])) for b in BANDS}
summ = {f'fwhm{fw}_rap{ra}_beta2.5': pick(fw, ra) for fw, ra in ((1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (2.5, 4.0), (1.5, 4.0), (2.5, 2.0))}
worst = C.sort_values('c_total_mira_bright_halo').iloc[-1].to_dict()
# differential dilution between bands for the depth-ratio test (same aperture): worst case over the grid
dd = []
for (fw, be, ra), s in C.groupby(['fwhm', 'beta', 'r_ap']):
    s = s.set_index('band')
    for b1, b2 in (('g', 'z'), ('V', 'Ic'), ('r', 'i'), ('g', 'i'), ('ASTEP_blue', 'ASTEP_red'), ('i', 'I+z')):
        dd.append(dict(fwhm=fw, beta=be, r_ap=ra, pair=f'{b1}/{b2}', ratio_bias=float((1 + s.loc[b2, 'c_total_mira_bright_halo'])/(1 + s.loc[b1, 'c_total_mira_bright_halo']) - 1)))
DD = pd.DataFrame(dd)
jdump(dict(geometry={nm: dict(sep_2027p5=sepof[nm], pa_deg=float(GEO[GEO.source_id == sid].pa.iloc[0])) for sid, nm in names.items()},
           mira_pa_from_m3v_deg=mira_pa_from_m3v, m3v_mags=m3v, neighbour_mags=NB, mira_mag_ranges=MIRA,
           summary=summ, worst_case=worst,
           diff_dilution_ratio_bias_max_abs={p: float(DD[DD.pair == p].ratio_bias.abs().max()) for p in DD.pair.unique()},
           diff_dilution_ratio_bias_rap_le2_fwhm_le1p5={p: float(DD[(DD.pair == p) & (DD.r_ap <= 2) & (DD.fwhm <= 1.5)].ratio_bias.abs().max()) for p in DD.pair.unique()},
           notes=['c = contaminating flux in the aperture / M3V flux in the aperture; depth is diluted by 1/(1+c)',
                  'pair B magnitudes are estimated from Gaia G (15.74) with M4V-like colours; the 3.1"/5.5" J values from z - 0.9 (VHS values flagged)',
                  'halo variant: 5% of light in an r^-3 aureole beyond 2 FWHM (pessimistic assumption, not a measured PSF)']),
      f'{OUT}/contamination_summary.json')
for k, v in summ.items():
    print(k, {b: (v[b]['total'], v[b]['star3p1']) for b in ('g', 'i', 'z', 'Ic', 'J', 'ExTrA')})
print('worst', {k: worst[k] for k in ('fwhm', 'beta', 'r_ap', 'band', 'c_total_mira_bright_halo')})
print('max |ratio bias|', DD.groupby('pair').ratio_bias.apply(lambda x: x.abs().max()).to_dict())
print('Mira PA from M3V', mira_pa_from_m3v)
for k, v in summ.items():
    print(k, 'Mira bright', {b: f"{v[b]['mira_bright']:.1e}" for b in ('g', 'i', 'z', 'J')}, 'pair', {b: f"{v[b]['pair']:.1e}" for b in ('g', 'i', 'z', 'J')}, 'halo tot', {b: v[b]['total_halo'] for b in ('i', 'z')})
