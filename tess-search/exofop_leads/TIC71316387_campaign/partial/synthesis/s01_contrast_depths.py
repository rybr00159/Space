# Synthesis s01: predicted contrasts (pair vs G star) and multi-band eclipse depths
# from the Mamajek (Pecaut & Mamajek 2013, v2022.04.16) dwarf sequence, for the adopted
# EB components.  Two routes: (a) main-sequence pair at the G star's distance (bound);
# (b) pair anchored to the measured TESS light fraction f_EB, colours from the table.
import numpy as np, json
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
cols = "SpT Teff logT BCv logL Mbol R Mv BV BtVt GV BpRp GRp MG by UB VRc VIc VKs JH HKs MJ MKs KsW1 W1W2 W1W3 W1W4 gr iz zY Msun SpT2".split()
rows = []
for line in open(B + 'astro/data/ext/EEM_dwarf.txt'):
    if line.startswith('#') or not line.strip():
        continue
    p = line.split()
    if len(p) < len(cols) or not (p[0][0] in 'FGKM' and p[0].endswith('V')):
        continue
    def f(x):
        try:
            return float(x)
        except ValueError:
            return np.nan
    rows.append({c: (p[i] if c in ('SpT', 'SpT2') else f(p[i])) for i, c in enumerate(cols)})
T = {c: np.array([r[c] for r in rows]) for c in cols if c not in ('SpT', 'SpT2')}
spt = [r['SpT'] for r in rows]

def tic_T(G, c):  # TIC-8 T(G, BP-RP)
    return G - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473

def absmags_from_mass(M):
    """Interpolate the table in Msun (M dwarfs only, monotonic)."""
    m = np.array([s.startswith('M') or s in ('K7V', 'K8V', 'K9V') for s in spt]) & np.isfinite(T['Msun']) & np.isfinite(T['MKs'])
    o = np.argsort(T['Msun'][m])
    def I(key):
        return float(np.interp(M, T['Msun'][m][o], T[key][m][o]))
    Mv = I('Mv'); MKs = I('MKs'); VIc = I('VIc'); VRc = I('VRc'); JH = I('JH'); HKs = I('HKs')
    MG = I('MG'); BpRp = I('BpRp'); GRp = I('GRp'); KsW1 = I('KsW1'); iz = I('iz'); Teff = I('Teff'); R = I('R')
    out = dict(V=Mv, Rc=Mv - VRc, Ic=Mv - VIc, Ks=MKs, H=MKs + HKs, J=MKs + HKs + JH, G=MG,
               RP=MG - GRp, BP=MG - GRp + BpRp, W1=MKs - KsW1, Teff=Teff, R=R)
    out['T'] = tic_T(MG, BpRp)
    # SDSS via Jordi et al. 2006 (stars): i-Ic = 0.247(R-I)+0.329 ; r-Rc = 0.77(V-R)-0.37 for V-R>=0.93 else 0.267(V-R)+0.088
    RI = VIc - VRc
    out['i'] = out['Ic'] + 0.247 * RI + 0.329
    out['z'] = out['i'] - iz if np.isfinite(iz) else np.nan
    out['r'] = out['Rc'] + (0.77 * VRc - 0.37 if VRc >= 0.93 else 0.267 * VRc + 0.088)
    out['B'] = Mv + I('BV')
    out['g'] = Mv + 0.630 * I('BV') - 0.124
    out['m562'] = out['V'] - 0.10 * (out['V'] - out['Rc'])
    out['m832'] = out['Ic'] - 0.30 * (out['Ic'] - out['z'])
    return out

bands = ['B', 'g', 'V', 'm562', 'r', 'Rc', 'G', 'BP', 'RP', 'i', 'Ic', 'm832', 'T', 'z', 'J', 'H', 'Ks', 'W1']
# observed blend
obs = dict(B=12.493, V=11.857, g=12.141, r=11.745, i=11.627, G=11.775, BP=12.066, RP=11.280, T=11.358,
           J=10.728, H=10.389, Ks=10.341, W1=10.304)
plx = 3.3622; mu = 5 * np.log10(1000 / plx / 10)
AV = 0.031  # E(B-V) ~0.010
Aratio = dict(B=1.32, g=1.20, V=1.0, m562=0.98, r=0.86, Rc=0.81, G=0.84, BP=1.08, RP=0.63, i=0.66, Ic=0.59, m832=0.55,
              T=0.60, z=0.48, J=0.28, H=0.18, Ks=0.12, W1=0.07)

def combine(m1, m2):
    return -2.5 * np.log10(10**(-0.4 * m1) + 10**(-0.4 * m2))

# G-star template colours (Teff 5880 K: interpolate G0V-G1V) relative to Ks
g = np.array([s.startswith('G') or s.startswith('F9') for s in spt])
def gstar_abs(Teff=5880.):
    m = np.array([s in ('F8V', 'F9V', 'F9.5V', 'G0V', 'G1V', 'G2V', 'G3V', 'G4V') for s in spt])
    o = np.argsort(T['Teff'][m])
    I = lambda k: float(np.interp(Teff, T['Teff'][m][o], T[k][m][o]))
    Mv = I('Mv'); MKs = I('MKs'); VIc = I('VIc'); VRc = I('VRc')
    o2 = dict(V=Mv, Rc=Mv - VRc, Ic=Mv - VIc, Ks=MKs, H=MKs + I('HKs'), J=MKs + I('HKs') + I('JH'), G=I('MG'),
              RP=I('MG') - I('GRp'), BP=I('MG') - I('GRp') + I('BpRp'), W1=MKs - I('KsW1'), B=Mv + I('BV'))
    o2['T'] = tic_T(o2['G'], o2['BP'] - o2['RP'])
    RI = VIc - VRc
    o2['i'] = o2['Ic'] + 0.247 * RI + 0.329
    o2['z'] = o2['i'] - 0.20  # solar-type i-z ~0.03-0.2 (table lacks i-z for G); placeholder
    o2['r'] = o2['Rc'] + 0.267 * VRc + 0.088
    o2['g'] = Mv + 0.630 * I('BV') - 0.124
    o2['m562'] = o2['V'] - 0.10 * (o2['V'] - o2['Rc'])
    o2['m832'] = o2['Ic'] - 0.30 * (o2['Ic'] - o2['z'])
    return o2

GS = gstar_abs()
res = {'route_a_bound_MS_at_297pc': {}, 'route_b_fEB_anchored': {}, 'notes': []}
cases = {'adopted (0.46+0.20)': (0.46, 0.20), 'low (0.44+0.19)': (0.44, 0.19), 'high (0.48+0.22)': (0.48, 0.22)}
for lab, (M1, M2) in cases.items():
    a1 = absmags_from_mass(M1); a2 = absmags_from_mass(M2)
    pair_abs = {b: combine(a1[b], a2[b]) for b in bands}
    # route a: pair at 297.4 pc
    app_pair = {b: pair_abs[b] + mu + AV * Aratio[b] for b in bands}
    # G star apparent: anchor to observed Ks minus pair, then table colours
    Ks_G = -2.5 * np.log10(10**(-0.4 * obs['Ks']) - 10**(-0.4 * app_pair['Ks']))
    app_G = {b: Ks_G + (GS[b] - GS['Ks']) + AV * (Aratio[b] - Aratio['Ks']) for b in bands}
    d = {b: app_pair[b] - app_G[b] for b in bands}
    fT = 1 / (1 + 10**(0.4 * d['T']))
    chk = {b: combine(app_pair[b], app_G[b]) - obs[b] for b in obs}
    res['route_a_bound_MS_at_297pc'][lab] = dict(
        M1=M1, M2=M2, SpT1_Teff=a1['Teff'], SpT2_Teff=a2['Teff'], R1_table=a1['R'], R2_table=a2['R'],
        contrast={b: round(d[b], 2) for b in bands}, fEB_T=round(fT, 4),
        pair_frac={b: round(1 / (1 + 10**(0.4 * d[b])), 4) for b in bands},
        secondary_frac_of_pair={b: round(10**(-0.4 * a2[b]) / (10**(-0.4 * a1[b]) + 10**(-0.4 * a2[b])), 4) for b in bands},
        blend_check_model_minus_obs={b: round(chk[b], 3) for b in chk})
    # route b: anchor to f_EB(T)=2.7% (and 2.2, 3.4)
    rb = {}
    for fEB in (0.022, 0.027, 0.034):
        dT = -2.5 * np.log10(fEB / (1 - fEB))
        dd = {b: dT + (pair_abs[b] - pair_abs['T']) - (GS[b] - GS['T']) for b in bands}
        rb[f'f{fEB}'] = {b: round(dd[b], 2) for b in bands}
    res['route_b_fEB_anchored'][lab] = rb

# ---- multi-band eclipse depths (ratios to TESS), adopted case, near-central geometry
# deep eclipse: secondary (k=0.44) transits primary at b~0.14 R1; quadratic LD (approx. Claret-like for 3550 K dwarf)
LD = dict(B=(0.55, 0.28), g=(0.52, 0.30), V=(0.47, 0.32), m562=(0.46, 0.33), r=(0.40, 0.36), Rc=(0.38, 0.37), G=(0.35, 0.37),
          BP=(0.45, 0.33), RP=(0.28, 0.40), i=(0.27, 0.41), Ic=(0.26, 0.41), m832=(0.24, 0.41), T=(0.24, 0.41), z=(0.21, 0.42),
          J=(0.12, 0.40), H=(0.07, 0.38), Ks=(0.06, 0.36), W1=(0.05, 0.30))
def blocked_frac(k, b, u1, u2, n=400):
    # fraction of primary flux blocked by an opaque disc of radius k at distance b (units of R1)
    x = np.linspace(-1, 1, n); X, Y = np.meshgrid(x, x); r2 = X**2 + Y**2
    inside = r2 < 1
    mu_ = np.sqrt(np.clip(1 - r2, 0, 1))
    I = (1 - u1 * (1 - mu_) - u2 * (1 - mu_)**2) * inside
    occ = ((X - b)**2 + Y**2) < k**2
    return float((I * occ).sum() / I.sum())
M1, M2 = cases['adopted (0.46+0.20)']
ra = res['route_a_bound_MS_at_297pc']['adopted (0.46+0.20)']
k, b = 0.44, 0.14
dep = {}
for bnd in bands:
    fp = ra['pair_frac'][bnd]; s2 = ra['secondary_frac_of_pair'][bnd]
    deep = fp * (1 - s2) * blocked_frac(k, b, *LD[bnd])
    shallow = fp * s2  # total occultation of the secondary (b*R1 + R2 < R1)
    dep[bnd] = (deep, shallow)
deepT, shT = dep['T']
obs_deep, obs_sh = 5770., 2600.
res['depths'] = {bnd: dict(deep_ratio_to_TESS=round(dep[bnd][0] / deepT, 3), shallow_ratio_to_TESS=round(dep[bnd][1] / shT, 3),
                           deep_ppm_scaled=round(obs_deep * dep[bnd][0] / deepT), shallow_ppm_scaled=round(obs_sh * dep[bnd][1] / shT),
                           shallow_over_deep=round(dep[bnd][1] / dep[bnd][0], 3)) for bnd in bands}
res['model_TESS_depths_ppm_unscaled'] = dict(deep=round(deepT * 1e6), shallow=round(shT * 1e6), fEB_T=ra['fEB_T'])
res['notes'] = ['m562 ~ V-0.10(V-Rc); m832 ~ Ic-0.30(Ic-z); SDSS g,r,i,z via Jordi+2006 (g uncertain for M dwarfs)',
                'Deep-eclipse band ratios include approximate quadratic LD; +-5-10% systematic',
                'Shallow eclipse treated as total (b R1 + R2 < R1), so its depth = secondary light fraction']
json.dump(res, open(B + 'synthesis/s01_contrast_depths.json', 'w'), indent=1)
for lab in cases:
    r = res['route_a_bound_MS_at_297pc'][lab]
    print(lab, 'Teff', round(r['SpT1_Teff']), round(r['SpT2_Teff']), 'fEB_T', r['fEB_T'])
    print('  contrast', {b: r['contrast'][b] for b in ['V', 'm562', 'G', 'Ic', 'm832', 'T', 'J', 'H', 'Ks']})
    print('  blend check', r['blend_check_model_minus_obs'])
    print('  route b f=2.7%', {b: res['route_b_fEB_anchored'][lab]['f0.027'][b] for b in ['V', 'm562', 'G', 'Ic', 'm832', 'T', 'J', 'Ks']})
print('model TESS', res['model_TESS_depths_ppm_unscaled'])
for bnd in bands:
    print(bnd, res['depths'][bnd])
