"""Pixel-level injection-recovery for the attribution test (one image method per run).

For every sector: real image products of the actual transits, and NOISE draws = the same image products built
from 'fake' transit windows placed at random times away from the real transits (same number of windows as real;
windows within one draw >= 0.7 d apart), i.e. the real pixel noise of that sector incl. the Mira's systematics.
Both image methods are linear in the pixel data, so injecting a transit of total flux loss A on a source with PRF
image p into the cube adds exactly A * p to an image (checked at cube level by r03_cubecheck.py).  An injected
realisation is therefore  noise draw + A_j * PSF_j(position_j + systematic perturbation).
  A_j      = combined real-data amplitude of hypothesis j (weighted over sectors), i.e. 'the observed signal on star j'
  position = predicted apparent position (Gaia + PM + calibrated WCS + lateral-colour curve; Mira: measured), perturbed
             by its systematic uncertainty (WCS 0.02 px/axis; colour curve 0.15/0.3/0.5" radial for C<1.8/<2.3/>=2.3;
             pair B 1.0"; Mira 0.8"/axis)
  PSF      = 'prf' baseline (same model as the fit), or for model-mismatch tests: 'alt' (PRF model of a CCD location
             ~450 px away), 'blur' (PRF smeared by 0.15 px rms), 'emp' (PRF + the fit residual pattern of the red field
             star Gaia DR3 5925373341570270208 (BP-RP 2.41, T 10.5) in that sector, scaled to the injected flux)
For each realisation and sector we store chi2 of the 'a*PSF+b' fit at every hypothesis (93 catalogued stars incl. the
Mira at its measured position), and the free fit on a 1" sky grid (source colour BP-RP 2.76).  Sector sums use the
same realisation index (independent noise in different sectors).
Usage: r01_sim.py METHOD PSFMODE NNOISE NINJ TARGETS(all|main|ring) TAG
Output: out/sim_{TAG}.npz
"""
import sys, time
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
os_threads = 2
import os
os.environ.setdefault('OMP_NUM_THREADS', '2'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '2'); os.environ.setdefault('MKL_NUM_THREADS', '2')
from l2fit import *

meth, psfmode, NN, NI, targets, tag = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
rng = np.random.default_rng(12345)
t0 = time.time()
models = {}
for sec in SECTORS:
    models[sec] = SectorModel(sec, meth)
    print(f'S{sec} model built ({time.time()-t0:.0f}s): n_real {len(models[sec].rv)} n_fake {len(models[sec].fv)}', flush=True)
H0 = models[SECTORS[0]].hyp
for sec in SECTORS:
    assert np.array_equal(models[sec].hyp.source_id.values, H0.source_id.values)
nh = len(H0); names = H0.name.values
# ---------------- real data
real = {}
for sec, M in models.items():
    c2, a, b, sa = M.fit(M.d_real, 'H'); cg, ag, _, sag = M.fit(M.d_real, 'G')
    real[sec] = dict(chi2H=c2, aH=a, saH=sa, chi2G=cg, aG=ag, null=M.null_chi2(M.d_real))
aw = np.array([real[s]['aH'] / real[s]['saH'] ** 2 for s in SECTORS]).sum(0) / np.array([1 / real[s]['saH'] ** 2 for s in SECTORS]).sum(0)
# Injected total flux loss: the observed signal's total (= combined amplitude at the M3V, the best catalogued
# hypothesis), the same for every injected source ('the observed transit signal placed on star j').
# Variant AMPMODE=own: each hypothesis's own best-fit amplitude.
i3_ = int(np.where(H0.name.values == 'M3V')[0][0])
A_inj = np.full(nh, aw[i3_]) if os.environ.get('AMPMODE', 'ref') == 'ref' else aw.copy()
A_own = aw
# ---------------- noise draws (shared by all injections; paired design)
noise = {sec: np.array([M.fake_combo(M.fake_draw(rng)) for _ in range(NN)]) for sec, M in models.items()}
print(f'noise draws done ({time.time()-t0:.0f}s)', flush=True)
# ---------------- injection targets
if targets == 'all':
    tlist = [(i, names[i]) for i in range(nh)]
elif targets == 'mira':
    tlist = [(int(np.where(names == 'Mira')[0][0]), 'Mira')]
elif targets == 'main':
    want = ['M3V', 'pairA', 'pairB', 'Mira', '5925373272827260544', '5925373272827262208', '5925373272827256576']
    tlist = [(int(np.where(names == w)[0][0]), w) for w in want]
else:
    tlist = []
# ring positions around the M3V (uncatalogued-source test): offsets r at 8 PAs, M-dwarf colour
ring = []
if targets == 'ring':
    for r in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0):
        for pa in range(0, 360, 45):
            ring.append((r, pa))
i3 = int(np.where(names == 'M3V')[0][0])


def psf_image(M, x, y, mode, rng):
    p = M.psf(x, y).ravel()
    if mode == 'prf':
        return p
    if mode == 'alt':
        if not hasattr(M, '_altpsf'):
            import PRF
            S = M.S; sub = 's0001' if S.sec < 4 else 's0004'
            dc = -450 if S.col0 + S.xm0 > 1100 else 450; dr = -450 if S.row0 + S.ym0 > 1000 else 450
            pp = PRF.TESS_PRF(S.cam, S.ccd, S.sec, S.col0 + S.xm0 + dc, S.row0 + S.ym0 + dr, localdatadir=f'{PRFDIR}/{sub}/')
            M._altpsf = lambda xx, yy: pp.locate(xx, yy, (S.ny, S.nx))
        return M._altpsf(x, y).ravel()
    if mode == 'blur':
        off = rng.normal(0, 0.15, (16, 2))
        return np.mean([M.psf(x + o[0], y + o[1]).ravel() for o in off], 0)
    if mode == 'emp':
        if not hasattr(M, '_empres'):
            M._empres = empirical_residual(M)
        R, (ix, iy) = M._empres
        img = np.zeros((M.S.ny, M.S.nx))
        sx, sy = int(round(x)) - ix, int(round(y)) - iy
        for (yy, xx), v in np.ndenumerate(R):
            Y, X = yy + sy, xx + sx
            if 0 <= Y < M.S.ny and 0 <= X < M.S.nx: img[Y, X] += v
        return p + img.ravel()
    if mode == 'miraemp':
        # empirical Mira PSF: median image minus all other stars (catalogue fluxes; reference star and background
        # fitted), divided by the Mira's fitted flux.  Used only for injections on the Mira.
        if not hasattr(M, '_miraemp'):
            M._miraemp = empirical_mira(M)
        return M._miraemp.ravel()
    raise ValueError(mode)


def empirical_mira(M):
    from scipy.optimize import least_squares
    S = M.S; g = stars(); xy = S.star_xy(g)
    im = int(np.where(g.source_id == IDS['Mira'])[0][0])
    O = np.median(S.flux[S.good], 0)
    from l1img import background
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    mx, my = M.mira_xy
    mask = np.hypot(xx - mx, yy - my) <= 3.5
    fixed = np.zeros_like(O); ref = None
    for j in np.where(g['T'].values < 18)[0]:
        if j == im: continue
        img = M.psf(*xy[j])
        if g.source_id.values[j] == IDS['ref']: ref = img
        else: fixed += g.F.values[j] * img
    A = np.array([M.psf(mx, my), ref, np.ones_like(O)])
    w = 1 / np.sqrt(np.abs(O[mask]) + (0.02 * O[mask]) ** 2)
    c, *_ = np.linalg.lstsq((A[:, mask].T) * w[:, None], (O - fixed)[mask] * w, rcond=None)
    E = O - fixed - c[1] * ref - c[2]
    E = np.where(np.hypot(xx - mx, yy - my) <= 4.5, E, M.psf(mx, my) * c[0])   # outside 4.5 px: PRF wings
    return E / c[0]


def empirical_residual(M):
    """Residual (data - PRF fit) / fitted flux of the red field star around its peak pixel, 5x5 px."""
    from scipy.optimize import least_squares
    S = M.S; g = stars(); xy = S.star_xy(g)
    j = int(np.where(g.source_id == 5925373341570270208)[0][0]); x0, y0 = xy[j]
    O = np.median(S.flux[S.good], 0); yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    mm = np.hypot(xx - x0, yy - y0) <= 2.6
    nbs = [k for k in np.where(g['T'].values < 16)[0] if k != j and np.hypot(*(xy[k] - xy[j])) < 4.5]
    fixed = sum(g.F.values[k] * M.psf(*xy[k]) for k in nbs)
    def r1(p):
        A = np.column_stack([M.psf(x0 + p[0], y0 + p[1])[mm], np.ones(mm.sum())])
        c, *_ = np.linalg.lstsq(A, (O - fixed)[mm], rcond=None); return A @ c - (O - fixed)[mm]
    rr = least_squares(r1, [0, 0], diff_step=1e-3)
    A = np.column_stack([M.psf(x0 + rr.x[0], y0 + rr.x[1])[mm], np.ones(mm.sum())])
    c, *_ = np.linalg.lstsq(A, (O - fixed)[mm], rcond=None)
    res = np.where(mm, O - fixed - (M.psf(x0 + rr.x[0], y0 + rr.x[1]) * c[0] + c[1]), 0.0) / c[0]
    ix, iy = int(round(x0 + rr.x[0])), int(round(y0 + rr.x[1]))
    R = res[iy - 2:iy + 3, ix - 2:ix + 3]
    return R, (2, 2)


def sys_sigma_arcsec(row):
    if row['name'] == 'Mira': return 0.8, 0.8
    if row['name'] == 'pairB': return 1.0, 0.0
    C = row['C']
    return (0.5 if C >= 2.3 else 0.3 if C >= 1.8 else 0.15), 0.0


def run_target(xy_fun, A, label, ni, batch=200):
    """xy_fun(sec, rng) -> apparent pixel position for this realisation. Returns arrays (batched over realisations)."""
    chiH = np.zeros((ni, len(SECTORS), nh), np.float32)
    aH = np.zeros((ni, len(SECTORS), nh), np.float32)
    free = np.zeros((ni, len(SECTORS), 2), np.float32)
    cmin = np.zeros((ni, len(SECTORS)), np.float32)
    gcomb_min = np.zeros(ni, np.float32); gcomb_pos = np.zeros((ni, 2), np.float32)
    gs_keep = []
    for b0 in range(0, ni, batch):
        b1 = min(ni, b0 + batch); nb = b1 - b0
        gsum = None
        for k, sec in enumerate(SECTORS):
            M = models[sec]
            sig = np.array([psf_image(M, *xy_fun(sec, rng), psfmode, rng)[M.mask] for _ in range(nb)])
            D = noise[sec][b0:b1] + A * sig
            for i in range(nb):
                c2, a, _, _ = M.fit(D[i], 'H'); chiH[b0 + i, k] = c2; aH[b0 + i, k] = a
            p = M.pG; w = M.w
            Sd = D @ w; Spd = D @ p['Mw'].T
            det = p['Spp'] * p['Sw'] - p['Sp'] ** 2
            a = (Spd * p['Sw'] - p['Sp'][None] * Sd[:, None]) / det
            bb = (p['Spp'][None] * Sd[:, None] - p['Sp'][None] * Spd) / det
            Sdd = (D * D * w).sum(1)
            cg = Sdd[:, None] - a * Spd - bb * Sd[:, None]
            del a, bb, Spd
            for i in range(nb):
                pos, cm = M.grid_refine(cg[i]); free[b0 + i, k] = pos; cmin[b0 + i, k] = cm
            gsum = cg if gsum is None else gsum + cg
        for i in range(nb):
            pos, cm = models[SECTORS[0]].grid_refine(gsum[i]); gcomb_min[b0 + i] = cm; gcomb_pos[b0 + i] = pos
        gs_keep.append(GS_REDUCE(gsum))
    return dict(chiH=chiH, aH=aH, free=free, cmin=cmin, gcomb_min=gcomb_min, gcomb_pos=gcomb_pos, gsum=np.concatenate(gs_keep))


KEEP_PTS = [np.zeros(2)] + [r * np.array([np.sin(np.radians(pa)), np.cos(np.radians(pa))]) for r in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0) for pa in range(0, 360, 45)]


def GS_REDUCE(gs):
    return np.column_stack([interp_grid(gs, p) for p in KEEP_PTS])


MAIN = ['M3V', 'pairA', 'pairB', 'Mira', '5925373272827260544', '5925373272827262208', '5925373272827256576']
AX = models[SECTORS[0]].grid_axis


def interp_grid(gs, pt):
    """bilinear interpolation of (ni, ng) grid maps at sky offset pt (relative to the M3V)."""
    n = len(AX); h = AX[1] - AX[0]
    fx = (pt[0] - AX[0]) / h; fy = (pt[1] - AX[0]) / h
    ix, iy = int(np.floor(fx)), int(np.floor(fy)); ux, uy = fx - ix, fy - iy
    G = gs.reshape(len(gs), n, n)
    return ((1 - ux) * (1 - uy) * G[:, iy, ix] + ux * (1 - uy) * G[:, iy, ix + 1] + (1 - ux) * uy * G[:, iy + 1, ix] + ux * uy * G[:, iy + 1, ix + 1])


res = {}
# noise only (A = 0)
r = run_target(lambda sec, rng: models[sec].hyp.loc[i3, ['x', 'y']].values.astype(float), 0.0, 'noise', NN)
r['dchi_m3v'] = r['gsum'][:, 0] - r['gcomb_min']
res['noise'] = {k: v for k, v in r.items() if k != 'gsum'}
print(f'noise-only done ({time.time()-t0:.0f}s)', flush=True)
# catalogued targets
for j, nm in tlist:
    row = H0.loc[j]
    sr, st = sys_sigma_arcsec(row)
    def xyf(sec, rng, j=j, sr=sr, st=st):
        M = models[sec]; S = M.S
        x, y = M.hyp.loc[j, ['x', 'y']].values.astype(float)
        dEN = rng.normal(0, 0.02, 2) @ S.J.T * 0 + rng.normal(0, 1, 2) * 0.0
        pert_px = rng.normal(0, 0.02, 2)                                   # WCS
        rad = rng.normal(0, sr); tan = rng.normal(0, st)
        dEN = rad * S.u_rad + tan * np.array([S.u_rad[1], -S.u_rad[0]])
        return x + pert_px[0] + (S.Ji @ dEN)[0], y + pert_px[1] + (S.Ji @ dEN)[1]
    r = run_target(xyf, A_inj[j], nm, NI)
    if nm not in MAIN: r.pop('aH')
    r['dchi_m3v'] = r['gsum'][:, 0] - r['gcomb_min']
    res[f'inj_{nm}'] = {k: v for k, v in r.items() if k != 'gsum'}
    print(f'injected {nm:22s} A={A_inj[j]:6.1f} ({time.time()-t0:.0f}s)', flush=True)
# ring targets
for rr_, pa in ring:
    pr = np.radians(pa)
    def xyf(sec, rng, rr_=rr_, pr=pr):
        M = models[sec]; S = M.S
        x, y = M.hyp.loc[i3, ['x', 'y']].values.astype(float)
        dEN = rr_ * np.array([np.sin(pr), np.cos(pr)])
        pert_px = rng.normal(0, 0.02, 2)
        return x + (S.Ji @ dEN)[0] + pert_px[0], y + (S.Ji @ dEN)[1] + pert_px[1]
    r = run_target(xyf, A_inj[i3], f'ring_{rr_}_{pa}', NI)
    r.pop('aH')
    # combined grid delta-chi2 at the true offset (grid is co-moving with the M3V, colour 2.76)
    kk = 1 + [(a, b) for a in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0) for b in range(0, 360, 45)].index((rr_, pa))
    r['dchi_true'] = r['gsum'][:, kk] - r['gcomb_min']
    r['dchi_m3v'] = r['gsum'][:, 0] - r['gcomb_min']
    res[f'ring_{rr_}_{pa}'] = {k: v for k, v in r.items() if k != 'gsum'}
    print(f'ring r={rr_} pa={pa} ({time.time()-t0:.0f}s)', flush=True)

# real data summary + grid maps
gs_real = np.array([real[s]['chi2G'] for s in SECTORS]).sum(0)
real_out = dict(gsum=gs_real, chiH=np.array([real[s]['chi2H'] for s in SECTORS]), aH=np.array([real[s]['aH'] for s in SECTORS]),
                saH=np.array([real[s]['saH'] for s in SECTORS]), chiG=np.array([real[s]['chi2G'] for s in SECTORS]),
                aG=np.array([real[s]['aG'] for s in SECTORS]), null=np.array([real[s]['null'] for s in SECTORS]))
flat = {}
for k, v in res.items():
    for kk, vv in v.items():
        flat[f'{k}__{kk}'] = vv
for k, v in real_out.items():
    flat[f'real__{k}'] = v
flat['hyp_source_id'] = H0.source_id.values; flat['hyp_name'] = names.astype(str)
flat['hyp_T'] = H0['T'].values; flat['hyp_C'] = H0.C.values; flat['hyp_dE'] = H0.dE.values; flat['hyp_dN'] = H0.dN.values
flat['hyp_F'] = H0.F.values
flat['A_inj'] = A_inj; flat['A_own'] = A_own; flat['grid'] = models[SECTORS[0]].grid; flat['grid_axis'] = models[SECTORS[0]].grid_axis
flat['sectors'] = np.array(SECTORS)
flat['m3v_off'] = np.array([models[s].m3v_off for s in SECTORS])
flat['hyp_xy'] = np.array([models[s].hyp[['x', 'y']].values for s in SECTORS])
flat['infl'] = np.array([models[s].infl_med for s in SECTORS]); flat['n_real'] = np.array([len(models[s].rv) for s in SECTORS])
flat['n_fake'] = np.array([len(models[s].fv) for s in SECTORS])
np.savez_compressed(f'{BASE}/out/sim_{tag}.npz', **flat)
print(f'saved out/sim_{tag}.npz ({time.time()-t0:.0f}s)')
