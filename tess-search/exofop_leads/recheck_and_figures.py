"""Independent re-check and figures for the CTOI eclipsing-binary leads.

For each star: fetch TESS light curves fresh from MAST (SPOC > TESS-SPOC > QLP per sector),
detrend with eclipses masked, find the epoch at the proposed period, and measure the depths of
the eclipses at phase 0 and 0.5 event by event. Writes one PNG per star and a JSON summary.
"""
import json, sys, warnings, os
import numpy as np
import lightkurve as lk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
OUT = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(OUT, 'cache')

# P = proposed true (eclipsing-binary) period; dur in hours (approximate, refined below)
LEADS = [
    dict(ctoi='327712163.01', tic=327712163, P=0.7807704, dur=1.0, kind='alt', ctoiP=0.39038),
    dict(ctoi='206599822.01', tic=206599822, P=2.151597, dur=2.0, kind='alt', ctoiP=1.076),
    dict(ctoi='71316387.01', tic=71316387, P=1.303671, dur=1.5, kind='alt', ctoiP=0.651837),
    dict(ctoi='147285336.01', tic=147285336, P=3.36999, dur=2.4, kind='alt', ctoiP=1.684942),
    dict(ctoi='234503282.01', tic=234503282, P=2.876461, dur=2.8, kind='alt', ctoiP=1.4376),
    dict(ctoi='65450588.01', tic=65450588, P=1.2973383, dur=1.5, kind='alt', ctoiP=0.64872),
    dict(ctoi='382360740.01', tic=382360740, P=5.131133, dur=4.0, kind='sec', ctoiP=5.15),
    dict(ctoi='102064145.01', tic=102064145, P=4.392547, dur=4.3, kind='neb', ctoiP=4.39187),
    dict(ctoi='218856662.01', tic=218856662, P=1.1295000, dur=3.0, kind='neb', ctoiP=1.129424),
    dict(ctoi='100918213.01', tic=100918213, P=0.5128982, dur=1.3, kind='neb', ctoiP=0.25643),
]


def fetch(tic):
    sr = lk.search_lightcurve(f'TIC {tic}', mission='TESS')
    best = {}
    rank = {'SPOC': 0, 'TESS-SPOC': 1, 'QLP': 2}
    for i, row in enumerate(sr.table):
        a = str(row['author'])
        if a not in rank:
            continue
        sec = int(row['sequence_number'])
        exptime = float(row['exptime'])
        if exptime < 100:  # skip 20-s products; 120-s and FFI cadences are enough
            continue
        key = (rank[a], exptime)
        if sec not in best or key < best[sec][0]:
            best[sec] = (key, i, a)
    lcs = []
    for sec, (key, i, a) in sorted(best.items()):
        try:
            lc = sr[i].download(download_dir=CACHE, quality_bitmask='default')
        except Exception as e:  # noqa
            print('  download failed', tic, sec, a, e, flush=True)
            continue
        col = 'pdcsap_flux' if a in ('SPOC', 'TESS-SPOC') else 'sap_flux'
        if col not in lc.colnames:
            col = 'flux'
        f = np.asarray(lc[col].value, float)
        t = np.asarray(lc.time.value, float)
        ok = np.isfinite(f) & np.isfinite(t) & (f > 0)
        if ok.sum() < 200:
            continue
        lcs.append(dict(sector=sec, author=a, t=t[ok], f=f[ok] / np.nanmedian(f[ok])))
    return lcs


def running_median(t, f, w, mask):
    # median on a coarse grid (every 0.02 d) of unmasked points, then interpolate
    tt, ff = t[~mask], f[~mask]
    grid = np.arange(t.min(), t.max() + 0.02, 0.02)
    gv = np.full(grid.size, np.nan)
    lo = np.searchsorted(tt, grid - w / 2); hi = np.searchsorted(tt, grid + w / 2)
    for i in range(grid.size):
        if hi[i] - lo[i] > 5:
            gv[i] = np.median(ff[lo[i]:hi[i]])
    good = np.isfinite(gv)
    out = np.interp(t, grid[good], gv[good])
    # do not extrapolate across large gaps
    near = np.abs(t[:, None] - grid[good][None, :]).min(axis=1) if t.size < 30000 else np.zeros(t.size)
    out[near > w] = np.nan
    return out


def phase(t, P, T0):
    return ((t - T0) / P + 0.5) % 1 - 0.5


def detrend(lc, P, T0, dur_d):
    t, f = lc['t'], lc['f']
    w = max(0.75, 4 * dur_d)
    if P < 1.0:
        w = max(0.5, 4 * dur_d)
    mask = np.zeros(t.size, bool)
    if T0 is not None:
        ph = phase(t, P, T0)
        ph2 = phase(t, P, T0 + P / 2)
        mask = (np.abs(ph * P) < 0.75 * dur_d) | (np.abs(ph2 * P) < 0.75 * dur_d)
    trend = running_median(t, f, w, mask)
    return f / trend


def find_T0(t, f, P, dur_d):
    # deepest dip in a box search over phase
    grid = np.linspace(0, P, 2000, endpoint=False)
    best, bt = 0, 0
    for g in grid:
        ph = phase(t, P, t.min() + g)
        inn = np.abs(ph * P) < dur_d / 4
        if inn.sum() < 5:
            continue
        d = 1 - np.nanmedian(f[inn])
        if d > best:
            best, bt = d, t.min() + g
    return bt


def event_depths(t, f, P, T0, dur_d):
    """Depth of each individual eclipse at T0 + n*P: in-eclipse mean (|dt|<dur/4) vs local
    out-of-eclipse ring (0.75..1.5 dur)."""
    n0 = np.floor((t.min() - T0) / P) - 1
    n1 = np.ceil((t.max() - T0) / P) + 1
    deps, ns = [], []
    for n in np.arange(n0, n1 + 1):
        tc = T0 + n * P
        dt = t - tc
        inn = np.abs(dt) < dur_d / 4
        ring = (np.abs(dt) > 0.75 * dur_d) & (np.abs(dt) < 1.5 * dur_d)
        if inn.sum() < 3 or ring.sum() < 6:
            continue
        deps.append(np.nanmedian(f[ring]) - np.nanmean(f[inn]))
        ns.append(n)
    deps, ns = np.array(deps), np.array(ns)
    if deps.size > 4:  # clip single corrupted events (5 robust sigma)
        med = np.median(deps); mad = 1.4826 * np.median(np.abs(deps - med))
        keep = np.abs(deps - med) < 5 * max(mad, 1e-6)
        deps, ns = deps[keep], ns[keep]
    return deps, ns


def analyse(L):
    tic, P = L['tic'], L['P']
    dur_d = L['dur'] / 24
    lcs = fetch(tic)
    if not lcs:
        return dict(L, error='no data')
    # reject broken sectors (e.g. corrupted QLP background) by point-to-point scatter
    for lc in lcs:
        lc['p2p'] = np.nanmedian(np.abs(np.diff(lc['f']))) * 1.48 / np.sqrt(2)
    med = np.median([lc['p2p'] for lc in lcs])
    lcs = [lc for lc in lcs if lc['p2p'] < 5 * med]
    # pass 1 detrend, find epoch; pass 2 masked detrend
    for lc in lcs:
        lc['d'] = detrend(lc, P, None, dur_d)
    t = np.concatenate([lc['t'] for lc in lcs]); f = np.concatenate([lc['d'] for lc in lcs])
    ok = np.isfinite(f)
    T0 = find_T0(t[ok], f[ok], P, dur_d)
    for lc in lcs:
        lc['d'] = detrend(lc, P, T0, dur_d)
    t = np.concatenate([lc['t'] for lc in lcs]); f = np.concatenate([lc['d'] for lc in lcs])
    sec = np.concatenate([np.full(lc['t'].size, lc['sector']) for lc in lcs])
    ok = np.isfinite(f); t, f, sec = t[ok], f[ok], sec[ok]
    # refine epoch on cleaned data
    T0 = find_T0(t, f, P, dur_d)
    dA, nA = event_depths(t, f, P, T0, dur_d)
    dB, nB = event_depths(t, f, P, T0 + P / 2, dur_d)
    res = dict(ctoi=L['ctoi'], tic=tic, P=P, T0_BTJD=round(float(T0), 5), kind=L['kind'],
               sectors=[f"{lc['sector']}:{lc['author']}" for lc in lcs],
               nA=int(dA.size), nB=int(dB.size))
    for lab, d in (('A', dA), ('B', dB)):
        res[f'depth{lab}_ppm'] = float(np.mean(d) * 1e6) if d.size else None
        res[f'err{lab}_ppm'] = float(np.std(d, ddof=1) / np.sqrt(d.size) * 1e6) if d.size > 2 else None
    if dA.size > 2 and dB.size > 2:
        diff = np.mean(dA) - np.mean(dB)
        err = np.hypot(np.std(dA, ddof=1) / np.sqrt(dA.size), np.std(dB, ddof=1) / np.sqrt(dB.size))
        res['diff_sigma'] = float(diff / err)
        # permutation null: shuffle A/B labels
        allD = np.concatenate([dA, dB]); rng = np.random.default_rng(1); null = []
        for _ in range(5000):
            p = rng.permutation(allD.size)
            a, b = allD[p[:dA.size]], allD[p[dA.size:]]
            null.append(np.mean(a) - np.mean(b))
        res['perm_p'] = float((np.sum(np.abs(null) >= abs(diff)) + 1) / 5001)
        res['perm_sigma_equiv'] = float(abs(diff) / np.std(null))
        # per-sector sign consistency
        signs = []
        for s in np.unique(sec):
            ta = set(np.round(T0 + nA * P, 3))
            selA = [d for d, n in zip(dA, nA) if np.any((sec == s) & (np.abs(t - (T0 + n * P)) < dur_d))]
            selB = [d for d, n in zip(dB, nB) if np.any((sec == s) & (np.abs(t - (T0 + P / 2 + n * P)) < dur_d))]
            if len(selA) >= 2 and len(selB) >= 2:
                signs.append((int(s), round(float((np.median(selA) - np.median(selB)) * 1e6))))
        res['per_sector_A_minus_B_ppm'] = signs
    # figure: same y-scale on both panels, clipped points, binned medians
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    win = 2.5 * dur_d
    dmax = max(res['depthA_ppm'] or 0, res['depthB_ppm'] or 0) / 1e6
    oot = np.abs(phase(t, P, T0) * P) > 1.5 * dur_d
    oot &= np.abs(phase(t, P, T0 + P / 2) * P) > 1.5 * dur_d
    sig = 1.4826 * np.median(np.abs(f[oot] - np.median(f[oot])))
    clip = (f > 1 - 2.0 * dmax - 6 * sig) & (f < 1 + 6 * sig + 0.3 * dmax)
    nb = 36
    for k, (off, lab) in enumerate(((0.0, 'phase 0 (deeper eclipse)'), (0.5, 'phase 0.5'))):
        ph = phase(t, P, T0 + off * P) * P
        sel = (np.abs(ph) < win) & clip
        ax[k].plot(ph[sel] * 24, (f[sel] - 1) * 1e3, '.', ms=1.0, color='0.75', rasterized=True)
        bins = np.linspace(-win, win, nb + 1)
        idx = np.digitize(ph[sel], bins)
        bc, bm, be = [], [], []
        for i in range(1, nb + 1):
            v = f[sel][idx == i]
            if v.size > 4:
                bc.append((bins[i - 1] + bins[i]) / 2); bm.append(np.median(v))
                be.append(1.2533 * 1.4826 * np.median(np.abs(v - np.median(v))) / np.sqrt(v.size))
        ax[k].errorbar(np.array(bc) * 24, (np.array(bm) - 1) * 1e3, yerr=np.array(be) * 1e3, fmt='o-', ms=3,
                       lw=1, color='C0' if k == 0 else 'C3')
        dd = res['depthA_ppm'] if k == 0 else res['depthB_ppm']
        ee = res['errA_ppm'] if k == 0 else res['errB_ppm']
        nn = res['nA'] if k == 0 else res['nB']
        ax[k].set_title(f"{lab}: {dd:.0f} ± {ee:.0f} ppm ({nn} events)" if dd is not None and ee is not None else lab,
                        fontsize=9)
        ax[k].set_xlabel('hours from eclipse centre')
        ax[k].axhline(0, color='k', lw=0.5)
        ax[k].axhline(-dmax * 1e3, color='0.5', lw=0.5, ls=':')
    ax[0].set_ylabel('relative flux - 1 (ppt)')
    ax[0].set_ylim(-(1.6 * dmax + 4 * sig) * 1e3, (0.35 * dmax + 4 * sig) * 1e3)
    extra = {'alt': ' (= 2 x CTOI period)', 'sec': '', 'neb': ''}[L['kind']]
    if L['kind'] == 'neb' and abs(P / L['ctoiP'] - 2) < 0.01:
        extra = ' (= 2 x CTOI period)'
    fig.suptitle(f"CTOI {L['ctoi']} (TIC {tic}) folded at P = {P} d{extra}; T0 = BTJD {T0:.4f}\n"
                 f"TESS sectors {', '.join(str(lc['sector']) for lc in lcs)}; per-event depths, single bad events clipped",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"CTOI_{L['ctoi']}.png"), dpi=130)
    plt.close(fig)
    return res


if __name__ == '__main__':
    only = sys.argv[1:]
    results = []
    for L in LEADS:
        if only and L['ctoi'] not in only:
            continue
        print('==', L['ctoi'], flush=True)
        try:
            r = analyse(L)
        except Exception as e:  # noqa
            r = dict(L, error=repr(e))
        print(json.dumps(r), flush=True)
        results.append(r)
    json.dump(results, open(os.path.join(OUT, 'recheck_' + ('_'.join(only) if only else 'all') + '.json'), 'w'), indent=1)
