"""Independent cross-check of the SuperWASP depths with the 'detrend with eclipses masked' recipe:
per night x camera, fit a line to points > 1.5 h from both eclipse centres, divide it out, clip at 4 sigma,
then the depth = template-weighted mean of the detrended flux deficit (matched filter, no simultaneous fit).
Null: same recipe at 1000 fake phase offsets (each re-masked at its own fake phases and at the real ones)."""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import wasp3 as W
T = json.load(open(L.V3 + '/results/tess3.json')); E3 = json.load(open(L.V3 + '/results/ephem3.json'))
shp = L.survey_shapes(T['shapes'], E3)
t, f, seg, cam, hjd = W.load(W.TARGET)
EX = 30 / 86400


def depth_masked(dphi):
    cen = []
    for which in 'ds':
        for d in ((dphi, 0.0) if dphi != 0 else (0.0,)):
            cen.append(np.abs(L.dphase(t, shp[which]['phc'] + d, shp[which]['per'], shp[which]['t0'])))
    far = np.min(np.array(cen), axis=0) > 1.5 / 24
    fd = np.full(len(t), np.nan)
    wt = np.zeros(len(t))
    for k in np.unique(seg):
        m = seg == k
        mm = m & far
        if mm.sum() < 10:
            continue
        c = np.polyfit(t[mm] - t[m].mean(), f[mm], 1)
        fd[m] = f[m] / np.polyval(c, t[m] - t[m].mean())
        r = fd[mm] - 1
        s = L.mad_std(r)
        bad = m & (np.abs(fd - 1) > 4 * s)
        fd[bad] = np.nan
        wt[m] = 1 / max(s, 2e-4) ** 2
    ok = np.isfinite(fd)
    out = []
    for which in 'ds':
        s = L.tmpl(t[ok], shp, which, dphi, EX)
        y = 1 - fd[ok]
        # per-point weights from the night's robust scatter
        ww = wt[ok]
        out.append(np.sum(ww * s * y) / np.sum(ww * s * s))
    return np.array(out)


obs = depth_masked(0.0)
rng = np.random.default_rng(4)
nl = np.array([depth_masked(d) for d in L.fake_offsets(400, rng)])
print('masked-detrending matched filter: deep %.0f +- %.0f ppm (null mean %.0f)  shallow %.0f +- %.0f ppm (null mean %.0f); null corr %.2f' % (
    obs[0] * 1e6, nl[:, 0].std() * 1e6, nl[:, 0].mean() * 1e6, obs[1] * 1e6, nl[:, 1].std() * 1e6, nl[:, 1].mean() * 1e6, np.corrcoef(nl.T)[0, 1]))
# injection: template-shaped 3000 ppm at fake phases
inj = []
for d in L.fake_offsets(60, rng):
    f0 = f.copy()
    f[:] = f0 * (1 - 0.003 * L.tmpl(t, shp, 'd', d, EX) - 0.0014 * L.tmpl(t, shp, 's', d, EX))
    a1 = depth_masked(d)
    f[:] = f0
    a0 = depth_masked(d)
    inj.append(((a1 - a0) / np.array([0.003, 0.0014])))
inj = np.array(inj)
print('injection gain (masked recipe): deep %.3f +- %.3f, shallow %.3f +- %.3f' % (inj[:, 0].mean(), inj[:, 0].std() / np.sqrt(len(inj)), inj[:, 1].mean(), inj[:, 1].std() / np.sqrt(len(inj))))
json.dump(dict(deep=float(obs[0]), shallow=float(obs[1]), deep_err=float(nl[:, 0].std()), shallow_err=float(nl[:, 1].std()),
               gain=inj.mean(0).tolist()), open(L.V3 + '/results/masked3.json', 'w'), indent=1)
