"""Diagnostics for the WASP fit: folded binned light curve, per-event amplitudes,
null distribution shape, jackknife over events."""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/scripts')
import chromlib as C
import survey_fit as SF
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = SF.OUT
tag = sys.argv[1] if len(sys.argv) > 1 else 'wasp_target'
EXPT = float(sys.argv[2]) / 86400 if len(sys.argv) > 2 else 30 / 86400
z = np.load(OUT + f'/results/{tag}_cleaned.npz')
t, f, w, sg = z['t'], z['f'], z['w'], z['seg']
R = json.load(open(OUT + f'/results/{tag}.json'))['default']
ref = SF.REF
Dd, Ds = ref['combined_deep']['mean'], ref['combined_shallow']['mean']

tem = [SF.tmpl(t, 'd', 0, EXPT), SF.tmpl(t, 's', 0, EXPT)]
a, ae, r, m = C.lsq(t, f, w, sg, tem, 1)
base = m + a[0] * (-tem[0]) * -1 + a[1] * (-tem[1]) * -1  # model = base - a0 s0 - a1 s1 -> base = model + a0 s0 + a1 s1
base = m + a[0] * tem[0] + a[1] * tem[1]
fn = f / base
ph = C.phase(t)

# per-event amplitudes (each event = one segment containing in-eclipse points); fit only that segment
ev = {'d': [], 's': []}
for which, c in (('d', 0.0), ('s', 0.5)):
    s_ = tem[0] if which == 'd' else tem[1]
    for k in np.unique(sg[s_ > 0.05]):
        mm = sg == k
        if (s_[mm] > 0.3).sum() < 3:
            continue
        # need out-of-eclipse points too
        if (s_[mm] == 0).sum() < 10:
            continue
        aa, ee, rr, _ = C.lsq(t[mm], f[mm], w[mm], np.zeros(mm.sum(), int), [s_[mm]], 1)
        ev[which].append((np.median(t[mm]), aa[0], ee[0] * np.sqrt(max(1, np.sum(w[mm] * rr ** 2) / (mm.sum() - 3))), int((s_[mm] > 0.05).sum())))
ev = {k: np.array(v) for k, v in ev.items()}
out = {}
for k in ('d', 's'):
    v = ev[k]
    wm = np.sum(v[:, 1] / v[:, 2] ** 2) / np.sum(1 / v[:, 2] ** 2)
    chi2 = np.sum((v[:, 1] - wm) ** 2 / v[:, 2] ** 2)
    out[k] = dict(n_events=len(v), wmean=wm, chi2=chi2, dof=len(v) - 1, median=float(np.median(v[:, 1])),
                  events=v.tolist())
    print(k, 'events', len(v), 'weighted mean %.0f ppm chi2 %.1f/%d median %.0f' % (wm * 1e6, chi2, len(v) - 1, np.median(v[:, 1]) * 1e6))
    # jackknife: drop one event at a time from global fit
jk = {'d': [], 's': []}
for which in ('d', 's'):
    s_ = tem[0] if which == 'd' else tem[1]
    evs = np.unique(sg[s_ > 0.05])
    for k in evs:
        mm = sg != k
        aa, _, _, _ = C.lsq(t[mm], f[mm], w[mm], sg[mm], [tem[0][mm], tem[1][mm]], 1)
        jk[which].append(aa[0 if which == 'd' else 1])
    v = np.array(jk[which]); n = len(v)
    out[which]['jackknife_sd'] = float(np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2)))
    out[which]['jackknife_range'] = [float(v.min()), float(v.max())]
    print(which, 'jackknife sd %.0f ppm; range %.0f..%.0f' % (out[which]['jackknife_sd'] * 1e6, v.min() * 1e6, v.max() * 1e6))
# null shape
nd = np.array(R['null_samples_deep']); ns = np.array(R['null_samples_shallow'])
for nm, v, obs in (('deep', nd, R['deep']), ('shallow', ns, R['shallow'])):
    zz = (v - v.mean()) / v.std()
    print(nm, 'null: sd %.0f, mad %.0f, frac |z|>2 %.3f (gauss .046), >3 %.4f; obs z %.2f' % (
        v.std() * 1e6, C.mad_std(v) * 1e6, np.mean(np.abs(zz) > 2), np.mean(np.abs(zz) > 3), (obs - v.mean()) / v.std()))
    out[nm + '_null_z_obs'] = float((obs - v.mean()) / v.std())
json.dump(out, open(OUT + f'/results/{tag}_diag.json', 'w'), indent=1, default=float)

# plot: folded binned
fig, ax = plt.subplots(2, 2, figsize=(11, 7))
for k, (c, nm, A, Aref) in enumerate(((0.0, 'deep', a[0], Dd), (0.5, 'shallow', a[1], Ds))):
    dp = (((ph - c) + 0.5) % 1 - 0.5) * C.P * 24
    mm = np.abs(dp) < 4
    bins = np.arange(-4, 4.01, 0.25)
    idx = np.digitize(dp[mm], bins)
    bc, bm, be = [], [], []
    for i in range(1, len(bins)):
        q = idx == i
        if q.sum() < 5:
            continue
        ww = w[mm][q]
        bc.append(0.5 * (bins[i - 1] + bins[i]))
        bm.append(np.sum(ww * (fn[mm][q] - 1)) / np.sum(ww))
        be.append(1 / np.sqrt(np.sum(ww)))
    bc, bm, be = map(np.array, (bc, bm, be))
    axx = ax[0, k]
    axx.errorbar(bc, bm * 1e3, be * 1e3, fmt='o', color='k', label='WASP (15-min bins)')
    x = np.linspace(-4, 4, 800)
    shp = ref['shapes']['d' if k == 0 else 's']
    prof = C.vshape(x / 24 / (shp['T_d'] / 2), shp['p'], shp['q'])
    axx.plot(x, -Aref * prof * 1e3, 'r--', label='TESS depth (ratio 1)')
    axx.plot(x, -A * prof * 1e3, 'b-', label='WASP fit')
    axx.axhline(0, color='gray', lw=0.5)
    axx.set_title(f'{tag}: {nm} eclipse'); axx.set_xlabel('hours from eclipse centre'); axx.set_ylabel('ppt')
    axx.legend(fontsize=8)
    v = ev['d' if k == 0 else 's']
    ax[1, k].errorbar(v[:, 0] - 2450000, v[:, 1] * 1e3, v[:, 2] * 1e3, fmt='o')
    ax[1, k].axhline(Aref * 1e3, color='r', ls='--'); ax[1, k].axhline(0, color='gray')
    ax[1, k].set_xlabel('BJD-2450000'); ax[1, k].set_ylabel('per-event amplitude (ppt)')
plt.tight_layout(); plt.savefig(OUT + f'/plots/{tag}_diag.png', dpi=100)
