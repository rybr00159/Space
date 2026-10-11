"""Box least squares transit search (numba) and per-signal diagnostics.

Light curves are expected detrended and normalised (flux ~ 1). Weights are 1/sigma^2 with sigma the
robust white-noise scatter, so the BLS power at the best box equals the chi^2 improvement and its square
root is the white-noise S/N of the box.
"""
import numpy as np
from numba import njit
from wotan import flatten

DUR_H = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.5])  # trial durations, hours


def freq_grid(baseline, pmin=0.5, pmax=None, os=3.0):
    """Frequency grid with spacing set by the shortest plausible duration at each period."""
    pmax = pmax or baseline / 2.0
    f, fmax, fs = 1.0 / pmax, 1.0 / pmin, []
    while f < fmax:
        fs.append(f)
        q = 0.5 * 0.0756 * f ** (2.0 / 3.0)  # half the duration fraction of a central transit on a Sun-like star
        f += q / (os * baseline)
    return np.array(fs)


@njit(cache=True)
def _bls(t, y, w, freqs, dur_d, binw):
    nf = len(freqs)
    power = np.zeros(nf); depth = np.zeros(nf); t0 = np.zeros(nf); dur = np.zeros(nf)
    W = w.sum()
    for i in range(nf):
        P = 1.0 / freqs[i]
        nb = int(P / binw) + 1
        bw = P / nb
        sw = np.zeros(nb); swy = np.zeros(nb)
        for j in range(len(t)):
            b = int((t[j] % P) / bw)
            if b >= nb: b = nb - 1
            sw[b] += w[j]; swy[b] += w[j] * y[j]
        best = 0.0
        for kd in range(len(dur_d)):
            k = int(round(dur_d[kd] / bw))
            if k < 1: k = 1
            if k > nb // 6: continue  # duty cycle < ~17 %
            s = 0.0; sy = 0.0
            for m in range(k):
                s += sw[m]; sy += swy[m]
            for st in range(nb):
                if sy < 0.0 and s > 0.0 and s < W:
                    p = sy * sy / (s * (1.0 - s / W))
                    if p > best:
                        best = p; depth[i] = -sy / s / (1.0 - s / W); t0[i] = (st + 0.5 * k) * bw; dur[i] = k * bw
                e = (st + k) % nb
                s += sw[e] - sw[st]; sy += swy[e] - swy[st]
        power[i] = best
    return power, depth, t0, dur


def robust_std(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def detrend(t, f, window=0.75):
    trend = flatten(t, f, method="biweight", window_length=window, return_trend=True)[1]
    ok = np.isfinite(trend)
    return t[ok], f[ok] / trend[ok]


def clean(t, f, sig_hi=4.0):
    r = f - np.median(f); s = robust_std(r)
    ok = r < sig_hi * s  # drop upward outliers only
    return t[ok], f[ok], s


def search(t, f, pmin=0.5, pmax=None, os=3.0, binw=10.0 / 1440):
    """Return the strongest box signal and the power spectrum summary."""
    t, f, s = clean(t, f)
    y = f - 1.0
    w = np.full(len(t), 1.0 / s ** 2)
    y = y - np.sum(w * y) / np.sum(w)
    base = t.max() - t.min()
    freqs = freq_grid(base, pmin, pmax, os)
    pw, dp, t0, du = _bls(t - t.min(), y, w, freqs, DUR_H / 24.0, binw)
    i = int(np.argmax(pw))
    # refine around the peak on a 10x finer grid
    df = (freqs[min(i + 1, len(freqs) - 1)] - freqs[max(i - 1, 0)])
    fine = np.linspace(freqs[i] - df, freqs[i] + df, 41)
    fine = fine[fine > 0]
    pw2, dp2, t02, du2 = _bls(t - t.min(), y, w, fine, DUR_H / 24.0, binw / 2)
    j = int(np.argmax(pw2))
    snr_w = np.sqrt(pw2[j])
    # signal detection efficiency on sqrt(power), with a running median removed
    sp = np.sqrt(pw)
    sde = (snr_w - np.median(sp)) / (np.std(sp) + 1e-12)
    return dict(P=1.0 / fine[j], t0=t02[j] + t.min(), dur=du2[j], depth=dp2[j], snr_w=snr_w, sde=sde,
                sigma=s, npts=len(t), baseline=base), (t, f)


def diagnostics(t, f, P, t0, dur, sigma):
    """Per-signal vetting numbers on a detrended light curve."""
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5  # phase in [-0.5, 0.5), transit at 0
    half = 0.5 * dur / P
    intr = np.abs(ph) < half
    oot = np.abs(ph) > 1.5 * half
    base = np.median(f[oot]) if oot.sum() > 10 else 1.0
    out = {}
    epoch = np.round((t - t0) / P).astype(int)
    ev = np.unique(epoch[intr])
    dep, snr, cov = [], [], []
    n_exp = max(1, dur / np.median(np.diff(t)))
    for e in ev:
        m = intr & (epoch == e)
        d = base - np.mean(f[m]); dep.append(d); snr.append(d / (sigma / np.sqrt(m.sum()))); cov.append(m.sum() / n_exp)
    dep, snr, cov = np.array(dep), np.array(snr), np.array(cov)
    full = cov >= 0.5
    out["n_tr"] = int(full.sum())
    out["n_tr_any"] = len(ev)
    tot = np.sum(snr[snr > 0] ** 2) if np.any(snr > 0) else 0.0
    out["max_event_frac"] = float(np.max(snr ** 2) / tot) if tot > 0 else 1.0
    din = base - np.mean(f[intr]); ein = sigma / np.sqrt(max(intr.sum(), 1))
    out["depth"] = din; out["depth_err"] = ein
    # red-noise S/N: scatter of out-of-transit flux binned on the transit duration
    tb = np.floor((t[oot] - t[oot].min()) / dur)
    _, inv, cnt = np.unique(tb, return_inverse=True, return_counts=True)
    mb = np.bincount(inv, f[oot]) / cnt
    good = cnt >= max(2, 0.5 * n_exp)
    sig_d = robust_std(mb[good]) if good.sum() > 5 else sigma / np.sqrt(n_exp)
    out["sigma_dur"] = sig_d
    out["snr_red"] = din / sig_d * np.sqrt(max(out["n_tr"], 1))
    out["red_ratio"] = sig_d / (sigma / np.sqrt(n_exp))
    # odd / even
    odd = intr & (epoch % 2 == 1); even = intr & (epoch % 2 == 0)
    if odd.sum() > 2 and even.sum() > 2:
        do, de = base - f[odd].mean(), base - f[even].mean()
        eo, ee = sigma / np.sqrt(odd.sum()), sigma / np.sqrt(even.sum())
        out["oddeven_sig"] = abs(do - de) / np.hypot(eo, ee)
    else:
        out["oddeven_sig"] = np.nan
    # secondary eclipse at phase 0.5 and the strongest dip anywhere outside the primary
    sec = np.abs(np.abs(ph) - 0.5) < half
    out["sec_depth"] = base - f[sec].mean() if sec.sum() > 2 else np.nan
    out["sec_sig"] = out["sec_depth"] / (sigma / np.sqrt(max(sec.sum(), 1))) if sec.sum() > 2 else np.nan
    # phase-curve uniqueness: other duration-wide phase bins should be flat (rejects pulsators, EBs, spots)
    nbin = max(6, int(round(P / dur)))
    bi = (((t - t0) / P + 0.5 / nbin) % 1.0 * nbin).astype(int) % nbin
    cnt = np.bincount(bi, minlength=nbin); sm = np.bincount(bi, f, minlength=nbin)
    ok_b = cnt[2:nbin - 1] > 2
    if ok_b.sum() >= 3:
        mb = (sm[2:nbin - 1] / np.maximum(cnt[2:nbin - 1], 1))[ok_b]
        eb_ = sigma / np.sqrt(cnt[2:nbin - 1][ok_b])
        dmin = base - mb.min()
        out["oot_dip_ratio"] = dmin / din if din > 0 else np.nan
        out["oot_dip_sig"] = dmin / eb_[np.argmin(mb)]
        out["oot_rms_ratio"] = np.std(mb - base) / np.median(eb_)
    else:
        out["oot_dip_ratio"] = out["oot_dip_sig"] = out["oot_rms_ratio"] = np.nan
    # transit shape: depth in the central half vs the outer half of the box
    cen = np.abs(ph) < 0.5 * half; edge = intr & ~cen
    if cen.sum() > 2 and edge.sum() > 2:
        out["shape_ratio"] = (base - f[edge].mean()) / (base - f[cen].mean() + 1e-12)
    else:
        out["shape_ratio"] = np.nan
    return out


def _local_peaks(pw, freqs, k=10, sep=0.01):
    order = np.argsort(pw)[::-1]; chosen = []
    for i in order:
        if pw[i] <= 0: break
        P = 1.0 / freqs[i]
        if all(abs(P - 1.0 / freqs[j]) / P > sep for j in chosen):
            chosen.append(i)
        if len(chosen) >= k: break
    return chosen


def search_robust(t, f, pmin=0.5, pmax=None, os=3.0, binw=10.0 / 1440, min_tr=2, npeaks=10, max_mask=3):
    """BLS with single-event masking and multi-peak selection.

    Isolated single events (one transit carrying > 60 % of the S/N) are masked and the search rerun, up
    to `max_mask` times. Among the top `npeaks` periods, the best by red-noise S/N with at least
    `min_tr` full transits and no dominant single event is returned.
    """
    t, f, s = clean(t, f)
    base = t.max() - t.min()
    freqs = freq_grid(base, pmin, pmax, os)
    keep = np.ones(len(t), bool); masked = 0
    while True:
        tt, ff = t[keep], f[keep]
        y = ff - 1.0; w = np.full(len(tt), 1.0 / s ** 2); y = y - np.sum(w * y) / np.sum(w)
        pw, dp, t0, du = _bls(tt - tt.min(), y, w, freqs, DUR_H / 24.0, binw)
        sp = np.sqrt(pw); med, sd = np.median(sp), np.std(sp) + 1e-12
        cands = []
        for i in _local_peaks(pw, freqs, npeaks):
            df = freqs[min(i + 1, len(freqs) - 1)] - freqs[max(i - 1, 0)]
            fine = np.linspace(freqs[i] - df, freqs[i] + df, 15); fine = fine[fine > 0]
            durs = np.array([du[i] / 1.4, du[i], du[i] * 1.4])  # around the coarse best duration
            p2, d2, t2, u2 = _bls(tt - tt.min(), y, w, fine, durs, binw)
            j = int(np.argmax(p2))
            c = dict(P=1.0 / fine[j], t0=t2[j] + tt.min(), dur=u2[j], snr_w=np.sqrt(p2[j]), sde=(np.sqrt(p2[j]) - med) / sd)
            c.update(diagnostics(tt, ff, c["P"], c["t0"], c["dur"], s))
            cands.append(c)
        top = cands[0]
        if top["max_event_frac"] > 0.6 and masked < max_mask:
            # mask the dominant event of the top peak and search again
            ph = ((t - top["t0"]) / top["P"] + 0.5) % 1.0 - 0.5
            ep = np.round((t - top["t0"]) / top["P"])
            intr = np.abs(ph) < top["dur"] / top["P"]
            evs, cnt = np.unique(ep[intr & keep], return_counts=True)
            dmax, emax = -1, None
            for e in evs:
                m = intr & keep & (ep == e); d = 1 - np.mean(f[m])
                if d * np.sqrt(m.sum()) > dmax: dmax, emax = d * np.sqrt(m.sum()), e
            keep &= ~(intr & (ep == emax)); masked += 1
            continue
        break
    ok = [c for c in cands if c["n_tr"] >= min_tr and c["depth"] > 0
          and c["max_event_frac"] <= (0.7 if c["n_tr"] == 2 else 0.5)]
    best = max(ok, key=lambda c: c["snr_red"]) if ok else max(cands, key=lambda c: c["snr_red"])
    rank = cands.index(best)
    best = check_aliases(t[keep], f[keep], best, s)
    best.update(sigma=s, npts=int(keep.sum()), baseline=base, n_masked=masked, passed=bool(ok), rank=rank)
    return best, (t[keep], f[keep])


def check_aliases(t, f, best, sigma, min_tr=2):
    """Try 2P and 3P (every epoch offset); keep a multiple if its red-noise S/N is higher."""
    out = best
    for k in (2, 3):
        for j in range(k):
            P = best["P"] * k; t0 = best["t0"] + j * best["P"]
            d = diagnostics(t, f, P, t0, best["dur"], sigma)
            if d["n_tr"] >= min_tr and d["max_event_frac"] <= 0.6 and d["snr_red"] > out["snr_red"] * 1.05:
                out = dict(best); out.update(d); out.update(P=P, t0=t0, alias_of=best["P"])
    return out
