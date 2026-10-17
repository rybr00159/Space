"""Referee r19: cross-talk matrix of the referee local-window photometry (r05) from pixel-level injections, and an
independent forward-model host test.
Injections: in each sector, NF fake windows (r05 fake centres) receive a 15 e-/s trapezoid (r06 global shape, at the
window centre) using the image of one source: M3V PRF, pair PRF (A + 0.367 B), Mira empirical b-image, Mira SPOC PRF
at its Gaia position (a template the fit does NOT use), bg3 PRF (Gaia DR3 5925373272827260544, 3" from the M3V), and
three deliberately mismatched M3V images (PRF shifted 0.1 px towards bg3; shifted 0.1 px away from bg3; blurred by a
0.15 px Gaussian) to see whether a PRF/WCS error of that size mimics the bg3 response.
Each injected window is re-processed exactly like the real data; the response of the four observables
  aM (M3V-only fit), aT_M3V, aT_pair, aT_Mira (three-star fit)
is (depth_injected - depth_not_injected)/15 per window (noise cancels).
Host test: per sector, the observed four-vector (mean single-window depth over the real strict windows) is compared
with X * response_h, X free and common to all sectors, using the 4x4 covariance of single-window depths at 120+
fake windows divided by the number of real windows. chi2 over 6 sectors x 4 observables (23 dof). For the mismatched
M3V images the expected (noise-free) chi2 of the M3V and bg3 hypotheses is also given (signal 15 e-/s), to compare
with the observed bg3-M3V chi2 difference. Windows with < 2 in-transit cadences are skipped.
Output: out/crosstalk.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import Sector, trap_shape, real_centres, fake_centres
from scipy.stats import chi2 as C2
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]; T0, P, T14, FL = G["T0"], G["P"], G["T14_h"] / 24, G["flat"]
NF = int(sys.argv[1]) if len(sys.argv) > 1 else 40
OBS = ("aM", "T_M3V", "T_pair", "T_Mira")


def depths(r, tc, ex):
    y = trap_shape(r["t"], tc, T14, FL, ex)
    A = np.column_stack([np.ones(len(y)), -y]); out = []
    for v in (r["aM"], r["aT"][:, 0], r["aT"][:, 1], r["aT"][:, 2]):
        c, *_ = np.linalg.lstsq(A, v, rcond=None); out.append(c[1])
    return np.array(out)


resp = {}; obs = {}; cov = {}; nreal = {}
rng = np.random.default_rng(404)
for s in SECTORS:
    S = Sector(s); C = S.C
    m3 = np.array(S.pos["M3V"]); u = np.array(S.pos["bg3"]) - m3; u /= np.hypot(*u)
    n_ = np.array([-2, -1, 0, 1, 2]) * 0.12; wg = np.exp(-0.5 * (n_ / 0.15) ** 2); wg /= wg.sum()
    blur = sum(wg[i] * wg[j] * C.psf(m3[0] + n_[i], m3[1] + n_[j]) for i in range(5) for j in range(5))
    imgs = dict(M3V=S.P_M3V, pair=S.P_pair, Mira_b=S.b_Mira.reshape(C.ny, C.nx), Mira_prf=S.P_Mira, bg3=C.psf(*S.pos["bg3"]),
                M3V_to_bg3=C.psf(*(m3 + 0.1 * u)), M3V_away_bg3=C.psf(*(m3 - 0.1 * u)), M3V_blur=blur)
    # observed
    cen, ks = real_centres(S, P, T0); dd = []
    for tc, k in zip(cen, ks):
        if s == 12 and k == -388: continue
        r = S.window(tc)
        if r is None: continue
        y = trap_shape(r["t"], tc, T14, FL, S.exptime)
        if np.sum(y > 0.5) < 2: continue
        dd.append(depths(r, tc, S.exptime))
    dd = np.array(dd); obs[s] = dd.mean(0); nreal[s] = len(dd)
    # null covariance and responses
    fk = fake_centres(S, 150, rng, P, T0)
    nul = []; R_ = {h: [] for h in imgs}
    for i, tc in enumerate(fk):
        r0 = S.window(tc)
        if r0 is None or np.sum(trap_shape(r0["t"], tc, T14, FL, S.exptime) > 0.5) < 2: continue
        d0 = depths(r0, tc, S.exptime); nul.append(d0)
        if i < NF:
            for h, im in imgs.items():
                r = S.window(tc, inject=[(im, 15.0, tc, T14, FL)])
                if r is None or len(r["t"]) != len(r0["t"]): continue
                R_[h].append((depths(r, tc, S.exptime) - d0) / 15.0)
    nul = np.array(nul); cov[s] = np.cov(nul.T) / nreal[s]
    resp[s] = {h: (np.mean(v, 0), np.std(v, 0)) for h, v in R_.items()}
    print(f"S{s}: n_real {nreal[s]}, obs " + ", ".join(f"{o} {x:.1f}+-{np.sqrt(cov[s][j,j]):.1f}" for j, (o, x) in enumerate(zip(OBS, obs[s]))), flush=True)
    for h, (m, sd) in resp[s].items():
        print(f"   inject {h:8s}: response " + ", ".join(f"{o} {a:+.3f}({b:.3f})" for o, a, b in zip(OBS, m, sd)), flush=True)
# host test
host = {}
for h in ("M3V", "pair", "Mira_b", "Mira_prf", "bg3", "M3V_to_bg3", "M3V_away_bg3", "M3V_blur"):
    num = 0.0; den = 0.0
    for s in SECTORS:
        Ci = np.linalg.inv(cov[s]); rv = resp[s][h][0]
        num += rv @ Ci @ obs[s]; den += rv @ Ci @ rv
    X = num / den; Xe = 1 / np.sqrt(den)
    c2 = float(sum((obs[s] - X * resp[s][h][0]) @ np.linalg.inv(cov[s]) @ (obs[s] - X * resp[s][h][0]) for s in SECTORS))
    # T-only (3 observables) variant
    num3 = 0.0; den3 = 0.0
    for s in SECTORS:
        Ci = np.linalg.inv(cov[s][1:, 1:]); rv = resp[s][h][0][1:]
        num3 += rv @ Ci @ obs[s][1:]; den3 += rv @ Ci @ rv
    X3 = num3 / den3
    c23 = float(sum((obs[s][1:] - X3 * resp[s][h][0][1:]) @ np.linalg.inv(cov[s][1:, 1:]) @ (obs[s][1:] - X3 * resp[s][h][0][1:]) for s in SECTORS))
    host[h] = dict(X=float(X), X_err=float(Xe), chi2=c2, dof=23, p=float(C2.sf(c2, 23)), X_Tonly=float(X3), chi2_Tonly=c23, dof_Tonly=17, p_Tonly=float(C2.sf(c23, 17)))
    print(f"host {h:8s}: flux loss on host {X:.1f}+-{Xe:.1f} e/s; chi2 {c2:.1f}/23 (p {C2.sf(c2,23):.3g}); three-star-only chi2 {c23:.1f}/17 (p {C2.sf(c23,17):.3g})", flush=True)
null_c2 = float(sum(obs[s] @ np.linalg.inv(cov[s]) @ obs[s] for s in SECTORS))
mis = {}
for truth in ("M3V_to_bg3", "M3V_away_bg3", "M3V_blur", "bg3"):
    row = {}
    for h in ("M3V", "bg3", "pair"):
        num = sum(resp[s][h][0] @ np.linalg.inv(cov[s]) @ (15 * resp[s][truth][0]) for s in SECTORS)
        den = sum(resp[s][h][0] @ np.linalg.inv(cov[s]) @ resp[s][h][0] for s in SECTORS); X = num / den
        row[h] = float(sum((15 * resp[s][truth][0] - X * resp[s][h][0]) @ np.linalg.inv(cov[s]) @ (15 * resp[s][truth][0] - X * resp[s][h][0]) for s in SECTORS))
    mis[truth] = row
    print(f"expected noise-free chi2 if the truth is {truth:12s}: M3V {row['M3V']:.1f}, bg3 {row['bg3']:.1f}, pair {row['pair']:.1f}  (observed bg3-M3V: {host['bg3']['chi2']-host['M3V']['chi2']:.1f})")
save_json(dict(NF=NF, nreal=nreal, obs={str(s): obs[s] for s in SECTORS}, cov={str(s): cov[s] for s in SECTORS},
               response={str(s): {h: dict(mean=v[0], std=v[1]) for h, v in resp[s].items()} for s in SECTORS},
               host=host, chi2_no_signal=null_c2, expected_chi2_mismatch=mis), f"{R}/out/crosstalk.json")
print(f"no-signal chi2 {null_c2:.1f}/24")
