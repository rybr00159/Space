"""Derived quantities along the geometric L3 profile (best solution at each fixed L3, either branch):
 rho1 from Kepler: rho1 = 3 pi / (G P^2) / r1^3 / (1+q), with q from the main-sequence mass-radius relation implied by k
 (q such that R(M1 q)/R(M1) = k at the M1 implied by rho1 -- solved iteratively with the smooth MS relation and fR=1.04),
 and alternatively for q = 0.2 and q = 1 (bracketing);
 predicted ellipsoidal semi-amplitude in total light (Morris & Naftilan 1993 leading term, alpha = 0.15 (15+u)(1+tau)/(3-u),
 u = 0.5, tau = 0.4 -> alpha = 1.30):  A_ell = (1-L3) alpha sin^2 i [ L1 q r1^3 + L2 r2^3 / q ];
 compared with the measured OOE limit (s05_ooe.py)."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import ms_smooth as MS, lclib as L
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final'
tag = sys.argv[1] if len(sys.argv) > 1 else 'default_recent'
res = json.load(open(f'{FD}/results/ml_{tag}.json'))
ooe = json.load(open(f'{FD}/results/ooe.json'))['default']
G = 6.674e-8; P = L.P0 * 86400; RHO_SUN = 1.408
K = 3 * np.pi / (G * P ** 2)
ALPHA = 0.15 * (15 + 0.5) * (1 + 0.4) / (3 - 0.5)
def q_from(r1, k, fR=1.04):
    q = max(0.1, k)
    for _ in range(50):
        rho = K / r1 ** 3 / (1 + q) / RHO_SUN
        M1 = MS.mass_of_rho(rho, fR)
        mg = np.linspace(0.08, M1, 400); kk = MS.radius(mg) / MS.radius(M1)
        qn = float(np.interp(k, kk, mg / M1)) if k < 1 else 1.0
        if abs(qn - q) < 1e-5: break
        q = 0.5 * (q + qn)
    return q, M1, rho
def aell(p, q):
    r1 = p['r1']; r2 = p['rsum'] - r1; k2J = p['k'] ** 2 * p['J']; L2 = k2J / (1 + k2J); L1 = 1 - L2
    s2 = 1 - p['cosi'] ** 2
    return (1 - p['L3']) * ALPHA * s2 * (L1 * q * r1 ** 3 + L2 * r2 ** 3 / q)
gb = min(r['tot'] for r in res)
out = []
print(f"OOE limit on A_ell (default): joint {ooe['Aell_joint']:+.1f} +- {ooe['Aell_joint_err']:.1f} ppm (UL95 {ooe['Aell_joint_UL95']:.0f}); "
      f"sector-mean {ooe['Aell_sectors']:+.1f} +- {ooe['Aell_sectors_err']:.1f} (UL95 {ooe['Aell_sectors_UL95']:.0f})")
grid = sorted(set(r['fixed']['L3'] for r in res if r['kind'] == 'L3prof'))
for v in grid:
    rr = [r for r in res if r['kind'] == 'L3prof' and abs(r['fixed']['L3'] - v) < 1e-9]
    b = min(rr, key=lambda r: r['tot'])
    q, M1, rho = q_from(b['r1'], b['k'])
    rho_q02 = K / b['r1'] ** 3 / 1.2 / RHO_SUN; rho_q1 = K / b['r1'] ** 3 / 2.0 / RHO_SUN
    A = aell(b, q) * 1e6
    z = (A - ooe['Aell_joint']) / ooe['Aell_joint_err']; z2 = (A - ooe['Aell_sectors']) / ooe['Aell_sectors_err']
    o = dict(L3=v, dchi2_lc=b['tot'] - gb, branch=b['branch'], r1=b['r1'], k=b['k'], b1=b['b1'], J=b['J'], q_MS=q, M1_MS=M1, rho1_sun=rho,
             rho1_sun_q02=rho_q02, rho1_sun_q1=rho_q1, Aell_pred_ppm=A, z_joint=z, z_sectors=z2, dchi2_ell_joint=z ** 2 if z > 0 else 0.0,
             dchi2_ell_sectors=z2 ** 2 if z2 > 0 else 0.0)
    out.append(o)
    print(f"L3={v:.3f} ({b['branch']:7s}) dchi2_LC={o['dchi2_lc']:6.2f} r1={b['r1']:.4f} k={b['k']:.3f} b1={b['b1']:.2f} J={b['J']:.3f} | q_MS={q:.2f} "
          f"rho1={rho:.2f} rho_sun (q=0.2: {rho_q02:.2f}, q=1: {rho_q1:.2f}) -> M1_MS~{M1:.2f} ({MS.spt_of_mass(M1)}) | A_ell pred {A:5.0f} ppm "
          f"(z={z:+.1f} joint, {z2:+.1f} sector-mean)")
json.dump(out, open(f'{FD}/results/profile_derived_{tag}.json', 'w'), indent=1)
