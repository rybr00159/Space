"""TESS-band surface-brightness ratio <-> temperature ratio.
(1) Blackbody convolved with the TESS response (SVO TESS/TESS.Red, photon-counting weighting lambda*S(lambda)).
(2) Empirical dwarf sequence (ms_smooth.SB_T: 10^(-0.4 M_T)/R^2 from Pecaut & Mamajek, M_T from Gaia G, Bp-Rp).
J is the disc-averaged surface-brightness ratio star2/star1 in the TESS band."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import ms_smooth as MS
_d = np.loadtxt('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/ext_dl/tess_svo.dat')
LAM = _d[:, 0] * 1e-10; RESP = _d[:, 1]
HC_K = 1.438777e-2
def bb_band(T):
    T = np.atleast_1d(np.asarray(T, float))
    B = 1 / LAM ** 5 / np.expm1(HC_K / (LAM[None, :] * T[:, None]))
    return np.trapezoid(B * RESP * LAM, LAM, axis=1)
_Tg = np.linspace(1500, 8000, 1301); _Bg = bb_band(_Tg)
def T2_from_J_bb(J, T1):
    """T2 such that S_bb(T2)/S_bb(T1) = J"""
    J = np.asarray(J, float); T1 = np.asarray(T1, float)
    return np.interp(J * np.interp(T1, _Tg, _Bg), _Bg, _Tg)
_mg = np.exp(np.linspace(np.log(0.08), np.log(1.25), 4000)); _sb = MS.SB_T(_mg)
def M2_from_J_ms(J, M1):
    J = np.asarray(J, float); M1 = np.asarray(M1, float)
    return np.interp(J * MS.SB_T(M1), _sb, _mg)
def T2T1_from_J_ms(J, M1):
    m2 = M2_from_J_ms(J, M1); return MS.teff(m2) / MS.teff(M1)
if __name__ == '__main__':
    for T1 in (3300, 3600, 3900, 4200, 4600, 5000):
        print(T1, [round(float(T2_from_J_bb(j, T1)) / T1, 3) for j in (0.4, 0.45, 0.5, 0.55, 0.6)])
    for M1 in (0.35, 0.45, 0.55, 0.65, 0.8):
        print(M1, round(float(MS.teff(M1))), [round(float(T2T1_from_J_ms(j, M1)), 3) for j in (0.4, 0.45, 0.5, 0.55, 0.6)],
              [round(float(M2_from_J_ms(j, M1)), 3) for j in (0.4, 0.5, 0.6)])
