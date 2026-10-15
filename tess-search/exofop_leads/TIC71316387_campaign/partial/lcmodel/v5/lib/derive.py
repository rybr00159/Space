"""Distance-independent derived quantities from the eclipse parameters.
 * stellar density of star 1 from Kepler's law: rho1 = 3 pi / (G P^2) * (a/R1)^3 / (1+q)
 * MS mass/radius/Teff of star 1 from rho1 (Pecaut & Mamajek dwarf sequence, lib/ms.py)
 * T2/T1 from the TESS-band surface-brightness ratio J (blackbody x TESS response, photon weighting) and from the
   empirical MS surface-brightness sequence
 * bound-hypothesis predictions: light fraction of a MS pair at the G star's distance, and the photometric
   distance of the pair implied by the eclipse light fraction f_EB"""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import ms
G_CGS = 6.6743e-8; RHO_SUN = 1.40894   # g/cm^3 (IAU nominal)
P_D = 1.3036712
_d = np.loadtxt('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/ext/tess_svo.dat')
LAM = _d[:, 0] * 1e-10; RESP = _d[:, 1]
HC_K = 1.438777e-2
T_TOT = 11.358      # TIC-8 Tmag of the blend (includes any unresolved companion)
MU_G = 5 * np.log10(1000 / 3.3622) - 5   # G-star distance modulus (Gaia DR3 parallax, no zero-point correction)
A_T = 0.02          # adopted TESS-band extinction (E(B-V) ~ 0.01-0.02 at b = -52 deg)

def rho1(r1, q, P=P_D):
    """density of star 1 in solar units"""
    Ps = P * 86400
    return 3 * np.pi / (G_CGS * Ps ** 2) * r1 ** -3 / (1 + q) / RHO_SUN

def bb_band(T):
    T = np.atleast_1d(np.asarray(T, float))
    B = 1 / LAM ** 5 / np.expm1(HC_K / (LAM[None, :] * T[:, None]))
    return np.trapezoid(B * RESP * LAM, LAM, axis=1)
_Tg = np.linspace(1500, 9000, 3001); _Bg = bb_band(_Tg)
def T2_from_J_bb(J, T1):
    J = np.asarray(J, float); T1 = np.asarray(T1, float)
    return np.interp(J * np.interp(T1, _Tg, _Bg), _Bg, _Tg)

def pair_MT(m1, m2):
    return -2.5 * np.log10(10 ** (-0.4 * ms.MT(m1)) + 10 ** (-0.4 * ms.MT(m2)))

def fEB_if_bound(m1, m2, mu=MU_G, a_t=A_T, ttot=T_TOT):
    """light fraction of a MS pair (m1, m2) placed at distance modulus mu, relative to the observed blend"""
    return 10 ** (-0.4 * (pair_MT(m1, m2) + mu + a_t - ttot))

def distance_from_fEB(fEB, m1, m2, a_t=A_T, ttot=T_TOT):
    """photometric distance (pc) of a MS pair contributing fraction fEB of the blend's TESS flux"""
    mT = ttot - 2.5 * np.log10(fEB)
    return 10 ** ((mT - pair_MT(m1, m2) - a_t + 5) / 5)
