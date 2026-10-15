"""Main-sequence-constrained parametrisation of the EB fit.
Free: M1, q = M2/M1, dR1, dR2 (fractional radius offsets from the MS radius, prior N(0, SIG_R)), dJ (fractional
offset of the TESS surface-brightness ratio from the MS sequence, prior N(0, SIG_J)), b0, L3 (or, in 'bound' mode,
dmag = offset of the pair's absolute TESS magnitude from the MS sequence, prior N(0, SIG_MAG), with
1 - L3 = F_EB/F_blend fixed by placing the pair at the G star's distance), ecosw, esinw, LD, T0, P, depth scales.
Kepler: a = 4.2081 Rsun (M1+M2)^(1/3) P[d]^(2/3)."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import ms, derive, fit as Fi
SIG_R, SIG_J, SIG_MAG = 0.07, 0.10, 0.15
A1D = 4.2081

def geom(M1, q, dR1, dR2, dJ, P):
    M2 = q * M1
    R1 = ms.radius(M1) * (1 + dR1); R2 = ms.radius(M2) * (1 + dR2)
    a = A1D * (M1 + M2) ** (1 / 3) * P ** (2 / 3)
    J = ms.SB_T(M2) / ms.SB_T(M1) * (1 + dJ)
    return (R1 + R2) / a, R2 / R1, J, a

class MSPosterior:
    def __init__(self, model, fixed, free, mode='free', mmin=0.12, sig_r=SIG_R, sig_j=SIG_J, sig_mag=SIG_MAG):
        self.m = model; self.fixed = dict(fixed); self.free = list(free); self.mode = mode; self.mmin = mmin
        self.sr, self.sj, self.sm = sig_r, sig_j, sig_mag
    def full(self, x):
        p = dict(self.fixed); p.update(zip(self.free, x)); return p
    def to_geom(self, p):
        rsum, k, J, a = geom(p['M1'], p['q'], p['dR1'], p['dR2'], p['dJ'], p['P'])
        g = dict(p); g.update(rsum=rsum, k=k, J=J, a_rsun=a)
        if self.mode == 'bound':
            fe = derive.fEB_if_bound(p['M1'], p['q'] * p['M1']) * 10 ** (-0.4 * p['dmag'])
            g['L3'] = 1 - fe
        return g
    def __call__(self, x):
        p = self.full(x)
        if not (self.mmin < p['M1'] < ms.MMAX and 0.15 < p['q'] <= 1.0 and p['q'] * p['M1'] > ms.MMIN and abs(p['dR1']) < 0.3
                and abs(p['dR2']) < 0.3 and abs(p['dJ']) < 0.4):
            return -np.inf
        g = self.to_geom(p)
        if self.mode == 'bound' and not (0 < g['L3'] < 0.999): return -np.inf
        lp = Fi.log_prior(g)
        if not np.isfinite(lp): return -np.inf
        lp += -0.5 * ((p['dR1'] / self.sr) ** 2 + (p['dR2'] / self.sr) ** 2 + (p['dJ'] / self.sj) ** 2)
        if self.mode == 'bound': lp += -0.5 * (p['dmag'] / self.sm) ** 2
        return lp - 0.5 * self.m.chi2(g)
