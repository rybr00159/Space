"""rwd: independent (non-batman) model of a small uniform hot source (WD) totally/partially occulted by a large
dark disk (the M dwarf), with exposure integration; secondary eclipse (WD transiting the M dwarf) with LD and
microlensing; and Kepler geometry. Used for the white-dwarf false-positive scenario."""
import numpy as np, batman
from rlib import G, MSUN, RSUN, DAY, C, SUB, r_nauenberg, kepler_a_rsun

def overlap_frac(d, k):
    """fraction of a uniform unit disk (source, radius 1) covered by a disk of radius k at centre separation d"""
    d = np.asarray(d, float); out = np.zeros_like(d)
    full = d <= k - 1                      # source fully covered (k > 1)
    out[full] = 1.0
    inside = d <= 1 - k                    # occulter fully inside source (k < 1)
    out[inside] = k**2
    part = (~full) & (~inside) & (d < 1 + k)
    dd = d[part]
    a1 = np.arccos(np.clip((dd**2 + 1 - k**2)/(2*dd), -1, 1))
    a2 = np.arccos(np.clip((dd**2 + k**2 - 1)/(2*dd*k), -1, 1))
    area = a1 + k**2*a2 - 0.5*np.sqrt(np.clip((-dd + 1 + k)*(dd + 1 - k)*(dd - 1 + k)*(dd + 1 + k), 0, None))
    out[part] = area/np.pi
    return out

class Super:
    """supersampled time grid for points within `near` of a transit"""
    def __init__(self, lc):
        self.lc = lc
        idx = np.where(lc.near)[0]
        ts, own = [], []
        for c in np.unique(lc.cad[idx]):
            ii = idx[lc.cad[idx] == c]; ns = SUB[int(c)]
            off = ((np.arange(ns) + 0.5)/ns - 0.5)*c/DAY
            ts.append((lc.t[ii][:, None] + off[None, :]).ravel()); own.append(np.repeat(ii, ns))
        self.ts = np.concatenate(ts); self.own = np.concatenate(own); self.idx = idx
        self.cnt = np.bincount(self.own, minlength=lc.n)
    def average(self, vals):
        s = np.bincount(self.own, vals, minlength=self.lc.n)
        out = np.ones(self.lc.n); m = self.cnt > 0
        out[m] = s[m]/self.cnt[m]
        return out

def wd_occultation(sup, T0, P, a_rsun, inc_cos, RA, RB, fA):
    """flux (normalised to 1 out of eclipse) when a uniform source of radius RA supplying fraction fA of the light
    is occulted by a dark disk of radius RB on a circular orbit of semi-major axis a (Rsun), cos i = inc_cos"""
    ph = 2*np.pi*(sup.ts - T0)/P
    x = a_rsun*np.sin(ph); y = a_rsun*inc_cos*np.cos(ph)
    front = np.cos(ph) > 0         # occulter in front of the source at phase 0
    d = np.sqrt(x**2 + y**2)/RA
    fr = np.where(front, overlap_frac(d, RB/RA), 0.)
    return sup.average(1 - fA*fr)

def secondary_depth(MA, RA, MB, RB, P, bB_in_RB, fA, u=(0.125, 0.533)):
    """depth (fraction of total light) when the WD (radius RA, mass MA) transits the M dwarf (RB) at phase 0.5:
    occultation of the LD M dwarf minus the microlensing excess 2 (R_E/R_B)^2 x local intensity (Agol 2003),
    diluted by 1 - fA. bB_in_RB = impact parameter in units of RB."""
    a = kepler_a_rsun(MA + MB, P)
    RE = np.sqrt(4*G*MA*MSUN*a*RSUN/C**2)/RSUN
    if bB_in_RB >= 1 + RA/RB: return 0., RE
    mu = np.sqrt(max(0., 1 - min(bB_in_RB, 1)**2))
    u1, u2 = u
    Iloc = (1 - u1*(1 - mu) - u2*(1 - mu)**2)/(1 - u1/3 - u2/6)
    p = batman.TransitParams(); p.t0 = 0; p.per = P; p.rp = RA/RB; p.a = a/RB; p.inc = np.degrees(np.arccos(min(1, bB_in_RB*RB/a)))
    p.ecc = 0; p.w = 90; p.limb_dark = 'quadratic'; p.u = list(u)
    occ = 1 - batman.TransitModel(p, np.array([0.0])).light_curve(p)[0]
    lens = 2*(RE/RB)**2*Iloc if bB_in_RB < 1 else 0.
    return (1 - fA)*(occ - lens), RE
