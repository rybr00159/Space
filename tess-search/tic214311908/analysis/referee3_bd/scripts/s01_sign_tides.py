"""s01: (a) independent check of the photo-eccentric relation used by fpp/f04 and referee2_fpp/r05 (claimed sign error):
numerical T14 from my own Kepler-orbit sky-projection code (no batman) for several (e, w) -> rho_circ/rho_true vs g^+-3;
(b) transit-probability-weighted tail P(rho_c/rho* <= x) for several e priors (own draws);
(c) tidal eccentricity damping for a BD at P = 5.6915 d around the 0.37 Msun M3V: Jackson+08 small-e constant-Q rate and
the Hut (1981) eccentricity enhancement f(e) at e = 0.3-0.6 (constant time lag; ratio to the e->0 rate), star spin
Omega/n = 0 and 0.57 (P_rot 10 d).
Output: out/s01_sign_tides.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
from scipy.optimize import brentq
out = {}
# ---------- (a) own orbit code: sky separation (in R*) vs time; T14 = time with separation < 1 + k
def kepler_E(M, e):
    E = M.copy()
    for _ in range(60): E = E - (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
    return E
def T14_num(e, w_deg, aR=25.4, b_circ=0.0, k=0.12, P=PC):
    """observer along +z; planet transits when its true anomaly f = 90 deg - w. inclination chosen so that the impact
    parameter at conjunction is b (in R*)"""
    w = np.radians(w_deg)
    rt = aR*(1 - e**2)/(1 + e*np.sin(w))            # star-planet distance at conjunction (R*)
    cosi = b_circ/rt
    t = np.linspace(-0.4, 0.4, 400001)              # days around conjunction
    # time of conjunction: true anomaly fc = pi/2 - w -> mean anomaly
    fc = np.pi/2 - w
    Ec = 2*np.arctan2(np.sqrt(1 - e)*np.sin(fc/2), np.sqrt(1 + e)*np.cos(fc/2))
    Mc = Ec - e*np.sin(Ec)
    M = Mc + 2*np.pi*t/P
    E = kepler_E(M, e)
    f = 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
    r = aR*(1 - e*np.cos(E))
    # orbit-plane coords with periapsis along x; rotate by w; sky: X = -r cos(w+f), Y = -r sin(w+f) cos i, Z = r sin(w+f) sin i
    X = -r*np.cos(w + f); Y = -r*np.sin(w + f)*cosi; Z = r*np.sin(w + f)
    d = np.hypot(X, Y)
    inn = (d < 1 + k) & (Z > 0)
    return (t[inn][-1] - t[inn][0])*24
chk = []
Tc = T14_num(0.0, 90.)
for e, wd in [(0.3, 90.), (0.3, 270.), (0.5, 90.), (0.5, 270.), (0.6, 250.), (0.4, 30.), (0.55, 300.)]:
    g = (1 + e*np.sin(np.radians(wd)))/np.sqrt(1 - e**2)
    T = T14_num(e, wd)
    chk.append(dict(e=e, w_deg=wd, T14_h=T, T14_circ_h=Tc, rho_c_over_rho=(Tc/T)**3, g_plus3=g**3, g_minus3=g**-3,
                    ptr_factor=(1 + e*np.sin(np.radians(wd)))/(1 - e**2)))
    print(chk[-1], flush=True)
out['sign_check'] = dict(rows=chk, verdict='rho_c/rho = g^+3 (matches to <1%); r05/f04 use g^-3 -> sign error confirmed' if all(abs(c['rho_c_over_rho']/c['g_plus3'] - 1) < 0.02 for c in chk) else 'MISMATCH')
# ---------- (b) transit-weighted tails
rng = np.random.default_rng(3)
N = 2_000_000
w = rng.uniform(0, 2*np.pi, N)
pri = {'kipping13 Beta(0.867,3.03)': rng.beta(0.867, 3.03, N),
       'vaneylen19 half-normal 0.32 (singles)': np.abs(rng.normal(0, 0.32, N)),
       'BD short-period Beta(0.426,2.354) (Ferreira+26, central)': rng.beta(0.426, 2.354, N),
       'BD long-period Beta(1.879,2.470) (Ferreira+26, central)': rng.beta(1.879, 2.470, N),
       'uniform 0-0.9': rng.uniform(0, 0.9, N)}
tails = {}
for k, e in pri.items():
    e = np.clip(e, 0, 0.9)
    g = (1 + e*np.sin(w))/np.sqrt(1 - e**2); geo = (1 + e*np.sin(w))/(1 - e**2)
    wt = geo/geo.mean()
    r = g**3
    tails[k] = {f'P_tw(rho_c/rho<={x})': float(np.mean(wt*(r <= x))) for x in (0.1, 0.15, 0.2, 0.25, 0.3, 0.5)}
    tails[k]['P_tw(0.12<=rho_c/rho<=0.30)'] = float(np.mean(wt*(r >= 0.12)*(r <= 0.30)))
    tails[k]['mean_ptr_boost'] = float(geo.mean())
    # wrong sign (r05/f04): rho_c = rho g^-3 with the same weights
    tails[k]['WRONGSIGN_P_tw(rho_c/rho<=0.2)'] = float(np.mean(wt*(g**-3 <= 0.2)))
    print(k, tails[k], flush=True)
out['rho_ratio_tails'] = tails
# ---------- (c) tides
Ms, Rs = MSTAR*MSUN, RSTAR*RSUN
def a_of(Mbd_MJ): return (G*(Ms + Mbd_MJ*MJUP)*(PC*DAY)**2/(4*np.pi**2))**(1/3.)
def f3(e): return 1 + 15/4*e**2 + 15/8*e**4 + 5/64*e**6
def f4(e): return 1 + 3/2*e**2 + 1/8*e**4
def hut_ratio(e, om):    # (de/dt / e) at e relative to e->0, constant time lag (Hut 1981 eq. 10)
    num = (f3(e) - 11/18*(1 - e**2)**1.5*f4(e)*om)/(1 - e**2)**6.5
    den = 1 - 11/18*om
    return num/den
T = []
for Mbd in (13, 30, 50, 65, 80):
    a = a_of(Mbd)
    for Rbd in (0.8, 1.0):
        for lQs in (5, 6, 7, 8):
            for lQb in (5, 6, 7, 8):
                Mb, Rb = Mbd*MJUP, Rbd*RJUP
                rb = 63/4*np.sqrt(G*Ms**3)*Rb**5/(10**lQb*Mb)*a**-6.5
                rs = 171/16*np.sqrt(G/Ms)*Rs**5*Mb/(10**lQs)*a**-6.5
                tau0 = 1/(rb + rs)/YR/1e9
                row = dict(M_bd=Mbd, R_bd=Rbd, logQs=lQs, logQbd=lQb, a_au=a/AU, tau_e_small_e_Gyr=tau0, tau_star_only_Gyr=1/rs/YR/1e9, tau_bd_only_Gyr=1/rb/YR/1e9)
                for e in (0.3, 0.5, 0.6):
                    for om in (0.0, 0.57):
                        row[f'tau_e_at_e{e}_om{om}_Gyr'] = tau0/hut_ratio(e, om)
                T.append(row)
import pandas as pd
tdf = pd.DataFrame(T)
sel = tdf[(tdf.M_bd == 65) & (tdf.R_bd == 0.8) & (tdf.logQbd == 7)]
print(sel[['logQs', 'a_au', 'tau_e_small_e_Gyr', 'tau_e_at_e0.5_om0.0_Gyr', 'tau_e_at_e0.5_om0.57_Gyr', 'tau_e_at_e0.6_om0.0_Gyr']].to_string(), flush=True)
# integrate de/dt = -e/tau0 * hut_ratio(e, om) backward: time to damp from e0 -> 0.1 and the e reached from 0.6 after t
def e_after(e0, tau0, t, om=0.0, n=20000):
    e = e0; dt = t/n
    for _ in range(n):
        e = max(e - dt*e/tau0*hut_ratio(e, om), 0.0)
    return e
evo = []
for lQs in (5, 6, 7):
    tau0 = float(sel[sel.logQs == lQs].tau_e_small_e_Gyr.iloc[0])
    for e0 in (0.55, 0.7):
        evo.append(dict(logQs=lQs, tau0_Gyr=tau0, e0=e0, e_after_1Gyr=e_after(e0, tau0, 1.0), e_after_3Gyr=e_after(e0, tau0, 3.0), e_after_5Gyr=e_after(e0, tau0, 5.0)))
        print(evo[-1], flush=True)
out['tides'] = dict(formula='Jackson+08 eq.2 small-e rate x Hut (1981) eccentricity enhancement (1-e^2)^-6.5 [f3 - 11/18 (1-e^2)^1.5 f4 Omega/n]/(1 - 11/18 Omega/n)',
                    hut_ratio={f'e={e}': {f'om={om}': hut_ratio(e, om) for om in (0.0, 0.57)} for e in (0.3, 0.5, 0.6)},
                    table=T, evolution_65MJ_0p8RJ_Qbd1e7=evo,
                    bd_over_star_tide_ratio_65MJ={f'Qbd/Qs=1': float(sel[sel.logQs == 7].tau_star_only_Gyr.iloc[0]/sel[sel.logQs == 7].tau_bd_only_Gyr.iloc[0])})
jdump(dict(method=__doc__, **out), f'{R3}/out/s01_sign_tides.json')
