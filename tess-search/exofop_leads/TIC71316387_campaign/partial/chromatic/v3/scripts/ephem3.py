"""Linear ephemeris from the per-sector TESS eclipse times of tess3.py, and its propagation to the
ground-survey epochs.  Errors per sector: max(formal, event-bootstrap); the fit error is scaled by
sqrt(chi2/dof) when > 1.  Deep and shallow fitted separately and jointly (common P, separate zero points)."""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L

T = json.load(open(L.V3 + '/results/tess3.json'))
shp = T['shapes']
rows = []
for sec, r in T['times'].items():
    for w in 'ds':
        off = (shp[w]['phc'] - (0 if w == 'd' else 0.5)) * L.P * 1440   # template centre offset (min)
        e = max(r[w]['err_formal_min'], r[w]['err_boot_min'])
        rows.append((int(sec), w, r[w]['E'] + (0 if w == 'd' else 0.5), r[w]['dt_min'] + off, e))
rows = np.array(rows, dtype=object)
out = {}
print('sector  type  E        O-C(min, vs nominal)  err')
for r in rows:
    print('%4d   %s  %8.1f   %+6.2f   %.2f' % tuple(r))


def linfit(E, y, e):
    A = np.vstack([np.ones_like(E), E]).T
    W = 1 / e ** 2
    C = np.linalg.inv(A.T @ (A * W[:, None]))
    p = C @ (A.T @ (W * y))
    chi = float(np.sum(W * (y - A @ p) ** 2)); dof = len(y) - 2
    s = np.sqrt(max(1, chi / dof))
    return p, C * s ** 2, chi, dof


for nm, sel in (('deep', rows[:, 1] == 'd'), ('shallow', rows[:, 1] == 's')):
    E = rows[sel, 2].astype(float); y = rows[sel, 3].astype(float); e = rows[sel, 4].astype(float)
    p, C, chi, dof = linfit(E, y, e)
    dP = p[1] / 1440; sP = np.sqrt(C[1, 1]) / 1440
    out[nm] = dict(oc0_min=float(p[0]), oc0_err=float(np.sqrt(C[0, 0])), dP_d=float(dP), dP_err_d=float(sP), chi2=chi, dof=dof,
                   P=float(L.P + dP), P_err=float(sP), C=C.tolist(), p=p.tolist())
    print('%-8s P = %.8f +- %.8f d (dP = %+.2e), O-C at E=0: %+.2f +- %.2f min, chi2 %.1f/%d' % (
        nm, L.P + dP, sP, dP, p[0], np.sqrt(C[0, 0]), chi, dof))
# joint: common slope, separate intercepts
E = rows[:, 2].astype(float); y = rows[:, 3].astype(float); e = rows[:, 4].astype(float); isd = (rows[:, 1] == 'd').astype(float)
A = np.vstack([isd, 1 - isd, E]).T
W = 1 / e ** 2
C = np.linalg.inv(A.T @ (A * W[:, None])); p = C @ (A.T @ (W * y))
chi = float(np.sum(W * (y - A @ p) ** 2)); dof = len(y) - 3; C *= max(1, chi / dof)
out['joint'] = dict(oc0_deep=float(p[0]), oc0_shallow=float(p[1]), dP_d=float(p[2] / 1440), dP_err_d=float(np.sqrt(C[2, 2]) / 1440),
                    P=float(L.P + p[2] / 1440), chi2=chi, dof=dof, C=C.tolist(), p=p.tolist())
print('joint    P = %.8f +- %.8f d; O-C(E=0) deep %+.2f shallow %+.2f min (shallow - deep = %+.2f +- %.2f min); chi2 %.1f/%d' % (
    L.P + p[2] / 1440, np.sqrt(C[2, 2]) / 1440, p[0], p[1], p[1] - p[0], np.sqrt(C[0, 0] + C[1, 1] - 2 * C[0, 1]), chi, dof))
# propagation to survey epochs (deep eclipse), relative to the nominal ephemeris
pd_, Cd = np.array(out['deep']['p']), np.array(out['deep']['C'])
for nm, tmid in (('ASAS-3 (2001-2009)', 2453500.0), ('SuperWASP (2006-2008)', 2454300.0), ('ASAS-SN V (2014-2018)', 2457600.0),
                 ('NEOWISE (2014-2024)', 2458800.0)):
    Ec = (tmid - L.T0) / L.P
    oc = pd_[0] + pd_[1] * Ec
    err = np.sqrt(Cd[0, 0] + 2 * Ec * Cd[0, 1] + Ec ** 2 * Cd[1, 1])
    nominal_err = abs(Ec) * L.P_ERR * 1440
    out['prop_' + nm.split()[0]] = dict(E=float(Ec), oc_min=float(oc), err_min=float(err), nominal_P_err_min=float(nominal_err))
    print('%-24s E=%7.0f  predicted O-C (vs nominal) %+6.2f +- %.2f min; nominal-ephemeris 1-sigma from P_err alone %.2f min' % (
        nm, Ec, oc, err, nominal_err))
json.dump(out, open(L.V3 + '/results/ephem3.json', 'w'), indent=1)
