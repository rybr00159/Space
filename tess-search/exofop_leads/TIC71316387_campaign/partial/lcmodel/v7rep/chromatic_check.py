"""Predicted primary-eclipse depth ratios (band X / TESS) for an eclipsed star of given spectral type, from the
Pecaut & Mamajek dwarf table (M_V, M_G, Bp-Rp, M_Ks, Ks-W1; TESS mag from the TIC-8 G/Bp-Rp relation), relative to
the blend V = 11.9, T = 11.358, W1 = Ks(10.341) - 0.06.  R_X = 10^(-0.4 [(X-T)_star - (X-T)_blend]).  Ignores limb-
darkening differences between bands and the exact SuperWASP passband (broader/redder than V).  Compared with the
archival measurements of the chromatic analysis: SuperWASP 0.42 +- 0.15, W1 3.7 +- 1.1."""
import sys, numpy as np, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import ms
T = ms.T
V_b, T_b, W1_b = 11.9, 11.358, 10.341 - 0.06
out = {}
for s in ['K5V', 'K7V', 'M0V', 'M1V', 'M1.5V', 'M2V', 'M2.5V', 'M3V']:
    i = np.where(T['SpT'] == s)[0][0]
    c = T['Bp-Rp'][i]; MT = T['M_G'][i] - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    MV = T['Mv'][i]; MW1 = T['M_Ks'][i] - T['Ks-W1'][i]
    RV = 10 ** (-0.4 * ((MV - MT) - (V_b - T_b))); RW = 10 ** (-0.4 * ((MW1 - MT) - (W1_b - T_b)))
    chi2 = ((RV - 0.42) / 0.15) ** 2 + ((RW - 3.7) / 1.1) ** 2
    out[s] = dict(Teff=float(T['Teff'][i]), R_V=round(float(RV), 3), R_W1=round(float(RW), 3), chi2_2pts=round(float(chi2), 2))
    print(s, out[s])
json.dump(out, open('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v7rep/chromatic_check.json', 'w'), indent=1)
