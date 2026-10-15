import json
P='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
d=dict(status="DRAFT v3 2026-10-15 02:50 UTC: all checks done except H2 TRICERATOPS high-N run at P=1.3037 (deep-only fold), running in astro/tri2/",
summary=("SED: no IR excess (W1/W2/W3/Ks residuals within 1 sigma of 319 photometric twins). Single-star 21-band fit Teff 5800+-24 K, R 1.148+-0.011, chi2 27.4. "
 "MS pair at 297 pc improves chi2 by 9-16 (uncalibrated, 21 bands); 9-band twin null p=0.025-0.06; no pair type is excluded; 95% upper limits on red-pair light f_TESS<=5.8-11.6%. "
 "Gaia C*=+0.0204 vs twin median -0.0038 (robust sd 0.0019; 1/696 twins as high) and BP/RP blended in 47/48 of 47/48 transits: extra unresolved flux in BP/RP window. "
 "Spectroscopy: none in RAVE DR5/DR6, GALAH DR2/DR3/DR4, APOGEE DR17, Gaia-ESO DR5 (30 arcsec); SoS-I only Gaia RV. GSP-Spec MG Teff 6235 (6016-6334), ANN 5843 (5701-5998); flags_gspspec fluxNoise=1, all others 0. "
 "PMA/RV/RUWE: no anomaly; companion orbit exclusion >=90% only for a<=5-10 au; nil at the Gaia-IPD-favoured 150-360 au. "
 "Priors (arithmetic verified): bound:background 338 central, 25-8600 full range; colour-matched ~1000; IPD-annulus 147 (36-440). TRICERATOPS 1.1.0 prior ratio 158. "
 "TRICERATOPS high-N (13.2M draws, P=0.6518 fold): SEBx2P 0.973 vs BEBx2P 0.027 (68% 0.959-0.981); prior-driven (likelihood ratio ~0.22 favouring BEB)."),
key_numbers=[], conclusions=[], implications_bound_vs_background="pending final", caveats="pending final",
files=P+'astro/')
json.dump(d,open(P+'astro_result_draft.json','w'),indent=1); print('ok')
