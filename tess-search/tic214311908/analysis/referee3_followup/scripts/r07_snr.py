"""r07: independent per-event transit S/N and depth-ratio precision (own noise model; assumptions differ from the analyst's).
Count rate: AB zero point 3631 Jy -> photons s^-1 m^-2 = 1.51e7 * 3631 * (dlambda/lambda); effective area pi/4 D^2 (1-eps^2);
total system efficiency (telescope+instrument+QE, excluding atmosphere) ASSUMED per band; extinction k per band at airmass 1.3.
Sky (AB mag/arcsec^2, dark, ASSUMED): g 22.0, r 21.0, i 20.2, z 19.4, I+z 19.6. Aperture r = max(1.5", 1.1 FWHM), Moffat beta 3.
Scintillation: Osborn+15 eq. 7 with C_Y 1.5 (median-ish), h per site, x1.0 (ensemble error neglected; analyst x1.2).
Per-event depth error: sigma_d = sqrt( (1.42 sigma_2min sqrt(1/N_in+1/N_out))^2 + sigma_red^2 ), N_in = 43 (1.45 h), N_out = 90.
The red term is the controlling unknown; results are shown for 0.5, 1.0 and 1.5 ppt per event.
Depth-ratio sigma for two bands: sqrt(s1^2+s2^2)/depth (red noise independent between bands; if it were fully common-mode
and achromatic the ratio error would shrink).
Output: out/r07_snr.csv, out/r07_snr.json"""
import os, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd
R='tic214311908/analysis/referee3_followup'
DEPTH=15.2e-3
mag={'g':16.36,'r':14.98,'i':13.57,'z':12.81,'I+z':12.75}       # AB; I+z ~ between i and z weighted to z
band={'g':(475,140,22.0,0.18),'r':(625,140,21.0,0.10),'i':(770,150,20.2,0.06),'z':(900,120,19.4,0.05),'I+z':(880,280,19.6,0.05)}
tel={'MuSCAT4 2m FTS':dict(D=2.0,eps=0.35,h=1165,eff={'g':0.30,'r':0.35,'i':0.30,'z':0.20},RN=12,ro=6,fwhm=1.8,pix=0.27,FW=7e4),
     'LCOGT 1m':dict(D=1.0,eps=0.40,h=2207,eff={'g':0.30,'r':0.35,'i':0.30,'z':0.15},RN=8,ro=28,fwhm=2.0,pix=0.389,FW=1e5),
     'SPECULOOS 1m':dict(D=1.0,eps=0.30,h=2635,eff={'I+z':0.40},RN=6,ro=10,fwhm=1.3,pix=0.35,FW=1e5),
     'TRAPPIST 0.6m':dict(D=0.6,eps=0.35,h=2400,eff={'I+z':0.30,'z':0.2},RN=9,ro=6,fwhm=1.6,pix=0.64,FW=1e5),
     'amateur 0.4m':dict(D=0.4,eps=0.45,h=500,eff={'i':0.30,'r':0.30},RN=3,ro=2,fwhm=2.5,pix=0.6,FW=5e4,skyoff=-1.0)}
AM=1.3; beta=3.0
def enc(r,fwhm):
    a=fwhm/(2*np.sqrt(2**(1/beta)-1)); return 1-(1+(r/a)**2)**(1-beta)
def peakfrac(fwhm,pix):
    a=fwhm/(2*np.sqrt(2**(1/beta)-1)); return (beta-1)/(np.pi*a*a)*pix*pix
rows=[]
for tn,t in tel.items():
    A=np.pi/4*t['D']**2*(1-t['eps']**2)
    for b,eff in t['eff'].items():
        c,w,sky,kx=band[b]
        n0=1.51e7*3631*(w/c)*A*eff*10**(-0.4*kx*AM)
        star=n0*10**(-0.4*mag[b]); skyr=n0*10**(-0.4*(sky+t.get('skyoff',0)))
        rap=max(1.5,1.1*t['fwhm']); e=enc(rap,t['fwhm']); npix=np.pi*rap**2/t['pix']**2
        # exposure: 70% full well for in-focus FWHM 1.3" at best, cap 120 s
        tsat=0.7*t['FW']/(star*peakfrac(1.3,t['pix'])+skyr*t['pix']**2); texp=float(np.clip(0.8*tsat,5,120))
        S=star*e*texp; sig_ph=np.sqrt(S+skyr*np.pi*rap**2*texp+npix*t['RN']**2)/S
        sig_sc=np.sqrt(1e-5*1.5**2*t['D']**(-4/3)/texp*AM**3*np.exp(-2*t['h']/8000))
        nexp=120/(texp+t['ro']); s2=np.hypot(sig_ph,sig_sc)/np.sqrt(max(nexp,1))
        white=1.42*s2*np.sqrt(1/43+1/90)
        row=dict(telescope=tn,band=b,e_per_s=round(star*e),t_exp=round(texp,1),sig_phot_exp_ppt=round(sig_ph*1e3,2),sig_scint_exp_ppt=round(sig_sc*1e3,2),
                 sig_2min_ppt=round(s2*1e3,2),white_depth_ppt=round(white*1e3,3))
        for red in (0.5,1.0,1.5):
            sd=np.hypot(white,red*1e-3); row[f'SN_red{red}']=round(DEPTH/sd,1); row[f'sdepth_red{red}']=round(sd*1e3,3)
        rows.append(row)
S=pd.DataFrame(rows); S.to_csv(f'{R}/out/r07_snr.csv',index=False)
pd.set_option('display.width',250); print(S.to_string())
# depth ratio per event (MuSCAT4 simultaneous) and events for 5% vs red noise
out={}
for red in (0.5,1.0,1.5):
    m=S[S.telescope=='MuSCAT4 2m FTS'].set_index('band')
    for p,(a,b) in {'r/z':('r','z'),'g/z':('g','z')}.items():
        sr=np.hypot(m.loc[a,f'sdepth_red{red}'],m.loc[b,f'sdepth_red{red}'])/15.2
        out[f'MuSCAT4 {p} red{red}']=dict(sigma_R_event=round(float(sr),3), events_for_5pct=round(float((sr/0.05)**2),1))
    l=S[S.telescope=='LCOGT 1m'].set_index('band')
    sr=np.hypot(l.loc['r',f'sdepth_red{red}'],l.loc['z',f'sdepth_red{red}'])/15.2
    out[f'LCOGT1m r/z alternating red{red}']=dict(sigma_R_pair=round(float(sr),3), events_for_5pct=round(float(2*(sr/0.05)**2),1))
print(json.dumps(out,indent=0))
json.dump(dict(table=rows,ratio=out),open(f'{R}/out/r07_snr.json','w'),indent=1)
