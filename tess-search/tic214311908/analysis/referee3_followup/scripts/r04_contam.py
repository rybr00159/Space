"""r04: aperture contamination of the M3V by every Gaia source within 25" (own code).
Magnitudes: M3V and pair A from Gaia XP synthetic (data/gaia_gspc.csv); other stars from DECaPS2 (../star/data/decaps2_90.csv,
1" match after propagating Gaia positions to the DECaPS2 epoch) or, if absent, from Gaia G/BP/RP with a colour-dependent
offset fitted on the matched stars of this field (own fit). Mira: brightest/faintest band magnitudes from the analyst's
ranges are not used; own conservative values: g 13.3, r 12.4, i 10.2, z 9.2 at maximum (bright end of SkyMapper/Gaia
ranges, rounded brighter).
Positions at 2027.5 (Gaia PM). PSF: circular Moffat, beta 2.5 / 4, FWHM 1.0-2.5"; aperture radius 1.5-4"; enclosed flux by
direct 2-D integration on a 0.02" grid (independent of the analyst's polar integration).
Outputs: out/r04_contam.csv, out/r04_neighbours.csv"""
import os, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd
B='tic214311908/analysis'
R=f'{B}/referee3_followup'
g=pd.read_csv(f'{R}/data/gaia_cone25.csv'); s=pd.read_csv(f'{R}/data/gaia_gspc.csv')
dc=pd.read_csv(f'{B}/star/data/decaps2_90.csv')
t0=g[g.source_id==5925373272850784768].iloc[0]
def pos(row, ep):
    dt=ep-2016.0; pmra=row.pmra if np.isfinite(row.pmra) else 0; pmdec=row.pmdec if np.isfinite(row.pmdec) else 0
    return row.ra+pmra*dt/3.6e6/np.cos(np.radians(row.dec)), row.dec+pmdec*dt/3.6e6
rows=[]
for r in g.itertuples():
    ra27,de27=pos(r,2027.5); ra0,de0=pos(t0,2027.5)
    dx=(ra27-ra0)*np.cos(np.radians(de0))*3600; dy=(de27-de0)*3600
    rad,ded=pos(r,2018.6)
    d=np.hypot((dc.ra-rad)*np.cos(np.radians(ded))*3600,(dc.dec-ded)*3600); j=int(np.argmin(d))
    m=dict(source_id=r.source_id, G=r.phot_g_mean_mag, bp_rp=r.bp_rp, sep27=np.hypot(dx,dy), pa27=np.degrees(np.arctan2(dx,dy))%360, dx=dx, dy=dy, dmatch=d[j])
    for b in 'griz':
        v=dc[f'mean_mag_{b}'].iloc[j]; m[f'{b}_dec']=v if (d[j]<1.0 and np.isfinite(v)) else np.nan
    if r.source_id in set(s.source_id):
        q=s[s.source_id==r.source_id].iloc[0]
        for b in 'griz': m[f'{b}_xp']=q[f'{b}_sdss_mag']
    rows.append(m)
N=pd.DataFrame(rows).sort_values('sep27')
# colour fit G - band vs bp_rp from DECaPS-matched stars (excluding saturated/very bright)
fit={}
ok=N[N.G.between(15.5,20.5)&N.bp_rp.notna()]
for b in 'griz':
    k=ok[f'{b}_dec'].notna()
    x=ok.bp_rp[k].values; y=(ok[f'{b}_dec'][k]-ok.G[k]).values
    if len(x)>4:
        c=np.polyfit(x,y,2); fit[b]=c
        print(b,'G->band fit n',len(x),'rms',np.std(y-np.polyval(c,x)).round(3))
for b in 'griz':
    N[b]=N[f'{b}_xp'] if f'{b}_xp' in N else np.nan
    N[b]=N[b].fillna(N[f'{b}_dec'])
    est=N.G+np.polyval(fit[b],N.bp_rp.fillna(1.2))
    N[f'{b}_src']=np.where(N[f'{b}_xp'].notna(),'XP',np.where(N[f'{b}_dec'].notna(),'DECaPS2','G+fit'))
    N[b]=N[b].fillna(est)
MIRA=5925373268528684288
mira_max=dict(g=13.3,r=12.4,i=10.2,z=9.2); mira_min=dict(g=17.1,r=15.6,i=12.1,z=10.6)
N.to_csv(f'{R}/out/r04_neighbours.csv',index=False)
print(N[['source_id','sep27','pa27','G','bp_rp','g','r','i','z','g_src','dmatch']].head(14).to_string())
M=N[N.source_id==5925373272850784768].iloc[0]
def moffat_grid(fwhm,beta):
    a=fwhm/(2*np.sqrt(2**(1/beta)-1)); return a
xs=np.arange(-6,6,0.02)+0.01; X,Y=np.meshgrid(xs,xs)
def frac_in(d, rap, fwhm, beta):
    a=moffat_grid(fwhm,beta)
    # star at (d,0); aperture at origin; integrate on the grid within rap
    inap=(X**2+Y**2)<=rap**2
    I=(beta-1)/(np.pi*a*a)*(1+((X-d)**2+Y**2)/(a*a))**(-beta)
    return float(I[inap].sum()*0.02**2)
out=[]
for beta in (2.5,4.0):
    for fwhm in (1.0,1.5,2.0,2.5):
        for rap in (1.5,2.0,3.0,4.0):
            f0=frac_in(0,rap,fwhm,beta)
            row=dict(beta=beta,fwhm=fwhm,rap=rap,enc_M3V=f0)
            for b in 'griz':
                tot=0; parts={}
                for n in N.itertuples():
                    if n.source_id==5925373272850784768 or n.sep27>24: continue
                    if n.source_id==MIRA:
                        fr=frac_in(n.sep27,rap,fwhm,beta) if n.sep27<6+rap else None
                        if fr is None:   # analytic far-wing: Moffat surface brightness at d times aperture area
                            a=moffat_grid(fwhm,beta); fr=(beta-1)/(np.pi*a*a)*(1+(n.sep27/a)**2)**(-beta)*np.pi*rap**2
                        for lab,mm in (('mira_max',mira_max[b]),('mira_min',mira_min[b])):
                            row[f'c_{b}_{lab}']=10**(-0.4*(mm-getattr(M,b)))*fr/f0
                        continue
                    if n.sep27<6+rap: fr=frac_in(n.sep27,rap,fwhm,beta)
                    else:
                        a=moffat_grid(fwhm,beta); fr=(beta-1)/(np.pi*a*a)*(1+(n.sep27/a)**2)**(-beta)*np.pi*rap**2
                    c=10**(-0.4*(getattr(n,b)-getattr(M,b)))*fr/f0
                    tot+=c; parts[n.sep27]=c
                row[f'c_{b}_stars']=tot
                s31=[v for k,v in parts.items() if 3.0<k<3.4]; row[f'c_{b}_3p1']=s31[0] if s31 else np.nan
                s43=[v for k,v in parts.items() if 4.2<k<4.5]; row[f'c_{b}_4p3']=s43[0] if s43 else np.nan
            out.append(row)
C=pd.DataFrame(out); C.to_csv(f'{R}/out/r04_contam.csv',index=False)
pd.set_option('display.width',250)
show=['beta','fwhm','rap']+[f'c_{b}_{t}' for b in 'gz' for t in ('stars','3p1','4p3','mira_max')]
print(C[show].round(5).to_string())
C['dg_z']=(C.c_g_stars-C.c_z_stars)
print('differential g-z dilution (stars only): ', C[['beta','fwhm','rap','dg_z']].round(4).to_string())
