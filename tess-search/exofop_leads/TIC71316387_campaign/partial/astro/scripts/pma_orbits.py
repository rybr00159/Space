"""Monte Carlo limits on a bound companion (mass Mc, semimajor axis a) from:
 (1) long-baseline vs Gaia DR3 proper-motion difference (TGAS-like: mean velocity between Tycho-2 mean epoch and 2015.0, vs DR3 velocity averaged over 2014.6-2017.4),
 (2) Gaia DR3 RV constancy (15 transits over 945 d, chi2 p=0.275),
 (3) Gaia DR3 RUWE=1.026 (photocentre wobble within the 34-month window).
Exclusion fraction = fraction of random orbits (isotropic, e~U(0,0.9), random phase) incompatible at 95%."""
import numpy as np, json
from scipy import stats
rng=np.random.default_rng(42)
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/'
MG=1.03; plx=3.3622; d_pc=1000/plx; beta=0.03   # companion light fraction in G
mas_per_au=plx
def kepler(M,e):
    E=M.copy()
    for _ in range(40): E=E-(E-e*np.sin(E)-M)/(1-e*np.cos(E))
    return E
def orbit_xyz(t,P,e,w,Om,inc,tp,a):
    """relative orbit (companion wrt G star) in au; returns (dRA,dDec) on sky [au] and z [au]; t in yr arrays (n,) ; params (m,1)"""
    M=2*np.pi*(t-tp)/P; E=kepler(np.mod(M,2*np.pi),e*np.ones_like(M))
    X=a*(np.cos(E)-e); Y=a*np.sqrt(1-e**2)*np.sin(E)
    A=np.cos(w)*np.cos(Om)-np.sin(w)*np.sin(Om)*np.cos(inc); B=np.cos(w)*np.sin(Om)+np.sin(w)*np.cos(Om)*np.cos(inc)
    F=-np.sin(w)*np.cos(Om)-np.cos(w)*np.sin(Om)*np.cos(inc); G=-np.sin(w)*np.sin(Om)+np.cos(w)*np.cos(Om)*np.cos(inc)
    C=np.sin(w)*np.sin(inc); H=np.cos(w)*np.sin(inc)
    dDec=A*X+F*Y; dRA=B*X+G*Y; z=C*X+H*Y
    return dRA,dDec,z
# observed constraints
dmu_obs=np.array([-0.495,-0.338]); dmu_err=np.array([0.538,0.794])     # TGAS - DR3 (mas/yr), formal (field-calibrated z-scale 0.85/0.75 -> formal errors conservative)
tm=(1973.1,1959.7); t_dr1=2015.0
t_dr3=np.linspace(2014.57,2017.40,60)
t_rv=np.sort(rng.uniform(2014.57,2017.17,15)); sig_rv=2.65; chi2_obs_rv=stats.chi2.isf(0.275,14)
sig_al=0.12  # per-epoch AL precision proxy (mas) for G~11.8 incl. excess noise floor; RUWE limit 1.4 -> rms excess
Ms=[0.3,0.5,0.7,1.0,1.3]; As=np.logspace(0,3,31)
out={}
n=4000
for Mc in Ms:
    Mt=MG+Mc; B=Mc/Mt; fac=-(B-beta)          # photocentre = fac * relative orbit
    row=[]
    for a in As:
        P=np.sqrt(a**3/Mt)
        e=rng.uniform(0,0.9,n)[:,None]; w=rng.uniform(0,2*np.pi,n)[:,None]; Om=rng.uniform(0,2*np.pi,n)[:,None]
        inc=np.arccos(rng.uniform(-1,1,n))[:,None]; tp=rng.uniform(0,P,n)[:,None]
        def phot(t):
            dRA,dDec,z=orbit_xyz(np.atleast_1d(t)[None,:],P,e,w,Om,inc,tp,a)
            return fac*dRA*mas_per_au, fac*dDec*mas_per_au, z
        # DR3 pm: linear fit slope over window
        ra3,de3,_=phot(t_dr3); tt=t_dr3-t_dr3.mean()
        mu3_ra=(ra3*tt).sum(1)/(tt**2).sum(); mu3_de=(de3*tt).sum(1)/(tt**2).sum()
        res_ra=ra3-ra3.mean(1,keepdims=True)-mu3_ra[:,None]*tt; res_de=de3-de3.mean(1,keepdims=True)-mu3_de[:,None]*tt
        rms_wob=np.sqrt(np.mean(res_ra**2+res_de**2,axis=1)/2)
        # long-baseline (TGAS-like) pm
        ra_m,_,_=phot(tm[0]); _,de_m,_=phot(tm[1]); ra15,de15,_=phot(t_dr1)
        mul_ra=(ra15[:,0]-ra_m[:,0])/(t_dr1-tm[0]); mul_de=(de15[:,0]-de_m[:,0])/(t_dr1-tm[1])
        dmu=np.stack([mul_ra-mu3_ra, mul_de-mu3_de],1)
        chi2_pm=(((dmu-dmu_obs)/dmu_err)**2).sum(1); chi2_pm0=((dmu_obs/dmu_err)**2).sum()
        excl_pm=(chi2_pm-chi2_pm0)>stats.chi2.isf(0.05,2)
        # RV of G star: vz = -B * d(z_rel)/dt (au/yr -> km/s 4.74)
        _,_,z=orbit_xyz(t_rv[None,:],P,e,w,Om,inc,tp,a); _,_,z2=orbit_xyz(t_rv[None,:]+1e-3,P,e,w,Om,inc,tp,a)
        vrv=-B*(z2-z)/1e-3*4.74047
        lam=(((vrv-vrv.mean(1,keepdims=True))/sig_rv)**2).sum(1)
        p_rv=stats.ncx2.cdf(chi2_obs_rv,14,np.maximum(lam,1e-9))
        excl_rv=p_rv<0.05
        # RUWE: excess rms vs per-epoch precision; RUWE~sqrt(1+(rms/sig)^2); observed 1.026 -> exclude if predicted RUWE>1.4
        ruwe_pred=np.sqrt(1+(rms_wob/sig_al)**2)
        excl_ruwe=ruwe_pred>1.4
        excl=excl_pm|excl_rv|excl_ruwe
        row.append(dict(a_au=round(float(a),2),rho_arcsec=round(float(a/d_pc),4),P_yr=round(float(P),2),f_excl_pm=round(float(excl_pm.mean()),3),
                        f_excl_rv=round(float(excl_rv.mean()),3),f_excl_ruwe=round(float(excl_ruwe.mean()),3),f_excl_any=round(float(excl.mean()),3),
                        median_abs_dmu=round(float(np.median(np.hypot(dmu[:,0],dmu[:,1]))),3)))
    out[str(Mc)]=row
    print('Mc',Mc,[(r['a_au'],r['f_excl_any']) for r in row[::3]],flush=True)
json.dump(out,open(D+'results/pma_orbit_limits.json','w'),indent=1)
# summary: a below which >=50% / >=90% orbits excluded
summ={}
for Mc,row in out.items():
    a=np.array([r['a_au'] for r in row]); f=np.array([r['f_excl_any'] for r in row])
    def lim(th):
        ok=np.where(f>=th)[0]
        return float(a[ok].max()) if len(ok) else None
    summ[Mc]=dict(a_max_excl50=lim(0.5),a_max_excl90=lim(0.9),f_excl_at_30au=float(np.interp(np.log(30),np.log(a),f)),f_excl_at_50au=float(np.interp(np.log(50),np.log(a),f)),f_excl_at_100au=float(np.interp(np.log(100),np.log(a),f)))
print(json.dumps(summ,indent=1)); json.dump(summ,open(D+'results/pma_orbit_summary.json','w'),indent=1)
