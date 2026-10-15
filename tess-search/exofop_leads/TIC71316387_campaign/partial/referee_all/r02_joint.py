# Referee: re-do the localisation joint fit from per-sector offsets (own least squares)
import numpy as np, json
from scipy.stats import chi2 as C2
R='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v3/results/'
d=json.load(open(R+'s05_baseline.json'))
obs=np.array(d['obs']); cov=np.array(d['cov']); u=np.array(d['u_to_boresight']); K=np.array(d['KR_rel'])
lab=d['labels']; typ=np.array([0 if 'deep' in l else 1 for l in lab])
W=np.array([np.linalg.inv(c) for c in cov])
def solve(cols):
    # linear model obs_i = sum_j p_j * cols_j(i) (2-vectors); GLS
    A=np.zeros((len(cols),len(cols))); b=np.zeros(len(cols))
    for i in range(len(obs)):
        X=np.array([c[i] for c in cols]).T  # 2 x np
        A+=X.T@W[i]@X; b+=X.T@W[i]@obs[i]
    p=np.linalg.solve(A,b); cv=np.linalg.inv(A)
    r=[obs[i]-np.array([c[i] for c in cols]).T@p for i in range(len(obs))]
    ch=sum(r[i]@W[i]@r[i] for i in range(len(obs)))
    return p,np.sqrt(np.diag(cv)),ch
E=[np.array([1,0])]*12; N=[np.array([0,1])]*12
away=[-K[i]*u[i] for i in range(12)]
cd=[away[i] if typ[i]==0 else np.zeros(2) for i in range(12)]
cs=[away[i] if typ[i]==1 else np.zeros(2) for i in range(12)]
res={}
for name,cols in [('M0_fixed_offset',[E,N]),('M1_colour_only',[cd,cs]),('M2_offset_plus_colour',[E,N,cd,cs])]:
    p,e,ch=solve(cols); dof=24-len(p)
    res[name]=dict(p=p.tolist(),err_formal=e.tolist(),chi2=ch,dof=dof,p_chi2=float(C2.sf(ch,dof)))
    print(name,np.round(p,3),np.round(e,3),round(ch,2),dof,C2.sf(ch,dof))
# radial/tangential decomposition of raw offsets
for i in range(12):
    ra=obs[i]@(-u[i]); ta=obs[i]@np.array([-u[i][1],u[i][0]])
    sr=np.sqrt((-u[i])@cov[i]@(-u[i])); st=np.sqrt(np.array([-u[i][1],u[i][0]])@cov[i]@np.array([-u[i][1],u[i][0]]))
    print(lab[i],'away=%.2f+-%.2f tan=%.2f+-%.2f (%.1f sig)'%(ra,sr,ta,st,ta/st))
# drop S105 in M2
keep=[i for i in range(12) if 'S105' not in lab[i]]
obs2=obs[keep];cov2=cov[keep]
o0,c0,W0=obs,cov,W
obs=obs2;cov=cov2;W=np.array([np.linalg.inv(c) for c in cov2])
cols=[[E[i] for i in keep],[N[i] for i in keep],[cd[i] for i in keep],[cs[i] for i in keep]]
p,e,ch=solve(cols); print('M2 drop S105',np.round(p,3),np.round(e,3),round(ch,2))
res['M2_dropS105']=dict(p=p.tolist(),err=e.tolist(),chi2=ch)
json.dump(res,open('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all/r02_joint.json','w'),indent=1)
