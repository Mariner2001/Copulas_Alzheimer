"""Evidencia formal (DeltaBIC/Vuong) DENTRO de los analisis estratificados."""
import numpy as np, pandas as pd, vinelib as V, pyvinecopulib as pv
from scipy.stats import norm, rankdata, multivariate_normal as mvn
d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d); g=d.DX4.values
idx={l:i for i,l in enumerate(V.LAB)}
ALLFAM=[f for f in V.FAMS if f!=pv.BicopFamily.indep]
SYM=[pv.BicopFamily.gaussian,pv.BicopFamily.frank]
PAIRS=[('Memory','FDG'),('Memory','Amyloid'),('PCC','Precuneus')]
R=np.random.default_rng(20240607)

def evid(u,tag,pair):
    n=len(u)
    best=pv.Bicop.from_data(u,controls=pv.FitControlsBicop(family_set=ALLFAM,selection_criterion="bic"))
    sym =pv.Bicop.from_data(u,controls=pv.FitControlsBicop(family_set=SYM,selection_criterion="bic"))
    lo,hi=V._tail(best)
    l1=np.log(np.maximum(best.pdf(u),1e-300)); l0=np.log(np.maximum(sym.pdf(u),1e-300))
    k1=len(np.atleast_1d(best.parameters).ravel()); k0=len(np.atleast_1d(sym.parameters).ravel())
    dif=l1-l0; sd=dif.std(ddof=1)
    vu=(dif.sum()-(k1-k0)*np.log(n)/2)/(np.sqrt(n)*sd) if sd>1e-12 else np.nan
    dB=(-2*l0.sum()+k0*np.log(n))-(-2*l1.sum()+k1*np.log(n))
    q=0.2; LL=np.mean((u[:,0]<q)&(u[:,1]<q)); UR=np.mean((u[:,0]>1-q)&(u[:,1]>1-q))
    rho=np.sin(np.pi*best.tau/2)
    gb=mvn(mean=[0,0],cov=[[1,rho],[rho,1]]).cdf([norm.ppf(q),norm.ppf(q)])
    return dict(analisis=tag,pair=f"{pair[0]}-{pair[1]}",n=n,fam=str(best.family).split('.')[-1],
        tau=best.tau,lamL=lo,dBIC=dB,vuong=vu,p=2*(1-norm.cdf(abs(vu))) if np.isfinite(vu) else np.nan,
        LLobs=100*LL,URobs=100*UR,ratio=LL/UR if UR>0 else np.nan,LLgauss=100*gb,
        z_gauss=(LL-gb)/np.sqrt(gb*(1-gb)/n))

rows=[]
U0=V.pseudo(X,R); Ua=V.pseudo(X,R,groups=g)
for p in PAIRS:
    u=U0[:,[idx[p[0]],idx[p[1]]]]; rows.append(evid(u,'agrupado',p))
    u=Ua[:,[idx[p[0]],idx[p[1]]]]; rows.append(evid(u,'rangos intra-estrato',p))
for lab,m in [('CN+EMCI',np.isin(g,['CN','EMCI'])),('LMCI+AD',np.isin(g,['LMCI','AD']))]:
    Us=V.pseudo(X[m].reset_index(drop=True),R)
    for p in PAIRS: rows.append(evid(Us[:,[idx[p[0]],idx[p[1]]]],lab,p))
E=pd.DataFrame(rows); E.to_csv('output/strat_evidence.csv',index=False)
pd.set_option('display.width',260)
print(E.round(3).to_string(index=False))
