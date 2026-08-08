"""Evidencia de dependencia de cola: DeltaBIC, Vuong, cuotas de esquina, lambda_L(u)."""
import numpy as np, pandas as pd, vinelib as V, pyvinecopulib as pv
from scipy.stats import norm, rankdata, multivariate_normal as mvn

d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d); n=len(d)
rng=np.random.default_rng(20240607)
U=V.pseudo(X,rng)
idx={l:i for i,l in enumerate(V.LAB)}
PAIRS=[('Memory','FDG'),('Memory','Amyloid'),('PCC','Precuneus'),
       ('Memory','Hippocampus'),('ExecFun','Memory')]
SYMFAM=[pv.BicopFamily.gaussian,pv.BicopFamily.frank]
ALLFAM=[f for f in V.FAMS if f!=pv.BicopFamily.indep]

def gauss_box(q,tau):
    rho=np.sin(np.pi*tau/2)
    return mvn(mean=[0,0],cov=[[1,rho],[rho,1]]).cdf([norm.ppf(q),norm.ppf(q)])

rows=[]; curves={}
for a,b in PAIRS:
    u=U[:,[idx[a],idx[b]]]
    best=pv.Bicop.from_data(u,controls=pv.FitControlsBicop(family_set=ALLFAM,selection_criterion="bic"))
    bsym=pv.Bicop.from_data(u,controls=pv.FitControlsBicop(family_set=SYMFAM,selection_criterion="bic"))
    lo,hi=V._tail(best)
    # log-verosimilitudes puntuales para Vuong
    l1=np.log(np.maximum(best.pdf(u),1e-300)); l0=np.log(np.maximum(bsym.pdf(u),1e-300))
    dif=l1-l0
    k1=len(np.atleast_1d(best.parameters).ravel()); k0=len(np.atleast_1d(bsym.parameters).ravel())
    corr=(k1-k0)*np.log(n)/2                       # correccion Schwarz
    vuong=(dif.sum()-corr)/(np.sqrt(n)*dif.std(ddof=1))
    dBIC=(-2*l0.sum()+k0*np.log(n))-(-2*l1.sum()+k1*np.log(n))   # >0 favorece la de cola
    # cuotas de esquina observadas
    q=0.2
    LL=np.mean((u[:,0]<q)&(u[:,1]<q)); UR=np.mean((u[:,0]>1-q)&(u[:,1]>1-q))
    tau=best.tau
    rows.append(dict(pair=f"{a}-{b}", fam=str(best.family).split('.')[-1], rot=best.rotation,
        tau=tau, lamL=lo, lamU=hi, sym_fam=str(bsym.family).split('.')[-1],
        dBIC=dBIC, vuong=vuong, p_vuong=2*(1-norm.cdf(abs(vuong))),
        LL_obs=100*LL, UR_obs=100*UR, ratio=LL/UR if UR>0 else np.nan,
        LL_gauss=100*gauss_box(q,tau), LL_indep=4.0,
        z_vs_gauss=(LL-gauss_box(q,tau))/np.sqrt(gauss_box(q,tau)*(1-gauss_box(q,tau))/n)))
    # curva lambda_L(u) empirica vs gaussiana con la misma tau
    us=np.linspace(0.02,0.20,25)
    emp=[np.mean((u[:,0]<t)&(u[:,1]<t))/t for t in us]
    rho=np.sin(np.pi*tau/2)
    gau=[mvn(mean=[0,0],cov=[[1,rho],[rho,1]]).cdf([norm.ppf(t),norm.ppf(t)])/t for t in us]
    curves[f"{a}-{b}"]=(us,np.array(emp),np.array(gau))

R=pd.DataFrame(rows); R.to_csv('output/tail_evidence.csv',index=False)
pd.set_option('display.width',260)
print(R.round(3).to_string(index=False))
np.savez('output/lambda_curves.npz',**{k:np.vstack(v) for k,v in curves.items()})
