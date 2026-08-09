"""Mixture benchmark: how much of the pooled corner asymmetry can diagnostic-stage
heterogeneity account for?  Preserves each group's observed margins exactly and imposes a
Gaussian (tail-independent) copula within group, then pools at observed group sizes.
Writes output/mixture_benchmark.csv."""
import sys; sys.path.insert(0,'src')
import numpy as np, pandas as pd
from scipy.stats import rankdata, norm
d=pd.read_csv('data/adni_merged.csv'); d['AMY']=-d['AMY']      # orientado
G=['CN','EMCI','LMCI','AD']

def cr(x,y,q=.2):
    u,v=rankdata(x)/(len(x)+1), rankdata(y)/(len(y)+1)
    LL=np.mean((u<q)&(v<q)); UR=np.mean((u>1-q)&(v>1-q))
    return LL/UR if UR>0 else np.nan

def sim(c1,c2,tau_w,B=2000,seed=0):
    rho=np.sin(np.pi*tau_w/2); rng=np.random.default_rng(seed); out=[]
    grp=[(np.sort(d[d.DX4==g][c1].values), np.sort(d[d.DX4==g][c2].values)) for g in G]
    for b in range(B):
        xs,ys=[],[]
        for X,Y in grp:
            n=len(X)
            z=rng.multivariate_normal([0,0],[[1,rho],[rho,1]],n)
            u,v=norm.cdf(z[:,0]),norm.cdf(z[:,1])
            # cuantiles empiricos del grupo -> margenes intra-grupo EXACTOS
            xs.append(X[np.clip((u*n).astype(int),0,n-1)])
            ys.append(Y[np.clip((v*n).astype(int),0,n-1)])
        out.append(cr(np.concatenate(xs),np.concatenate(ys)))
    return np.array(out)


OBS={'Memory-FDG':(('ADNI_MEM','FDG'),0.142,2.00),
     'Memory-Amyloid':(('ADNI_MEM','AMY'),0.175,1.88),
     'PCC-Precuneus':(('PCC','Precuneus'),0.250,1.24)}
import csv, os
rows=[]
print(f"{'pair':18s}{'R obs':>8s}{'R under pure mixing (median [5-95%])':>40s}{'p':>8s}")
for lab,((c1,c2),tw,o) in OBS.items():
    R=sim(c1,c2,tw,B=4000,seed=7)
    lo,hi=np.percentile(R,[5,95]); p=(np.sum(R>=o)+1)/(len(R)+1)
    print(f"{lab:18s}{o:8.2f}{f'{np.median(R):.2f}  [{lo:.2f}, {hi:.2f}]':>40s}{p:8.3f}")
    rows.append(dict(pair=lab,tau_within=tw,R_observed=o,R_mixture_median=np.median(R),
                     R_mixture_lo=lo,R_mixture_hi=hi,p_value=p,
                     share_of_excess=(np.median(R)-1)/(o-1)))
os.makedirs('output',exist_ok=True)
with open('output/mixture_benchmark.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\noutput/mixture_benchmark.csv escrito")
