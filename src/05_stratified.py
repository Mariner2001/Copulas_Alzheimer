"""Heterogeneidad diagnostica: (a) rangos intra-estrato, (b) DX como nodo, (c) reajustes estratificados."""
import numpy as np, pandas as pd, vinelib as V, pyvinecopulib as pv
d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d)
g=d.DX4.values
PAIRS=[('Memory','FDG'),('Memory','Amyloid'),('PCC','Precuneus')]
idx={l:i for i,l in enumerate(V.LAB)}
ALLFAM=[f for f in V.FAMS if f!=pv.BicopFamily.indep]
out=[]

def pairfit(U,a,b,tag,nn):
    u=U[:,[idx[a],idx[b]]]
    bc=pv.Bicop.from_data(u,controls=pv.FitControlsBicop(family_set=ALLFAM,selection_criterion="bic"))
    lo,hi=V._tail(bc)
    out.append(dict(analisis=tag,pair=f"{a}-{b}",n=nn,fam=str(bc.family).split('.')[-1],
                    rot=bc.rotation,tau=bc.tau,lamL=lo))

R=np.random.default_rng(20240607)
# (0) referencia: agrupado
U0=V.pseudo(X,R)
for a,b in PAIRS: pairfit(U0,a,b,"(0) agrupado (referencia)",len(d))
# (a) rangos calculados DENTRO de cada estrato
Ua=V.pseudo(X,R,groups=g)
for a,b in PAIRS: pairfit(Ua,a,b,"(a) rangos intra-estrato",len(d))
# (c) reajustes estratificados
for lab,mask in [("(c1) CN+EMCI",np.isin(g,['CN','EMCI'])),("(c2) LMCI+AD",np.isin(g,['LMCI','AD']))]:
    Xs=X[mask].reset_index(drop=True); Us=V.pseudo(Xs,R)
    for a,b in PAIRS: pairfit(Us,a,b,lab,int(mask.sum()))
for lab in ['CN','EMCI','LMCI','AD']:
    m=g==lab; Xs=X[m].reset_index(drop=True); Us=V.pseudo(Xs,R)
    for a,b in PAIRS: pairfit(Us,a,b,f"(c3) solo {lab}",int(m.sum()))

# (b) diagnostico como nodo ordinal: lambda_L del par condicionado a DX
order={'CN':0,'EMCI':1,'LMCI':2,'AD':3}
dxv=np.array([order[x] for x in g],dtype=float)
Xb=X.copy(); Xb['DX']= -dxv                      # orientado: menor = mas patologia
VARS2=V.VARS+['DX']; LAB2=V.LAB+['DX']
n=len(Xb); Ub=np.empty((n,len(VARS2)))
from scipy.stats import rankdata
for j,c in enumerate(VARS2):
    v=Xb[c].values.astype(float).copy()
    if c in V.DISCRETE or c=='DX': v=v+R.uniform(0,1,n)
    Ub[:,j]=rankdata(v)/(n+1)
cop=pv.Vinecop.from_data(Ub,controls=V.ctrl())
i2={l:i for i,l in enumerate(LAB2)}
S=cop.structure; M=np.array(S.matrix); dd=S.dim
print("\n(b) diagnostico como nodo (p=15): aristas de los tres pares")
for t in range(S.trunc_lvl):
    for e in range(dd-1-t):
        bc=cop.pair_copulas[t][e]
        if str(bc.family)=='BicopFamily.indep': continue
        a=LAB2[int(M[dd-1-e,e])-1]; b=LAB2[int(M[t,e])-1]
        D=tuple(LAB2[int(M[k,e])-1] for k in range(t))
        if frozenset({a,b}) in [frozenset(p) for p in PAIRS]:
            lo,_=V._tail(bc)
            print(f"    {a}-{b} | D={D or '()'}  arbol T{t+1}  {str(bc.family).split('.')[-1]}  tau={bc.tau:.3f}  lamL={lo:.3f}")

D=pd.DataFrame(out); D.to_csv('output/stratified.csv',index=False)
pd.set_option('display.width',200)
print("\n=== (a) y (c): familia y lambda_L por analisis ===")
print(D.pivot_table(index='analisis',columns='pair',values='lamL').round(3).to_string())
print("\n=== familias seleccionadas ===")
print(D.pivot_table(index='analisis',columns='pair',values='fam',aggfunc='first').to_string())
