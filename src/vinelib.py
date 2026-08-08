import numpy as np, pandas as pd, pyvinecopulib as pv
from scipy.stats import rankdata

VARS=['AGE','SEX_bin','EDUC','ADNI_MEM','ADNI_EF','APOE4_bin','AMY','FDG',
      'Hippocampus','Caudate','Putamen','Thalamus','PCC','Precuneus']
LAB =['Age','Sex','Education','Memory','ExecFun','APOE4','Amyloid','FDG',
      'Hippocampus','Caudate','Putamen','Thalamus','PCC','Precuneus']
DISCRETE={'SEX_bin','EDUC','APOE4_bin'}          # margenes discretos
REFLECT ={'AGE','AMY','APOE4_bin'}                # orientar: menor = mas patologia

FAMS=[pv.BicopFamily.indep,pv.BicopFamily.gaussian,pv.BicopFamily.student,
      pv.BicopFamily.clayton,pv.BicopFamily.gumbel,pv.BicopFamily.frank,
      pv.BicopFamily.joe,pv.BicopFamily.bb1,pv.BicopFamily.bb7]
SYM  =[pv.BicopFamily.gaussian,pv.BicopFamily.frank]

def orient(df):
    X=df[VARS].astype(float).copy()
    for c in REFLECT: X[c]=-X[c]
    return X

def pseudo(X, rng, groups=None):
    """Pseudo-observaciones. Jitter en discretas. groups -> rangos intra-estrato."""
    n=len(X); U=np.empty((n,len(VARS)))
    for j,c in enumerate(VARS):
        v=X[c].values.astype(float).copy()
        if c in DISCRETE: v = v + rng.uniform(0,1,n)      # continuous extension
        if groups is None:
            U[:,j]=rankdata(v)/(n+1)
        else:
            u=np.empty(n)
            for g in np.unique(groups):
                m=groups==g; u[m]=rankdata(v[m])/(m.sum()+1)
            U[:,j]=u
    return U

def ctrl(trunc=None,fams=FAMS,crit="mbicv"):
    kw=dict(family_set=fams, selection_criterion=crit, num_threads=4)
    if trunc is None: kw['select_trunc_lvl']=True
    else: kw['trunc_lvl']=trunc
    return pv.FitControlsVinecop(**kw)

def fit(U,trunc=None,crit="mbicv"):
    return pv.Vinecop.from_data(U, controls=ctrl(trunc,crit=crit))

def edges(cop):
    """Devuelve lista de aristas activas: (i,j,D,arbol,familia,tau,lamL,lamU)."""
    out=[]
    order=np.array(cop.matrix)          # matriz R-vine
    for t,tree in enumerate(cop.pair_copulas):
        for e,bc in enumerate(tree):
            if bc.family==pv.BicopFamily.indep: continue
            s=cop.get_struct_array if False else None
            out.append(dict(tree=t+1,edge=e,fam=str(bc.family).split('.')[-1],
                            rot=bc.rotation,tau=bc.tau,
                            lamL=bc.tau and _tail(bc)[0], lamU=_tail(bc)[1],
                            bc=bc))
    return out

def _tail(bc):
    """Coeficientes de dependencia de cola ANALITICOS segun familia y rotacion."""
    import numpy as np
    from scipy.stats import t as tdist
    f=str(bc.family).split('.')[-1]; p=np.atleast_1d(bc.parameters).ravel().astype(float); r=int(bc.rotation)
    lo=hi=0.0
    if f=='clayton':
        th=p[0]; lo=2**(-1/th) if th>0 else 0.0; hi=0.0
    elif f=='gumbel':
        th=p[0]; hi=2-2**(1/th) if th>=1 else 0.0; lo=0.0
    elif f=='joe':
        th=p[0]; hi=2-2**(1/th) if th>=1 else 0.0; lo=0.0
    elif f=='bb1':
        th,de=p[0],p[1]; lo=2**(-1/(th*de)); hi=2-2**(1/de)
    elif f=='bb7':
        th,de=p[0],p[1]; lo=2**(-1/de); hi=2-2**(1/th)
    elif f=='student':
        rho,nu=p[0],p[1]
        lo=hi=2*tdist.cdf(-np.sqrt((nu+1)*(1-rho)/(1+rho)), nu+1)
    # gaussian, frank, indep -> 0,0
    if r==180: lo,hi=hi,lo          # copula de supervivencia: intercambia colas
    elif r in (90,270): lo,hi=0.0,0.0   # rotaciones de dependencia negativa
    return float(lo), float(hi)

def edge_sets(cop):
    """(tree, edge) -> (a, b, D) con indices 0-based sobre VARS."""
    S=cop.structure; d=S.dim; res={}
    M=np.array(S.matrix)          # matriz d x d, R-vine en formato triangular
    for t in range(S.trunc_lvl):
        for e in range(d-1-t):
            a=int(M[d-1-e, e])-1                     # variable "diagonal"
            b=int(M[t, e])-1                         # variable conditioned
            D=tuple(sorted(int(M[k,e])-1 for k in range(t)))
            res[(t,e)]=(a,b,D)
    return res

def named_edges(cop, min_tau=0.0):
    out=[]
    ES=edge_sets(cop)
    for t,tree in enumerate(cop.pair_copulas):
        for e,bc in enumerate(tree):
            if str(bc.family)=='BicopFamily.indep': continue
            a,b,D=ES[(t,e)]
            lo,hi=_tail(bc)
            out.append(dict(tree=t+1, a=LAB[a], b=LAB[b],
                            D=tuple(LAB[k] for k in D),
                            ia=a, ib=b, iD=D,
                            fam=str(bc.family).split('.')[-1], rot=bc.rotation,
                            tau=bc.tau, lamL=lo, lamU=hi,
                            pair=frozenset((LAB[a],LAB[b]))))
    return pd.DataFrame(out)
