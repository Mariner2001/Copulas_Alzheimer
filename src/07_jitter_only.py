"""Estabilidad SOLO-JITTER: muestra fija, solo cambia la extension continua.
Separa la ambiguedad de identificacion (Prop. 4) de la variabilidad muestral."""
import numpy as np, pandas as pd, vinelib as V, pyvinecopulib as pv, pickle, os, sys, time
from collections import Counter
d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d)
KEYS={'Sex-Thalamus':frozenset({'Sex','Thalamus'}),'Sex-Hippocampus':frozenset({'Sex','Hippocampus'}),
      'Sex-Memory':frozenset({'Sex','Memory'}),'Sex-ExecFun':frozenset({'Sex','ExecFun'}),
      'Sex-PCC':frozenset({'Sex','PCC'}),'Sex-Precuneus':frozenset({'Sex','Precuneus'}),
      'APOE4-Amyloid':frozenset({'APOE4','Amyloid'}),'Education-ExecFun':frozenset({'Education','ExecFun'}),
      'Memory-FDG':frozenset({'Memory','FDG'}),'PCC-Precuneus':frozenset({'PCC','Precuneus'})}
ST='output/jit_raw.pkl'; recs=pickle.load(open(ST,'rb')) if os.path.exists(ST) else []
rng=np.random.default_rng(555+len(recs))
ctl=pv.FitControlsVinecop(family_set=V.FAMS,selection_criterion="mbicv",trunc_lvl=3,num_threads=1)
for b in range(int(sys.argv[1])):
    U=V.pseudo(X,rng)                       # MISMA muestra, jitter distinto
    E=V.named_edges(pv.Vinecop.from_data(U,controls=ctl))
    seen={r['pair']:r for _,r in E.iterrows()}
    recs.append({k:(1,int(seen[v]['tree']),seen[v]['fam'],float(seen[v]['tau'])) if v in seen else (0,0,'indep',0.0)
                 for k,v in KEYS.items()})
    if (b+1)%10==0: pickle.dump(recs,open(ST,'wb'))
pickle.dump(recs,open(ST,'wb')); print(f"{len(recs)} extracciones de jitter")
