import numpy as np, pandas as pd, vinelib as V, pyvinecopulib as pv, sys, os, time, pickle
d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d); n=len(d)
CHUNK=int(sys.argv[1])
PAIRS={'Memory-FDG':('Memory','FDG'),'Memory-Amyloid':('Memory','Amyloid'),
       'PCC-Precuneus':('PCC','Precuneus'),'Memory-Hippocampus':('Memory','Hippocampus'),
       'Age-Hippocampus':('Age','Hippocampus'),'Sex-Memory':('Sex','Memory'),
       'Sex-Thalamus':('Sex','Thalamus'),'Sex-Hippocampus':('Sex','Hippocampus'),
       'Sex-ExecFun':('Sex','ExecFun'),'Sex-PCC':('Sex','PCC'),'Sex-Precuneus':('Sex','Precuneus'),
       'APOE4-Amyloid':('APOE4','Amyloid'),'Education-ExecFun':('Education','ExecFun'),
       'Age-Memory':('Age','Memory'),'ExecFun-Memory':('ExecFun','Memory'),
       'Caudate-Putamen':('Caudate','Putamen'),'Thalamus-Hippocampus':('Thalamus','Hippocampus')}
KEYS={k:frozenset(v) for k,v in PAIRS.items()}
STORE='output/boot_raw.pkl'
recs = pickle.load(open(STORE,'rb')) if os.path.exists(STORE) else []
done=len(recs)
rng=np.random.default_rng(90210+done)
ctl=pv.FitControlsVinecop(family_set=V.FAMS,selection_criterion="mbicv",trunc_lvl=3,num_threads=1)
t0=time.time()
for b in range(CHUNK):
    idx=rng.integers(0,n,n)
    U=V.pseudo(X.iloc[idx].reset_index(drop=True),rng)
    cop=pv.Vinecop.from_data(U,controls=ctl)
    E=V.named_edges(cop)
    seen={r['pair']:r for _,r in E.iterrows()}
    rec={'nact':len(E)}
    for name,pk in KEYS.items():
        if pk in seen:
            r=seen[pk]; rec[name]=(1,int(r['tree']),r['fam'],float(r['tau']),float(r['lamL']),tuple(sorted(r['D'])))
        else: rec[name]=(0,0,'indep',0.0,0.0,())
    recs.append(rec)
    if (b+1)%10==0: pickle.dump(recs,open(STORE,'wb'))   # guardado incremental
pickle.dump(recs,open(STORE,'wb'))
print(f"acumuladas {len(recs)} remuestras (+{CHUNK} en {time.time()-t0:.0f}s)")
