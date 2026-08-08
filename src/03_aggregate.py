import pickle, numpy as np, pandas as pd
from collections import Counter
recs=pickle.load(open('output/boot_raw.pkl','rb')); B=len(recs)
TAIL={'clayton','gumbel','joe','bb1','bb7','student'}
names=[k for k in recs[0] if k!='nact']
rows=[]
for nm in names:
    inc=np.array([r[nm][0] for r in recs]); tree=np.array([r[nm][1] for r in recs])
    fam=[r[nm][2] for r in recs]; tau=np.array([r[nm][3] for r in recs]); lL=np.array([r[nm][4] for r in recs])
    Ds=[r[nm][5] for r in recs]
    sel=inc==1; t1=sel&(tree==1)
    tailsel=np.array([f in TAIL for f in fam])&sel
    # lambda_L>0 requiere familia con cola inferior
    lLpos=np.array([f in {'clayton','bb1','bb7'} or (f in {'gumbel','joe'} ) for f in fam])&sel
    md=Counter([D for D,s in zip(Ds,sel) if s]).most_common(1)
    rows.append(dict(pair=nm,EIP=inc.mean(),MCSE=inc.std(ddof=1)/np.sqrt(B),
        pT1=(tree==1).mean(),
        tau_med=np.median(tau[sel]) if sel.any() else np.nan,
        tau_lo=np.percentile(tau[sel],5) if sel.any() else np.nan,
        tau_hi=np.percentile(tau[sel],95) if sel.any() else np.nan,
        p_tailfam=(np.array([f in TAIL for f in fam])&t1).sum()/max(t1.sum(),1),
        p_lamLpos=((lL>0)&t1).sum()/max(t1.sum(),1),
        kap_med=np.median(lL[t1]) if t1.any() else np.nan,
        kap_lo=np.percentile(lL[t1],5) if t1.any() else np.nan,
        kap_hi=np.percentile(lL[t1],95) if t1.any() else np.nan,
        modalD=str(md[0][0]) if md else '', modalD_p=md[0][1]/max(sel.sum(),1) if md else np.nan,
        fams=', '.join(f"{f}:{c/max(sel.sum(),1):.2f}" for f,c in Counter([f for f,s in zip(fam,sel) if s]).most_common(3))))
R=pd.DataFrame(rows).sort_values('EIP',ascending=False)
R.to_csv('output/bootstrap_agg.csv',index=False)
pd.set_option('display.width',300); pd.set_option('display.max_colwidth',44)
print(f"B={B}   pares activos por remuestra: media {np.mean([r['nact'] for r in recs]):.1f}  "
      f"(5-95%: {np.percentile([r['nact'] for r in recs],5):.0f}-{np.percentile([r['nact'] for r in recs],95):.0f})\n")
print(R[['pair','EIP','MCSE','pT1','tau_med','tau_lo','tau_hi','p_lamLpos','kap_med','kap_lo','kap_hi']].round(3).to_string(index=False))
print("\n=== familias mas frecuentes (entre remuestras que retienen el par) ===")
for _,r in R.iterrows(): print(f"  {r['pair']:22s} {r['fams']}")
