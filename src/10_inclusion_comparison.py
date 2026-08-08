"""Compare included (n=630) with excluded eligible participants, and print the nested
attrition flow.  Addresses complete-case selection bias."""
import pandas as pd, numpy as np
from scipy import stats
U='data/'
adsl=pd.read_csv('data/adsl.csv').rename(columns={'SUBJID':'RID'})
coh=adsl[adsl.ORIGPROT.isin(['ADNIGO','ADNI2'])].copy()
arm=pd.read_csv('data/arm.csv'); arm=arm[arm.ORIGPROT.isin(['ADNIGO','ADNI2'])]
def dx4(s):
    s=str(s)
    for k,v in [('EMCI','EMCI'),('LMCI','LMCI'),('MCI','LMCI'),('AD','AD'),('NL','CN'),('SMC','SMC')]:
        if s.startswith(k): return v
    return np.nan
arm['DX4']=arm.ARM.map(dx4); arm=arm.dropna(subset=['DX4']).drop_duplicates('RID')[['RID','DX4']]
coh=coh.merge(arm,on='RID',how='left')
coh=coh[coh.DX4.isin(['CN','EMCI','LMCI','AD'])]          # mismo filtro que el articulo
final=pd.read_csv('data/adni_merged.csv')
coh['incluido']=coh.RID.isin(final.RID)
print(f"elegibles (ADNI-GO/2, brazo CN/EMCI/LMCI/AD): {len(coh)}")
print(f"   incluidos {coh.incluido.sum()}   excluidos {(~coh.incluido).sum()}\n")
inc=coh[coh.incluido]; exc=coh[~coh.incluido]
print(f"{'variable':22s}{'incluidos (n=630)':>22s}{'excluidos (n=642)':>22s}{'p':>10s}")
for lab,col in [('Edad, media (DE)','AGE'),('Educacion, media (DE)','EDUC')]:
    a,b=inc[col].dropna(),exc[col].dropna()
    print(f"{lab:22s}{f'{a.mean():.1f} ({a.std():.1f})':>22s}{f'{b.mean():.1f} ({b.std():.1f})':>22s}{stats.ttest_ind(a,b,equal_var=False).pvalue:>10.3f}")
a=(inc.SEX=='Female'); b=(exc.SEX=='Female')
ct=np.array([[a.sum(),len(a)-a.sum()],[b.sum(),len(b)-b.sum()]])
print(f"{'Mujeres, n (%)':22s}{f'{a.sum()} ({100*a.mean():.0f}%)':>22s}{f'{b.sum()} ({100*b.mean():.0f}%)':>22s}{stats.chi2_contingency(ct)[1]:>10.3f}")
t=pd.crosstab(coh.DX4,coh.incluido).reindex(['CN','EMCI','LMCI','AD'])
print(f"\n{'estadio':10s}{'incluidos':>12s}{'excluidos':>12s}{'% incluido':>12s}")
for g,r in t.iterrows(): print(f"{g:10s}{r[True]:>12d}{r[False]:>12d}{100*r[True]/(r[True]+r[False]):>11.0f}%")
print(f"\nchi2 estadio x inclusion: p = {stats.chi2_contingency(t.values)[1]:.4f}")
