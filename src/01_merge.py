import pandas as pd, numpy as np, re, glob, os
def pick(pat):
    f=sorted(glob.glob('data/'+pat))
    assert f, f'No encuentro data/{pat}'
    return f[-1]
U='data/'
flow=[]
def log(s,n): flow.append((s,n)); print(f"{s:55s} n={n}")

# ---------- 1. Cohorte: ADNI-GO / ADNI-2 ----------
adsl=pd.read_csv('data/adsl.csv'); adsl=adsl.rename(columns={'SUBJID':'RID'})
log("ADSL total", len(adsl))
coh=adsl[adsl.ORIGPROT.isin(['ADNIGO','ADNI2'])].copy()
log("ORIGPROT in {ADNIGO, ADNI2}", len(coh))

# ---------- 2. Diagnostico 4 niveles desde ARM ----------
arm=pd.read_csv('data/arm.csv')
arm=arm[arm.ORIGPROT.isin(['ADNIGO','ADNI2'])]
def dx4(s):
    s=str(s)
    if s.startswith('EMCI'): return 'EMCI'
    if s.startswith('LMCI'): return 'LMCI'
    if s.startswith('MCI'):  return 'LMCI'
    if s.startswith('AD'):   return 'AD'
    if s.startswith('NL'):   return 'CN'
    if s.startswith('SMC'):  return 'SMC'
    return np.nan
arm['DX4']=arm['ARM'].map(dx4)
arm=arm.dropna(subset=['DX4']).drop_duplicates('RID')[['RID','DX4']]
coh=coh.merge(arm,on='RID',how='left')
log("con brazo de reclutamiento (DX4)", coh.DX4.notna().sum())
print("   ", coh.DX4.value_counts(dropna=False).to_dict())

# ---------- 3. APOE4 ----------
def apoe4(g):
    if pd.isna(g): return np.nan
    return str(g).count('4')
coh['APOE4']=coh['APOE'].map(apoe4)
coh['APOE4_bin']=(coh['APOE4']>=1).astype(float); coh.loc[coh.APOE4.isna(),'APOE4_bin']=np.nan

# ---------- 4. Cognicion (baseline) ----------
uw=pd.read_csv('data/uwnp.csv')
uw=uw[uw.VISCODE2=='bl'].drop_duplicates('RID')[['RID','ADNI_MEM','ADNI_EF']]
coh=coh.merge(uw,on='RID',how='left')
log("con ADNI-MEM/EF baseline", coh.ADNI_MEM.notna().sum())

# ---------- 5. FDG MetaROI (pons-normalizado) ----------
fdg=pd.read_csv(pick('UCBERKELEYFDG_8mm*.csv'),low_memory=False)
fdg=fdg[fdg.VISCODE2=='bl']
piv=fdg.pivot_table(index='RID',columns='ROINAME',values='MEAN')
piv['FDG']=piv['MetaROI']  # ya normalizado a pons en esta version
coh=coh.merge(piv[['FDG']].reset_index(),on='RID',how='left')
log("con FDG baseline", coh.FDG.notna().sum())

# ---------- 6. Amiloide (FBP, whole-cerebellum) ----------
amy=pd.read_csv(pick('UCBERKELEY_AMY_6MM*.csv'),low_memory=False)
amy=amy[(amy.VISCODE2=='bl')&(amy.TRACER=='FBP')]
assert np.allclose(amy.WHOLECEREBELLUM_SUVR.dropna(),1.0), "SUMMARY_SUVR no esta normalizado a cerebelo"
amy=amy.drop_duplicates('RID')[['RID','SUMMARY_SUVR']].rename(columns={'SUMMARY_SUVR':'AMY'})
coh=coh.merge(amy,on='RID',how='left')
log("con amiloide FBP baseline", coh.AMY.notna().sum())

# ---------- 7. Volumenes FreeSurfer 5.1 ----------
fs=pd.read_csv(pick('UCSFFSX51*.csv'),low_memory=False)
fs=fs[fs.VISCODE2.isin(['scmri','bl'])]
log("FreeSurfer filas baseline (pre-QC)", len(fs))
fs=fs[fs.OVERALLQC=='Pass']
log("FreeSurfer OVERALLQC = Pass", len(fs))
CODES={'Hippocampus':('ST29SV','ST88SV'),'Caudate':('ST16SV','ST75SV'),
       'Putamen':('ST53SV','ST112SV'),'Thalamus':('ST61SV','ST120SV'),
       'PCC':('ST50CV','ST109CV'),'Precuneus':('ST52CV','ST111CV')}
keep=['RID','ST10CV']+[c for lr in CODES.values() for c in lr]
fs=fs[keep].copy()
for k,(l,r) in CODES.items(): fs[k]=(fs[l]+fs[r])/2.0
fs=fs.rename(columns={'ST10CV':'ICV'})
fs=fs.dropna(subset=['ICV']).sort_values('ICV').drop_duplicates('RID')
for k in CODES: fs[k]=fs[k]/fs['ICV']*1000.0   # escalado por ICV
coh=coh.merge(fs[['RID','ICV']+list(CODES)],on='RID',how='left')
log("con volumenes FreeSurfer QC-pass", coh.Hippocampus.notna().sum())

VARS14=['AGE','SEX','EDUC','ADNI_MEM','ADNI_EF','APOE4_bin','AMY','FDG',
        'Hippocampus','Caudate','Putamen','Thalamus','PCC','Precuneus']
coh=coh[coh.DX4.isin(['CN','EMCI','LMCI','AD'])]
log("excluidos SMC / sin brazo", len(coh))
final=coh.dropna(subset=VARS14).copy()
log(">>> CASOS COMPLETOS (p=14)", len(final))
print()
print(final.DX4.value_counts().reindex(['CN','EMCI','LMCI','AD']).to_string())
final['SEX_bin']=(final.SEX=='Female').astype(int)
final=final.reset_index(drop=True); final.index.name='id'
final.to_csv('data/adni_merged.csv')
pd.DataFrame(flow,columns=['paso','n']).to_csv('output/flow.csv',index=False)
print("\nGuardado adni_merged.csv", final.shape)

# exportar la lista de aristas activas del vine (para gcgm_benchmark.R)
import sys; sys.path.insert(0,'src')
try:
    import vinelib as V, numpy as _np
    _cop = V.fit(V.pseudo(V.orient(final), _np.random.default_rng(20240607)))
    V.named_edges(_cop)[['tree','a','b','fam','tau','lamL']].to_csv('output/vine_edges.csv',index=False)
    print("output/vine_edges.csv escrito")
except Exception as e:
    print("aviso: no se pudo escribir vine_edges.csv:", e)
