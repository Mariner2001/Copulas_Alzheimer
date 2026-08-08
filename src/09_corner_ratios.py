"""Punto 3: IC bootstrap para log R.  Punto 26: sensibilidad del ratio a u."""
import sys; sys.path.insert(0,'src')
import numpy as np, pandas as pd, vinelib as V
d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d); g=d.DX4.values
idx={l:i for i,l in enumerate(V.LAB)}
rng=np.random.default_rng(20240607)
U0=V.pseudo(X,rng); Ua=V.pseudo(X,rng,groups=g)
PAIRS=[('Memory','FDG'),('Memory','Amyloid'),('PCC','Precuneus'),('ExecFun','Memory')]

def ratio(u,q):
    LL=np.mean((u[:,0]<q)&(u[:,1]<q)); UR=np.mean((u[:,0]>1-q)&(u[:,1]>1-q))
    return LL,UR,(LL/UR if UR>0 else np.nan)

def boot_logR(u,q,B=5000,seed=1):
    r=np.random.default_rng(seed); n=len(u); out=[]
    for _ in range(B):
        i=r.integers(0,n,n); LL,UR,R=ratio(u[i],q)
        if UR>0 and LL>0: out.append(np.log(R))
    return np.array(out)

print("=== PUNTO 3: cociente de esquinas con IC bootstrap del 90% (5000 remuestras) ===")
print(f"{'par':22s}{'analisis':16s}{'R':>7s}{'IC 90% de R':>20s}{'  ¿compatible con R=1?'}")
for a,b in PAIRS:
    for tag,UU in [('agrupado',U0),('intra-estrato',Ua)]:
        u=UU[:,[idx[a],idx[b]]]
        LL,UR,R=ratio(u,.2); lr=boot_logR(u,.2)
        lo,hi=np.exp(np.percentile(lr,[5,95]))
        comp="si" if lo<=1<=hi else "NO"
        print(f"{a+'-'+b:22s}{tag:16s}{R:7.2f}   [{lo:5.2f}, {hi:5.2f}]      {comp}")
    print()

print("=== PUNTO 26: sensibilidad del cociente al umbral u ===")
print(f"{'par':22s}{'analisis':16s}" + "".join(f"{f'u={q}':>12s}" for q in (.05,.10,.15,.20)))
for a,b in PAIRS[:3]:
    for tag,UU in [('agrupado',U0),('intra-estrato',Ua)]:
        u=UU[:,[idx[a],idx[b]]]
        row=""
        for q in (.05,.10,.15,.20):
            LL,UR,R=ratio(u,q)
            row+=f"{R:12.2f}" if np.isfinite(R) else f"{'inf':>12s}"
        print(f"{a+'-'+b:22s}{tag:16s}{row}")
    print()
