import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys; sys.path.insert(0,'src')
import vinelib as V
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.linewidth':.6,'pdf.fonttype':42})
BLUE='#1565C0'; RED='#C62828'; GREY='#9E9E9E'; DARK='#37474F'
d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d)
cop=V.fit(V.pseudo(X,np.random.default_rng(20240607))); E=V.named_edges(cop)
POS={'Age':(-3.0,1.5),'Sex':(-3.3,0),'Education':(-3.0,-1.5),
     'Memory':(-.9,2.6),'ExecFun':(1.0,2.6),
     'APOE4':(-1.0,-2.7),'Amyloid':(.6,-2.7),'FDG':(2.2,-2.7),
     'PCC':(1.4,1.75),'Precuneus':(2.6,1.0),'Hippocampus':(4.0,1.75),
     'Thalamus':(4.6,.35),'Caudate':(4.6,-1.0),'Putamen':(3.4,-1.9)}
SHORT={'Age':'Age','Sex':'Sex','Education':'Educ.','Memory':'ADNI-MEM','ExecFun':'ADNI-EF',
  'APOE4':'APOE$\\varepsilon$4','Amyloid':'Amyloid','FDG':'FDG','Hippocampus':'Hippocamp.',
  'Thalamus':'Thalamus','Caudate':'Caudate','Putamen':'Putamen','PCC':'PCC','Precuneus':'Precuneus'}
CLR={'Age':'#5C6BC0','Sex':'#5C6BC0','Education':'#5C6BC0','Memory':'#43A047','ExecFun':'#43A047',
 'APOE4':'#FB8C00','Amyloid':'#FB8C00','FDG':'#FB8C00','Hippocampus':'#8E24AA','Thalamus':'#8E24AA',
 'Caudate':'#8E24AA','Putamen':'#8E24AA','PCC':'#8E24AA','Precuneus':'#8E24AA'}
DISC={'Sex','Education','APOE4'}; REFL={'Age','APOE4'}   # signo a invertir al mostrar
fig,ax=plt.subplots(figsize=(7.4,5.0))
for _,r in E.iterrows():
    a,b=r['a'],r['b']
    if a not in POS or b not in POS: continue
    x1,y1=POS[a]; x2,y2=POS[b]
    tau=r['tau']*(-1 if (a in REFL)^(b in REFL) else 1)   # signo en escala natural
    tail=r['lamL']>1e-6
    col=BLUE if tail else (RED if tau<0 else GREY)
    lw=.45+3.2*abs(tau); al=.95 if r['tree']==1 else .30
    ls='-' if r['tree']==1 else ('--' if r['tree']==2 else ':')
    ax.plot([x1,x2],[y1,y2],color=col,lw=lw,alpha=al,ls=ls,zorder=1,solid_capstyle='round')
    if tail and r['tree']==1:
        ax.text((x1+x2)/2,(y1+y2)/2,f"$\\kappa$={r['lamL']:.2f}",fontsize=7.2,color=BLUE,
                ha='center',va='center',zorder=5,fontweight='bold',
                bbox=dict(fc='white',ec=BLUE,lw=.6,pad=1.5,alpha=.97))
for k,(x,y) in POS.items():
    ax.scatter([x],[y],s=210,c=CLR[k],ec='white',lw=1.3,zorder=3)
    if k in DISC: ax.scatter([x],[y],s=400,facecolors='none',ec=DARK,lw=1.0,ls=(0,(2,1.5)),zorder=3)
    dy=.42 if y>=0 else -.42
    ax.text(x,y+dy,SHORT[k],fontsize=7.0,ha='center',va='bottom' if y>=0 else 'top',
            color=DARK,fontweight='bold',zorder=4)
leg=[Line2D([],[],color=BLUE,lw=2.4,label='lower-tail ($\\kappa>0$)'),
     Line2D([],[],color=GREY,lw=2.4,label='symmetric ($\\kappa=0$)'),
     Line2D([],[],color=RED,lw=2.4,label='negative ($\\tau<0$)'),
     Line2D([],[],color=DARK,lw=1.1,ls='-',label='$T_1$: unconditional'),
     Line2D([],[],color=DARK,lw=1.1,ls='--',label='$T_2$: given 1 var.'),
     Line2D([],[],color=DARK,lw=1.1,ls=':',label='$T_3$: given 2 var.'),
     Line2D([],[],marker='o',mfc='none',mec=DARK,ls='',ms=9,label='discrete margin ($\\kappa$ n.i.)')]
ax.legend(handles=leg,loc='upper center',bbox_to_anchor=(.5,.045),ncol=4,frameon=False,
          fontsize=6.8,handlelength=1.9,columnspacing=1.4)
ax.set_xlim(-4.4,5.6); ax.set_ylim(-4.1,3.5); ax.axis('off')
plt.tight_layout(); plt.savefig('figures/fig1_network.pdf',bbox_inches='tight'); plt.close(); print("fig1 v2 OK")
import numpy as np, pandas as pd, pickle, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys; sys.path.insert(0,'src')
import vinelib as V, pyvinecopulib as pv
from scipy.stats import norm, multivariate_normal as mvn
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.linewidth':.6,'pdf.fonttype':42})
BLUE='#1565C0'; GREY='#616161'; RED='#C62828'
d=pd.read_csv('data/adni_merged.csv'); X=V.orient(d); g=d.DX4.values
U=V.pseudo(X,np.random.default_rng(20240607)); idx={l:i for i,l in enumerate(V.LAB)}
def gbox(q,tau):
    r=np.sin(np.pi*tau/2); return mvn(mean=[0,0],cov=[[1,r],[r,1]]).cdf([norm.ppf(q),norm.ppf(q)])

# ---------------- FIG 2: dispersogramas ----------------
P=[('Memory','FDG','A','Clayton',BLUE),('Memory','Amyloid','B','Clayton',BLUE),
   ('Memory','Hippocampus','C','Frank',GREY),('ExecFun','Memory','D','Gaussian',GREY)]
fig,axs=plt.subplots(1,4,figsize=(7.4,2.25))
for ax,(a,b,pan,fam,col) in zip(axs,P):
    u=U[:,[idx[a],idx[b]]]; q=.2
    LL=np.mean((u[:,0]<q)&(u[:,1]<q)); UR=np.mean((u[:,0]>1-q)&(u[:,1]>1-q))
    tau=pv.Bicop.from_data(u,controls=pv.FitControlsBicop(family_set=[f for f in V.FAMS if f!=pv.BicopFamily.indep])).tau
    ax.add_patch(Rectangle((0,0),q,q,fc=col,alpha=.16,ec=col,lw=.7,zorder=0))
    ax.add_patch(Rectangle((1-q,1-q),q,q,fc='none',ec=col,lw=.7,ls='--',zorder=0))
    ax.scatter(u[:,0],u[:,1],s=2.2,c='#37474F',alpha=.42,lw=0,zorder=1)
    ax.text(.02,.185,f"{100*LL:.1f}%",fontsize=6.8,color=col,fontweight='bold',va='top')
    ax.text(.98,.815,f"{100*UR:.1f}%",fontsize=6.8,color=col,ha='right',va='bottom')
    ax.text(.5,-.30,f"ratio {LL/UR:.2f}   (Gauss {100*gbox(q,tau):.1f}%)",transform=ax.transAxes,
            ha='center',fontsize=6.6,color=col if LL/UR>1.5 else '#555555')
    lab=lambda z: z.replace('Memory','ADNI-MEM').replace('ExecFun','ADNI-EF').replace('Hippocampus','Hipp. vol.')
    ax.set_title(f"$\\bf{{{pan}}}$  {lab(a)}–{lab(b)}\n{fam}, $\\tau$={tau:.2f}",fontsize=7,pad=4)
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_xticks([0,.5,1]); ax.set_yticks([0,.5,1])
    ax.tick_params(labelsize=6); ax.set_aspect('equal')
axs[0].set_ylabel('$\\hat U_j$',fontsize=8)
fig.text(.5,-.02,'$\\hat U_i$  (small values = joint pathology)',ha='center',fontsize=8)
plt.tight_layout(); plt.savefig('figures/fig2_scatters.pdf',bbox_inches='tight'); plt.close(); print("fig2 OK")

# ---------------- FIG 4: agrupado vs intra-estrato ----------------
Ua=V.pseudo(X,np.random.default_rng(20240607),groups=g)
PR=[('Memory','FDG','ADNI-MEM – FDG'),('Memory','Amyloid','ADNI-MEM – Amyloid'),('PCC','Precuneus','PCC – Precuneus')]
us=np.linspace(.02,.20,25)
fig,axs=plt.subplots(2,3,figsize=(7.4,4.3),sharex=True,sharey=True)
rngB=np.random.default_rng(7)
for row,(UU,tag) in enumerate([(U,'Pooled sample'),(Ua,'Within-stratum ranks')]):
    for c,(a,b,ttl) in enumerate(PR):
        ax=axs[row,c]; u=UU[:,[idx[a],idx[b]]]; n=len(u)
        emp=np.array([np.mean((u[:,0]<t)&(u[:,1]<t))/t for t in us])
        Bs=np.array([[np.mean((u[i,0]<t)&(u[i,1]<t))/t for t in us]
                     for i in [rngB.integers(0,n,n) for _ in range(300)]])
        lo,hi=np.percentile(Bs,[5,95],axis=0)
        tau=pv.Bicop.from_data(u,controls=pv.FitControlsBicop(family_set=[f for f in V.FAMS if f!=pv.BicopFamily.indep])).tau
        r=np.sin(np.pi*tau/2)
        gau=np.array([mvn(mean=[0,0],cov=[[1,r],[r,1]]).cdf([norm.ppf(t),norm.ppf(t)])/t for t in us])
        ax.fill_between(us,lo,hi,color=BLUE,alpha=.18,lw=0)
        ax.plot(us,emp,color=BLUE,lw=1.7,label='$\\hat\\lambda_L(u)$')
        ax.plot(us,gau,color=RED,lw=1.3,ls='--',label='Gaussian, same $\\tau$')
        q=.2; LL=np.mean((u[:,0]<q)&(u[:,1]<q)); UR=np.mean((u[:,0]>1-q)&(u[:,1]>1-q))
        ax.text(.97,.94,f"ratio {LL/UR:.2f}",transform=ax.transAxes,ha='right',va='top',fontsize=7.5,
                fontweight='bold',color=BLUE if LL/UR>1.5 else '#666666',
                bbox=dict(fc='white',ec=BLUE if LL/UR>1.5 else '#BBBBBB',lw=.6,pad=1.8))
        if row==0: ax.set_title(ttl,fontsize=8)
        if c==0: ax.set_ylabel(f"{tag}\n$\\hat C_n(u,u)/u$",fontsize=7.5)
        if row==1: ax.set_xlabel('$u$',fontsize=8)
        ax.tick_params(labelsize=6.5); ax.set_ylim(0,.62)
axs[0,0].legend(fontsize=6.4,frameon=False,loc='lower right')
plt.tight_layout(); plt.savefig('figures/fig4_stratified.pdf',bbox_inches='tight'); plt.close(); print("fig4 OK")
import numpy as np, pickle, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.linewidth':.6,'pdf.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False})
BLUE='#1565C0'; GREY='#BDBDBD'; RED='#C62828'; DARK='#37474F'
R=pickle.load(open('output/boot_raw.pkl','rb')); B=len(R)
TAILF={'clayton','gumbel','joe','bb1','bb7','student'}
fig,axs=plt.subplots(2,3,figsize=(7.4,4.5))
TOP=[('Memory-FDG','A','ADNI-MEM – FDG'),('Memory-Amyloid','B','ADNI-MEM – Amyloid'),
     ('PCC-Precuneus','C','PCC – Precuneus vol.')]
for ax,(k,pan,ttl) in zip(axs[0],TOP):
    v=[r[k] for r in R]; t1=[x for x in v if x[0]==1 and x[1]==1]
    lam=np.array([x[4] for x in t1]); fam=[x[2] for x in t1]
    tail=np.array([f in TAILF for f in fam])
    xj=np.random.default_rng(3).normal(0,.10,len(lam))
    ax.axhline(0,color=RED,ls='--',lw=1.0,zorder=0)
    ax.scatter(xj[~tail],lam[~tail],s=3.5,c=GREY,alpha=.55,lw=0,zorder=1)
    ax.scatter(xj[tail],lam[tail],s=3.5,c=BLUE,alpha=.55,lw=0,zorder=2)
    lo,hi,md=np.percentile(lam,5),np.percentile(lam,95),np.median(lam)
    ax.plot([.42,.42],[lo,hi],color=DARK,lw=2.2,zorder=3,solid_capstyle='round')
    ax.scatter([.42],[md],s=34,facecolor='white',ec=DARK,lw=1.2,zorder=4)
    ax.text(.5,.965,f"$\\bf{{{pan}}}$  {ttl}",transform=ax.transAxes,ha='center',va='top',fontsize=7.4)
    ax.text(.5,.855,f"tail family in {100*tail.mean():.0f}% of $T_1$ resamples",transform=ax.transAxes,
            ha='center',va='top',fontsize=6.5,color=BLUE if tail.mean()>.9 else '#666666')
    ax.text(.5,.055,f"$\\kappa$ = {md:.2f}  [{lo:.2f}, {hi:.2f}]",transform=ax.transAxes,ha="center",
            va='bottom',fontsize=6.8,fontweight='bold')
    ax.set_xlim(-.5,.75); ax.set_ylim(-.04,.72); ax.set_xticks([])
    ax.tick_params(labelsize=6.5)
axs[0,0].set_ylabel('co-severity $\\kappa=\\lambda_L$',fontsize=7.6)
BOT=[('Memory-Hippocampus','D','ADNI-MEM – Hipp. (Frank)'),
     ('Age-Hippocampus','E','Age – Hipp. (Frank, neg.)'),
     ('Sex-Memory','F','Sex – ADNI-MEM (discrete)')]
for ax,(k,pan,ttl) in zip(axs[1],BOT):
    v=[r[k] for r in R]; sel=[x for x in v if x[0]==1]
    tau=np.array([x[3] for x in sel]); eip=len(sel)/B
    if k=='Age-Hippocampus': tau=-tau
    allt=np.concatenate([tau,np.zeros(B-len(sel))])
    ax.hist(allt,bins=34,color=BLUE if eip>.9 else '#90A4AE',alpha=.85,lw=0)
    ax.axvline(0,color=RED,ls='--',lw=.9)
    ax.text(.5,.965,f"$\\bf{{{pan}}}$  {ttl}",transform=ax.transAxes,ha='center',va='top',fontsize=7.2)
    ax.text(.5,.85,f"EIP = {eip:.2f}",transform=ax.transAxes,ha='center',va='top',fontsize=6.8,
            color=BLUE if eip>.9 else RED,fontweight='bold')
    ax.set_xlabel("Kendall's $\\tau$",fontsize=7.4); ax.tick_params(labelsize=6.5)
axs[1,0].set_ylabel('bootstrap resamples',fontsize=7.6)
plt.tight_layout(); plt.savefig('figures/fig3_bootstrap_kappa.pdf',bbox_inches='tight'); plt.close()
print("fig3 OK")
