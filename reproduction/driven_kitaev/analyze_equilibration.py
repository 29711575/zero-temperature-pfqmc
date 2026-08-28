#!/usr/bin/env python3
"""Offline burn-in and autocorrelation analysis for driven PfQMC debug runs."""
import csv,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'/'equilibration_check'
EDFILE=HERE/'results'/'ed_benchmark_v2'/'ed_results.csv'
OBS=('S_pi','S_pi_dq','R_cdw')
SERIES=('sign','sign_S_pi_numerator','sign_S_pi_dq_numerator')

def meanerr(x):
    x=np.asarray(x,float); return float(x.mean()),float(x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else math.nan

def autocorr_tau(x):
    x=np.asarray(x,float); n=len(x); x=x-x.mean()
    if n<20 or np.dot(x,x)==0:return np.array([1.]),.5,float(n),False
    f=np.fft.rfft(x,n=2*n); ac=np.fft.irfft(f*np.conjugate(f))[:n]
    ac/=np.arange(n,0,-1); ac/=ac[0]
    tau=.5; cutoff=1
    for lag in range(1,min(n//2,1000)):
        if ac[lag]<=0:break
        tau+=ac[lag];cutoff=lag+1
        if lag>5*tau:break
    stable=cutoff>2 and cutoff<min(n//2,1000)-1
    return ac[:min(201,n)],float(tau),float(n/(2*tau)),stable

def write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    manifest=list(csv.DictReader(open(HERE/'manifest_equilibration.csv')))
    ed={}
    for x in csv.DictReader(open(EDFILE)):
        if int(x['L'])==6 and float(x['V0'])==0 and float(x['rate'])==1:
            ed[(float(x['Vf']),float(x['dt']))]={o:float(x[o]) for o in OBS}
    raw=[]; ts={}; groups=defaultdict(list); acrows=[]; blockrows=[]; cumrows=[]
    for m in manifest:
        d=HERE/m['output_dir']; result=d/'result.json'; series=d/'timeseries.csv'
        if not result.exists() or not series.exists():raise SystemExit(f'missing {d}')
        x=json.load(open(result)); key=(m['label'],float(m['Vf']),float(m['dt']),int(m['burn']))
        groups[key].append(x); raw.append(x)
        a=np.genfromtxt(series,delimiter=',',names=True); ts[x['seed']]=a
        for name in SERIES:
            ac,tau,ess,stable=autocorr_tau(a[name])
            acrows.append(dict(label=m['label'],Vf=m['Vf'],dt=m['dt'],burn=m['burn'],seed=x['seed'],series=name,tau_int=tau,ESS=ess,stable=stable))
            for lag,v in enumerate(ac):
                if lag<=100: blockrows.append(dict(record='autocorrelation',label=m['label'],Vf=m['Vf'],dt=m['dt'],burn=m['burn'],seed=x['seed'],series=name,index=lag,value=float(v)))
        n=len(a); bs=300
        for start in range(0,n,bs):
            sl=slice(start,min(start+bs,n)); den=a['sign'][sl].sum()
            if den==0:continue
            p=a['sign_S_pi_numerator'][sl].sum()/den; q=a['sign_S_pi_dq_numerator'][sl].sum()/den
            for name,val in [('S_pi',p),('S_pi_dq',q),('R_cdw',1-q/p)]:blockrows.append(dict(record='block_mean',label=m['label'],Vf=m['Vf'],dt=m['dt'],burn=m['burn'],seed=x['seed'],series=name,index=start,value=float(val)))
        den=np.cumsum(a['sign']); p=np.divide(np.cumsum(a['sign_S_pi_numerator']),den,out=np.full(n,np.nan),where=den!=0);q=np.divide(np.cumsum(a['sign_S_pi_dq_numerator']),den,out=np.full(n,np.nan),where=den!=0)
        for i in range(99,n,100):
            for name,val in [('S_pi',p[i]),('S_pi_dq',q[i]),('R_cdw',1-q[i]/p[i])]:cumrows.append(dict(label=m['label'],Vf=m['Vf'],dt=m['dt'],burn=m['burn'],seed=x['seed'],measurement=i,value=float(val),observable=name))
    summary=[]
    for (label,vf,dt,burn),xs in sorted(groups.items()):
        row=dict(label=label,L=6,V0=0,Vf=vf,rate=1,dt=dt,burn=burn,n_seeds=len(xs))
        for o in OBS:
            row[o],row[o+'_stderr']=meanerr([x[o] for x in xs]);row['ED_'+o]=ed[(vf,dt)][o];row['z_'+o]=(row[o]-row['ED_'+o])/row[o+'_stderr']
        for o in ('average_sign','acceptance'):row[o],row[o+'_stderr']=meanerr([x[o] for x in xs])
        for name in SERIES:
            vals=[r for r in acrows if r['label']==label and float(r['Vf'])==vf and float(r['dt'])==dt and int(r['burn'])==burn and r['series']==name]
            row['tau_'+name]=float(np.median([r['tau_int'] for r in vals]));row['ESS_'+name]=float(sum(r['ESS'] for r in vals));row['tau_'+name+'_stable']=all(r['stable'] for r in vals)
        summary.append(row)
    write(OUT/'burn_summary.csv',summary);write(OUT/'autocorrelation_summary.csv',acrows);write(OUT/'diagnostic_series.csv',blockrows);write(OUT/'cumulative_means.csv',cumrows)
    try:
        import matplotlib.pyplot as plt
        for label in ('A','B','C'):
            ss=[r for r in summary if r['label']==label]
            if not ss:continue
            fig,ax=plt.subplots(1,3,figsize=(10,3))
            for a,o in zip(ax,OBS):
                a.errorbar([r['burn'] for r in ss],[r[o] for r in ss],yerr=[r[o+'_stderr'] for r in ss],fmt='o-');a.axhline(ss[0]['ED_'+o],color='k',ls='--');a.set(xlabel='burn',title=o)
            fig.tight_layout();fig.savefig(OUT/f'observable_vs_burn_{label}.png',dpi=120);plt.close(fig)
        for label in ('A','B'):
            chosen=[m for m in manifest if m['label']==label and int(m['burn']) in (300,6000) and m['seed'].endswith('1')]
            fig,ax=plt.subplots(1,3,figsize=(10,3))
            for m in chosen:
                a=ts[int(m['seed'])];den=np.cumsum(a['sign']);p=np.divide(np.cumsum(a['sign_S_pi_numerator']),den,out=np.full(len(a),np.nan),where=den!=0);q=np.divide(np.cumsum(a['sign_S_pi_dq_numerator']),den,out=np.full(len(a),np.nan),where=den!=0)
                for aa,o,y in zip(ax,OBS,(p,q,1-q/p)):aa.plot(a['measurement'],y,label='burn='+m['burn']);aa.axhline(ed[(float(m['Vf']),float(m['dt']))][o],color='k',ls='--');aa.set(title=o,xlabel='measurement')
            ax[0].legend();fig.tight_layout();fig.savefig(OUT/f'cumulative_{label}.png',dpi=120);plt.close(fig)
        fig,ax=plt.subplots(1,3,figsize=(10,3))
        for aa,name in zip(ax,SERIES):
            for label in ('A','B'):
                r=next(m for m in manifest if m['label']==label and int(m['burn'])==6000);ac,*_=autocorr_tau(ts[int(r['seed'])][name]);aa.plot(ac[:101],label=label)
            aa.set(title=name,xlabel='lag');aa.legend()
        fig.tight_layout();fig.savefig(OUT/'autocorrelation.png',dpi=120);plt.close(fig)
    except ImportError:pass

if __name__=='__main__':main()
