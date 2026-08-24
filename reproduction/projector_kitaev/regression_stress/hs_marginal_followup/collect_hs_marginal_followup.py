#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division,print_function
import csv,json,math,os,re
from collections import defaultdict

BASE=os.path.dirname(os.path.abspath(__file__));BIN=25;OBS=('S_pi','S_pi_dq','R_cdw')
def loadj(p):
 s=open(p).read();s=re.sub(r'(?<=:)inf(?=[,}])','1e999',s);s=re.sub(r'(?<=:)-inf(?=[,}])','-1e999',s);s=re.sub(r'(?<=:)nan(?=[,}])','null',s);return json.loads(s)
def mean(x):return sum(x)/len(x)
def tau(x):
 n=len(x);m=mean(x);v=sum((z-m)**2 for z in x)/n
 if v<=0:return .5
 ans=.5
 for lag in range(1,min(n//2,2000)):
  c=sum((x[i]-m)*(x[i+lag]-m) for i in range(n-lag))/(n-lag)/v
  if c<=0:break
  ans+=c
 return ans
def ob(total,name):
 s,p,d=total
 if name=='S_pi':return p/s
 if name=='S_pi_dq':return d/s
 return 1-d/p
def jkerr(x):
 m=mean(x);return math.sqrt((len(x)-1)/len(x)*sum((z-m)**2 for z in x))
def write(p,rows,fields):
 with open(p,'w') as h:w=csv.DictWriter(h,fields);w.writeheader();w.writerows(rows)
def main():
 manifest=list(csv.DictReader(open(os.path.join(BASE,'manifest.csv'))));data={};meta={};missing=[]
 for m in manifest:
  out=os.path.join(BASE,m['output_dir']);rp=os.path.join(out,'result.json');mp=os.path.join(out,'measurements.csv')
  if not os.path.isfile(rp) or not os.path.isfile(mp):missing.append(m['task_id']);continue
  rows=list(csv.DictReader(open(mp)));vals=[(float(x['sign']),float(x['sign_S_pi_numerator']),float(x['sign_S_pi_dq_numerator'])) for x in rows]
  data[(int(m['hs_scheme']),int(m['seed']))]=vals;meta[(int(m['hs_scheme']),int(m['seed']))]=loadj(rp)
 if missing:
  open(os.path.join(BASE,'hs_marginal_followup.md'),'w').write('# HS marginal follow-up\n\nRunning: %d/24 tasks complete; missing task IDs: %s.\n'%(24-len(missing),','.join(missing)));print('missing',missing);return
 seeds=sorted(set(seed for hs,seed in data));units=defaultdict(dict);tot={0:[0.,0.,0.],1:[0.,0.,0.]};taurows=[]
 for hs in (0,1):
  for seed in seeds:
   x=data[(hs,seed)]
   for name,col in [('sign',0),('sign_S_pi',1),('sign_S_pi_dq',2)]:taurows.append({'hs_scheme':hs,'seed':seed,'series':name,'tau_int_retained_measurements':tau([z[col] for z in x])})
   for b in range(len(x)//BIN):
    q=x[b*BIN:(b+1)*BIN];block=[sum(z[i] for z in q) for i in range(3)];units[(seed,b)][hs]=block
    for i in range(3):tot[hs][i]+=block[i]
 keys=sorted(units);ed=loadj(os.path.join(BASE,'ed','result.json'));results=[]
 for name in OBS:
  full=[ob(tot[hs],name) for hs in (0,1)];rep=[[],[]];diff=[]
  for k in keys:
   v=[]
   for hs in (0,1):v.append(ob([tot[hs][i]-units[k][hs][i] for i in range(3)],name));rep[hs].append(v[-1])
   diff.append(v[0]-v[1])
  e=[jkerr(rep[0]),jkerr(rep[1])];de=jkerr(diff);dv=full[0]-full[1]
  results.append({'observable':name,'hs0':full[0],'hs0_error':e[0],'hs1':full[1],'hs1_error':e[1],'hs0_minus_hs1':dv,'difference_error':de,'difference_sigma':dv/de,'ED':ed[name],'hs0_minus_ED':full[0]-ed[name],'hs0_minus_ED_sigma':(full[0]-ed[name])/e[0],'hs1_minus_ED':full[1]-ed[name],'hs1_minus_ED_sigma':(full[1]-ed[name])/e[1],'jackknife_bin_size':BIN,'jackknife_units':len(keys)})
 write(os.path.join(BASE,'hs_marginal_followup.csv'),results,['observable','hs0','hs0_error','hs1','hs1_error','hs0_minus_hs1','difference_error','difference_sigma','ED','hs0_minus_ED','hs0_minus_ED_sigma','hs1_minus_ED','hs1_minus_ED_sigma','jackknife_bin_size','jackknife_units'])
 write(os.path.join(BASE,'tau_int.csv'),taurows,['hs_scheme','seed','series','tau_int_retained_measurements'])
 diagnostics=[]
 for hs in (0,1):
  ms=[meta[(hs,s)] for s in seeds];diagnostics.append({'hs_scheme':hs,'average_sign':mean([x['average_sign'] for x in ms]),'acceptance':mean([x['acceptance'] for x in ms]),'max_tau_int':max(x['tau_int_retained_measurements'] for x in taurows if x['hs_scheme']==hs),'direct_comparisons':sum(x['diagnostic_comparisons'] for x in ms),'pm1_sign_mismatches':sum(x['diagnostic_sign_mismatch_count'] for x in ms),'sign_corrections':sum(x['sign_corrections'] for x in ms),'max_sign_imag':max(x['max_sign_imag'] for x in ms),'max_fast_full_green_relative_error':max(x['diagnostic_relative_frobenius_max'] for x in ms)})
 write(os.path.join(BASE,'diagnostics.csv'),diagnostics,['hs_scheme','average_sign','acceptance','max_tau_int','direct_comparisons','pm1_sign_mismatches','sign_corrections','max_sign_imag','max_fast_full_green_relative_error'])
 sdq=next(x for x in results if x['observable']=='S_pi_dq');reproduced=abs(sdq['difference_sigma'])>=3
 lines=['# HS marginal follow-up\n\n','Parameters: L=10, V=2, theta=18, beta_trial=8, dt=0.1, PBC; guard OFF; 12 new matched seeds and 10000 measurements per scheme. Common delete-one `(seed, 25-measurement bin)` jackknife units: %d.\n\n'%len(keys),'## Observables\n\n','| observable | hs0 | hs1 | hs0-hs1 | sigma | ED | hs0-ED sigma | hs1-ED sigma |\n|---|---:|---:|---:|---:|---:|---:|---:|\n']
 for x in results:lines.append('| %s | %.9g ± %.2g | %.9g ± %.2g | %.3g ± %.2g | %+.2f | %.9g | %+.2f | %+.2f |\n'%(x['observable'],x['hs0'],x['hs0_error'],x['hs1'],x['hs1_error'],x['hs0_minus_hs1'],x['difference_error'],x['difference_sigma'],x['ED'],x['hs0_minus_ED_sigma'],x['hs1_minus_ED_sigma']))
 lines+=['\n## Diagnostics\n\n']
 for x in diagnostics:lines.append('* hs%d: average sign %.6g; acceptance %.6g; max tau_int %.3g retained measurements; direct comparisons %d; ±1 mismatches %d; sign corrections %d; max sign imaginary %.3e; max fast/full Green relative error %.3e.\n'%(x['hs_scheme'],x['average_sign'],x['acceptance'],x['max_tau_int'],x['direct_comparisons'],x['pm1_sign_mismatches'],x['sign_corrections'],x['max_sign_imag'],x['max_fast_full_green_relative_error']))
 lines+=['\n## Verdict\n\n','* Original S_pi_dq 3.43 sigma difference %s.\n'%('reproduces' if reproduced else 'does not reproduce')]
 if reproduced:lines.append('* ED deviations above identify which representation is displaced.\n')
 else:lines.append('* The previous marginal difference is consistent with a statistical multiple-comparison fluctuation.\n')
 if any(x['pm1_sign_mismatches'] for x in diagnostics) or any(x['max_fast_full_green_relative_error']>1e-6 for x in diagnostics):lines.append('* Sign/Green diagnostics show an anomaly requiring follow-up.\n')
 else:lines.append('* No ±1 sign mismatch or material fast/full Green anomaly is present.\n')
 open(os.path.join(BASE,'hs_marginal_followup.md'),'w').write(''.join(lines));print('complete')
if __name__=='__main__':main()
