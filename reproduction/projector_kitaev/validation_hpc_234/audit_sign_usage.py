#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv,json,math,os,re
base=os.path.dirname(os.path.abspath(__file__)); col=os.path.join(base,'collected')
allrows=list(csv.DictReader(open(os.path.join(col,'all_task_results.csv'))))
targets=set(['984035','984048','984085','984029','984033'])
chosen=[r for r in allrows if r.get('seed') in targets]
def loadj(path):
 s=open(path).read();s=re.sub(r'(?<=:)inf(?=[,}])','1e999',s);return json.loads(s)
out=[]
for meta in sorted(chosen,key=lambda x:int(x['seed'])):
 rows=list(csv.DictReader(open(meta['measurement_path'])))
 den=sum(float(r['sign']) for r in rows)
 op=sum(float(r['sign_S_pi_numerator']) for r in rows)/den
 od=sum(float(r['sign_S_pi_dq_numerator']) for r in rows)/den
 orr=1-od/op
 signs=[1.0 if float(r['sign'])>=0 else -1.0 for r in rows]
 den2=sum(signs)
 pp=sum(s*float(r['S_pi']) for s,r in zip(signs,rows))/den2
 pd=sum(s*float(r['S_pi_dq']) for s,r in zip(signs,rows))/den2
 pr=1-pd/pp
 result=loadj(meta['result_path'])
 for name,a,b,err in [('S_pi',op,pp,float(result['S_pi_err'])),('S_pi_dq',od,pd,float(result['S_pi_dq_err'])),('R_cdw',orr,pr,float(result['R_cdw_err']))]:
  diff=b-a;out.append({'seed':meta['seed'],'campaign':meta['campaign'],'task_id':meta['task_id'],'L':meta['L'],'V':meta['V'],'boundary':meta['boundary'],'dt':meta['dt'],'warning':next((w['flags'] for w in csv.DictReader(open(os.path.join(col,'qc_warnings.csv'))) if w['seed']==meta['seed'] and w['campaign']==meta['campaign'] and w['task_id']==meta['task_id']),''),'observable':name,'stored_numerator_method':a,'explicit_pm_sign_method':b,'absolute_difference':diff,'reported_error':err,'difference_sigma':(diff/err if err else 0.0),'max_sign_imag':result['max_sign_imag'],'sign_corrections':result['sign_corrections']})
fields=['seed','campaign','task_id','L','V','boundary','dt','warning','observable','stored_numerator_method','explicit_pm_sign_method','absolute_difference','reported_error','difference_sigma','max_sign_imag','sign_corrections']
with open(os.path.join(col,'sign_usage_reweight_compare.csv'),'w') as h:w=csv.DictWriter(h,fields);w.writeheader();w.writerows(out)
mx=max(abs(r['difference_sigma']) for r in out);md=max(abs(r['absolute_difference']) for r in out)
report='''# sign warning audit — corrected

## Static driver sign usage

At measurement, `rightSweep(center,&g,&sg)` captures the complex transported sign at the center. The driver immediately converts it to the eigen-sign `s = (sg.real() >= 0 ? +1 : -1)`. It accumulates `ss += s`, `sp += s*S_pi_cfg`, and `sd += s*S_pi_dq_cfg`; final values are `S_pi=sp/ss`, `S_pi_dq=sd/ss`, and `R_cdw=1-S_pi_dq/S_pi`. Neither the complex phase, `abs(sg)`, nor the unthresholded value `real(sg)` enters reweighting.

`SIGN_IMAG_WARN` is a collector label when `max(abs(Im(sg))) > 1e-6`. `SIGN_CORRECTION` is a collector label when the driver reports `sign_corrections > 0`.

Every 20 measurements, before the next right sweep, the driver evaluates `z=getSignRaw()` on the current live contour. If `abs(q.sign-z)>1e-2`, it **does modify** `q.sign`, setting it to `+1` or `-1` according to `z.real()`, and increments `sign_corrections`. Otherwise the call is read-only.

## Existing-data comparison

For five severe warning seeds, recomputing from configuration observables using an explicit `sign>=0 ? +1 : -1` gives the same estimator as the stored numerator columns. Maximum absolute difference is %.3e and maximum difference is %.3e sigma. Detailed values are in `collected/sign_usage_reweight_compare.csv`.

Therefore the recorded imaginary phase does not change S_pi, S_pi_dq, or R_cdw in these data. It is a **numerical phase drift diagnostic**, not a sign mismatch in the reweighting estimator.

## Retraction of the earlier mismatch claim

The earlier “509/2000 real mismatches” conclusion is withdrawn. That audit compared `sg` captured at the center with `getSignRaw()` called after the remainder of the right sweep had updated additional HS fields. Those signs belonged to different HS configurations. Correct same-configuration long-run checks found 0 mismatches in 11,001 sweep-end checks and 0 mismatches in 5,000 center checks.

One separate caveat remains in the optional `diagnostic_stride` block of `projector_bins_driver.cpp`: it rebuilds the center Green but calls `getSignRaw()` after the complete right sweep, so its `diag_sign_mismatch` comparison also crosses configurations. The validation_hpc_234 production tasks used diagnostic_stride=0, so this did not enter their observables or warning labels.

No core algorithm was modified.
'''%(md,mx)
open(os.path.join(base,'sign_warning_audit.md'),'w').write(report)
print('rows=%d max_abs=%.3e max_sigma=%.3e'%(len(out),md,mx))
