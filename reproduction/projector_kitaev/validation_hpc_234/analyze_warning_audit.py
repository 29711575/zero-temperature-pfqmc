#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import csv, glob, json, math, os, re

base=os.path.dirname(os.path.abspath(__file__))
def loadj(p):
 s=open(p).read(); s=re.sub(r'(?<=:)inf(?=[,}])','1e999',s); return json.loads(s)
def num(x):
 try:return float(x)
 except:return float('nan')
def avg(x): return sum(x)/len(x) if x else float('nan')

orig={(r['campaign'],r['task_id']):r for r in csv.DictReader(open(os.path.join(base,'collected','qc_warnings_enriched.csv')))}
out=[]
for path in sorted(glob.glob(os.path.join(base,'raw','warning_audit','*_retry1','result.json'))):
 d=loadj(path); folder=os.path.basename(os.path.dirname(path)); m=re.search(r'task_(\d+)_from_(.+)_(\d+)_s(\d+)_retry1',folder)
 aid,campaign,tid,seed=m.group(1),m.group(2),m.group(3),m.group(4)
 rows=list(csv.DictReader(open(os.path.join(os.path.dirname(path),'measurements.csv'))))
 mm=[r for r in rows if int(r['transport_direct_mismatch'])==1]
 ok=[r for r in rows if int(r['transport_direct_mismatch'])==0]
 def stat(rs,k):
  z=[num(r[k]) for r in rs]; return (min(z),max(z),avg(z)) if z else (float('nan'),)*3
 omin=orig[(campaign,tid)]
 x={'audit_id':aid,'campaign':campaign,'task_id':tid,'seed':seed,'L':d['L'],'V':d['V'],'boundary':d['boundary'],'dt':d['dt'],'original_warning':omin['flags'],'original_max_sign_imag':omin['max_sign_imag'],'original_sign_corrections':omin['sign_corrections'],'audit_burn':d['burn'],'audit_measurements':d['measurements'],'mismatch_count':len(mm),'mismatch_fraction':len(mm)/float(len(rows)),'direct_checks':d['transport_direct_comparisons'],'direct_sign_imag_max':d['max_direct_sign_imag'],'tracked_sign_imag_max':d['max_sign_imag'],'audit_sign_corrections':d['sign_corrections'],'average_sign':d['average_sign'],'full_contour_rel_green_max':d['diagnostic_relative_frobenius_max'],'full_contour_Spi_abs_diff_max':d['diagnostic_S_pi_abs_diff_max'],'full_contour_R_abs_diff_max':d['diagnostic_R_cdw_abs_diff_max']}
 for prefix,rs in [('all',rows),('mismatch',mm),('matched',ok)]:
  for key,label in [('S_pi','S_pi'),('S_pi_dq','S_pi_dq'),('R_cdw','R_cdw')]:
   lo,hi,mu=stat(rs,key); x[prefix+'_'+label+'_min']=lo; x[prefix+'_'+label+'_max']=hi; x[prefix+'_'+label+'_mean']=mu
 out.append(x)
fields=['audit_id','campaign','task_id','seed','L','V','boundary','dt','original_warning','original_max_sign_imag','original_sign_corrections','audit_burn','audit_measurements','direct_checks','mismatch_count','mismatch_fraction','direct_sign_imag_max','tracked_sign_imag_max','audit_sign_corrections','average_sign','full_contour_rel_green_max','full_contour_Spi_abs_diff_max','full_contour_R_abs_diff_max']
for p in ('all','mismatch','matched'):
 for n in ('S_pi','S_pi_dq','R_cdw'): fields += [p+'_'+n+'_min',p+'_'+n+'_max',p+'_'+n+'_mean']
with open(os.path.join(base,'collected','sign_warning_audit.csv'),'w') as h:
 w=csv.DictWriter(h,fields);w.writeheader();w.writerows(out)
total=sum(x['mismatch_count'] for x in out); checks=sum(x['direct_checks'] for x in out)
lines=['# sign-warning audit\n','## Original warning concentration\n','The 28 warnings are all at dt=0.1: L=10,V=4,PBC (5 SIGN_IMAG_WARN); L=6,V=4/6, with PBC dominant. None is from dt=0.2. Detailed grouping: `collected/qc_warning_parameter_summary.csv`.\n','## Targeted reruns\n','Five representative same-parameter/same-seed runs used guard OFF, 100 burn, 400 measurements, direct raw-sign stride=1, and full-contour Green diagnostic stride=1. `getSignRaw()` recomputes the complete contour Pfaffian; it is compared at the captured center boundary with the transported `q.sign`.\n','## Result: real transported/direct sign mismatches\n','**%d/%d (%.2f%%) measurement points have a ± transported-sign vs direct-raw-sign mismatch.** This is not classifiable as a harmless phase-only diagnostic. No core code was changed.\n' % (total,checks,100.0*total/checks),'The three largest original sign-imag cases reproduce the same scale in the audit: seed 984035: 5.793e-3 -> %.3e; 984048: 6.031e-4 -> %.3e; 984085: 1.265e-4 -> %.3e.\n' % (out[0]['tracked_sign_imag_max'],out[1]['tracked_sign_imag_max'],out[2]['tracked_sign_imag_max']),'Direct raw-sign imaginary parts remain much smaller (see CSV), while tracked max imaginary parts can be large. Configuration-level S_pi/S_pi_dq/R_cdw ranges and full-contour Green/observable discrepancies are retained in `collected/sign_warning_audit.csv`; they become large near mismatch episodes, so no physical interpretation is assigned.\n','## Per-seed\n']
for x in out:
 lines.append('* seed %s (L=%s,V=%s,BC=%s): mismatch %d/%d (%.1f%%); tracked/direct max imaginary %.3e / %.3e; original %s.\n' % (x['seed'],x['L'],x['V'],x['boundary'],x['mismatch_count'],x['direct_checks'],100*x['mismatch_fraction'],x['tracked_sign_imag_max'],x['direct_sign_imag_max'],x['original_warning']))
open(os.path.join(base,'sign_warning_audit.md'),'w').write(''.join(lines))
print('audit tasks=%d mismatch=%d/%d' % (len(out),total,checks))
