#!/usr/bin/env python
from __future__ import division
import csv, os
from collections import defaultdict

base=os.path.dirname(os.path.abspath(__file__))
allrows={(r['campaign'],r['task_id']):r for r in csv.DictReader(open(os.path.join(base,'collected','all_task_results.csv')))}
warnings=list(csv.DictReader(open(os.path.join(base,'collected','qc_warnings.csv'))))
out=[]; groups=defaultdict(list)
for w in warnings:
    r=allrows[(w['campaign'],w['task_id'])]
    x=dict(w)
    for k in ('L','V','boundary','dt','theta','beta_trial','burn','measurements','average_sign','result_path'):
        x[k]=r[k]
    out.append(x)
    groups[(x['campaign'],x['L'],x['V'],x['boundary'],x['dt'],x['flags'])].append(x)
fields=['campaign','task_id','seed','L','V','boundary','dt','theta','beta_trial','burn','measurements','flags','max_sign_imag','sign_corrections','diagnostic_sign_mismatch_count','average_sign','result_path']
with open(os.path.join(base,'collected','qc_warnings_enriched.csv'),'w') as h:
    w=csv.DictWriter(h,fields,extrasaction='ignore');w.writeheader();w.writerows(out)
summary=[]
for k,vs in sorted(groups.items()):
    summary.append({'campaign':k[0],'L':k[1],'V':k[2],'boundary':k[3],'dt':k[4],'warning_type':k[5],'n':len(vs),'seeds':';'.join(x['seed'] for x in vs),'max_sign_imag':max(float(x['max_sign_imag']) for x in vs),'sum_corrections':sum(int(x['sign_corrections']) for x in vs)})
with open(os.path.join(base,'collected','qc_warning_parameter_summary.csv'),'w') as h:
    w=csv.DictWriter(h,['campaign','L','V','boundary','dt','warning_type','n','seeds','max_sign_imag','sum_corrections']);w.writeheader();w.writerows(summary)
# Representative policy: top three imaginary amplitudes, and two correction cases.
top=sorted(out,key=lambda x:float(x['max_sign_imag']),reverse=True)[:3]
corr=[x for x in out if 'SIGN_CORRECTION' in x['flags']]
chosen=[]
for x in top+corr[:2]:
    if (x['campaign'],x['task_id']) not in [(y['campaign'],y['task_id']) for y in chosen]: chosen.append(x)
for i,x in enumerate(chosen): x['audit_id']=i; x['audit_measurements']=400; x['audit_burn']=100; x['direct_check_stride']=1
with open(os.path.join(base,'manifests','warning_audit.csv'),'w') as h:
    w=csv.DictWriter(h,['audit_id']+fields+['audit_burn','audit_measurements','direct_check_stride'],extrasaction='ignore');w.writeheader();w.writerows(chosen)
print('warnings=%d groups=%d audit_tasks=%d' % (len(out),len(summary),len(chosen)))
