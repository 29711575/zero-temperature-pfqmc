#!/usr/bin/env python3
import csv,json,math,os
from collections import defaultdict,Counter
ROOT='/home/sunxr/PfQMC-main/reproduction/driven_kitaev/results/driven_critical_rate_ratio'
OUT=ROOT+'/collected'; os.makedirs(OUT,exist_ok=True)
def txt(p):
 try:
  with open(p) as f:return f.read()
 except:return ''
def load(p):return json.loads(txt(p).replace(':inf',':null').replace(':-inf',':null'))
manifest=list(csv.DictReader(open(ROOT+'/manifests/rigorous_manifest.csv')))
rows=[]
for m in manifest:
 d=ROOT+'/'+m['output_dir']; st=txt(d+'/status.txt')
 if 'exit_code=3' not in st:continue
 x=load(d+'/result.json.tmp'); tol=float(x.get('sign_imag_tolerance',1e-7))
 si=float(x.get('max_sign_imag',0));oi=float(x.get('max_observable_imag',0));ge=float(x.get('max_green_skew_symmetry_error',0));sc=int(x.get('sign_corrections',0))
 flags=[]
 if si>tol:flags.append('SIGN_IMAG')
 if oi>tol:flags.append('OBSERVABLE_IMAG')
 if ge>tol:flags.append('GREEN_SKEW')
 if sc>0:flags.append('SIGN_CORRECTION')
 if not flags:flags=['OTHER_DIAGNOSTIC']
 ratios={'sign':si/tol,'observable':oi/tol,'green':ge/tol}
 rows.append({'task_id':int(m['task_id']),'L':int(m['L']),'Vf':float(m['Vf']),'rate':float(m['rate']),'seed':int(m['seed']),'failure_reason':x.get('failure_reason',''),'trigger_flags':';'.join(flags),'tolerance':tol,'max_sign_imag':si,'sign_over_tol':ratios['sign'],'max_observable_imag':oi,'observable_over_tol':ratios['observable'],'max_green_skew_symmetry_error':ge,'green_over_tol':ratios['green'],'sign_corrections':sc,'average_sign':x.get('average_sign',''),'acceptance':x.get('acceptance',''),'measurements_completed':x.get('measurements_completed',''),'max_excess_ratio':max(ratios.values()),'result_tmp':d+'/result.json.tmp'})
fields=list(rows[0])
with open(OUT+'/driven_failure_audit.csv','w',newline='') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(sorted(rows,key=lambda r:r['task_id']))
groups=defaultdict(list)
for r in rows:groups[(r['L'],r['Vf'],r['rate'])].append(r)
gf=['L','Vf','rate','planned','failures','failure_rate','sign_imag_count','observable_imag_count','green_skew_count','sign_correction_count','other_count','max_sign_over_tol','max_observable_over_tol','max_green_over_tol','median_max_excess_ratio']
grows=[]
for k,rs in sorted(groups.items()):
 vals=sorted(r['max_excess_ratio'] for r in rs)
 grows.append(dict(zip(gf[:3],k),planned=12,failures=len(rs),failure_rate=len(rs)/12,sign_imag_count=sum('SIGN_IMAG' in r['trigger_flags'] for r in rs),observable_imag_count=sum('OBSERVABLE_IMAG' in r['trigger_flags'] for r in rs),green_skew_count=sum('GREEN_SKEW' in r['trigger_flags'] for r in rs),sign_correction_count=sum(r['sign_corrections']>0 for r in rs),other_count=sum(r['trigger_flags']=='OTHER_DIAGNOSTIC' for r in rs),max_sign_over_tol=max(r['sign_over_tol'] for r in rs),max_observable_over_tol=max(r['observable_over_tol'] for r in rs),max_green_over_tol=max(r['green_over_tol'] for r in rs),median_max_excess_ratio=vals[len(vals)//2]))
with open(OUT+'/driven_failure_audit_by_LVR.csv','w',newline='') as f:w=csv.DictWriter(f,gf);w.writeheader();w.writerows(grows)
# Six representative missing compatible tasks, R=0.5, covering L=26/34 and V=3.8/4.0/4.2; force earliest and global maximum into selection.
c=[r for r in rows if r['rate']==0.5 and r['L'] in (26,34) and r['Vf'] in (3.8,4.0,4.2)]
selected=[]
def add(r,why):
 if r and r['task_id'] not in {z['task_id'] for z in selected}: selected.append(dict(r,selection_reason=why))
add(min(c,key=lambda r:r['task_id']),'earliest_compatible_failure')
add(max(c,key=lambda r:r['max_excess_ratio']),'largest_compatible_excess')
for L,V in ((26,3.8),(26,4.0),(26,4.2),(34,3.8),(34,4.0),(34,4.2)):
 if len(selected)>=6:break
 z=[r for r in c if r['L']==L and r['Vf']==V]
 add(max(z,key=lambda r:r['max_excess_ratio']) if z else None,f'representative_L{L}_V{V}')
sf=list(selected[0])
with open(OUT+'/driven_failure_diagnostic_selection.csv','w',newline='') as f:w=csv.DictWriter(f,sf);w.writeheader();w.writerows(selected)
trig=Counter()
for r in rows:
 for x in r['trigger_flags'].split(';'):trig[x]+=1
byrate=Counter((r['rate']) for r in rows)
with open(OUT+'/driven_failure_audit.md','w') as f:
 f.write('# Driven rigorous 153416 failure audit\n\n')
 f.write(f'Exit-code-3 tasks audited: {len(rows)}. Trigger counts are nonexclusive: {dict(trig)}. Failures by rate: {dict(sorted(byrate.items()))}.\n\n')
 f.write('A task is assigned a trigger when its recorded maximum exceeds its own JSON tolerance (normally 1e-7). Sign corrections are reported separately; they are not part of the driver exit-code-3 predicate.\n\n')
 f.write('## Representative deterministic reruns\n\nSelected tasks are listed in `driven_failure_diagnostic_selection.csv`. Same-configuration/same-boundary direct-sign and full-Green diagnostics are pending targeted rerun; no cross-configuration sign comparison is permitted.\n')
print('failures',len(rows),'triggers',dict(trig),'byrate',dict(byrate),'selected',[(r['task_id'],r['L'],r['Vf'],r['seed'],r['selection_reason']) for r in selected])
