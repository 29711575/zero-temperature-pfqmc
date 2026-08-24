#!/usr/bin/env python3
import csv, json, math, os, glob
from collections import defaultdict, Counter

ROOT='/home/sunxr/PfQMC-main/reproduction/driven_kitaev/results/driven_critical_rate_ratio'
MANIFEST=ROOT+'/manifests/rigorous_manifest.csv'
OUT=ROOT+'/collected'
os.makedirs(OUT,exist_ok=True)

def text(path):
    try:
        with open(path) as f: return f.read()
    except OSError: return ''
def jload(path):
    return json.loads(text(path).replace(':inf', ':null').replace(':-inf', ':null'))
def finite(x):
    return isinstance(x,(int,float)) and math.isfinite(x)

manifest=[]
with open(MANIFEST,newline='') as f:
    manifest=list(csv.DictReader(f))
rows=[]; counts=Counter(); warn=[]
for m in manifest:
    tid=int(m['task_id']); out=ROOT+'/'+m['output_dir']; status=text(out+'/status.txt')
    exit_code=''
    for line in status.splitlines():
        if line.startswith('exit_code='): exit_code=line.split('=',1)[1]
    err=text(out+'/stderr.log')
    r={'task_id':tid,'L':int(m['L']),'Vf':float(m['Vf']),'rate':float(m['rate']),'seed':int(m['seed']),
       'burn':int(m['burn']),'measurements_expected':int(m['measurements']),'output_dir':m['output_dir'],
       'exit_code':exit_code,'classification':'other_failure'}
    result_path=out+'/result.json'
    if os.path.getsize(result_path) if os.path.exists(result_path) else False:
        try:
            x=jload(result_path); r['classification']='success'
            for k in ['measurements_completed','S_pi','S_pi_err','average_sign','acceptance','max_sign_imag','max_observable_imag','sign_corrections','sign_recomputes','negative_signs','runtime_seconds','n_bins_used','min_update_denominator']:
                r[k]=x.get(k,'')
            r['bin_records_path']=out+'/bin_records.csv'
            flags=[]
            if x.get('measurements_completed') != int(m['measurements']): flags.append('MEASUREMENT_INCOMPLETE')
            if not finite(x.get('S_pi')) or not finite(x.get('S_pi_err')): flags.append('SPI_NONFINITE')
            if not finite(x.get('average_sign')): flags.append('SIGN_NONFINITE')
            if not finite(x.get('acceptance')) or not (0 <= x.get('acceptance') <= 1): flags.append('ACCEPTANCE_INVALID')
            if finite(x.get('max_sign_imag')) and x['max_sign_imag']>1e-6: flags.append('SIGN_IMAG_WARN')
            if finite(x.get('max_observable_imag')) and x['max_observable_imag']>1e-6: flags.append('OBS_IMAG_WARN')
            if isinstance(x.get('sign_corrections'),(int,float)) and x['sign_corrections']>0: flags.append('SIGN_CORRECTION_WARN')
            if not os.path.exists(r['bin_records_path']): flags.append('BIN_RECORDS_MISSING')
            r['qc_flags']=';'.join(flags)
            if flags: warn.append(r.copy())
        except Exception as e:
            r['classification']='other_failure'; r['parse_error']=str(e)
    elif exit_code=='2' and 'must be an integer' in err:
        r['classification']='expected_incompatible_grid'
    r['stderr_tail']=' | '.join(err.strip().splitlines()[-2:])
    counts[r['classification']]+=1; rows.append(r)
exit_counts=Counter(r['exit_code'] or 'missing' for r in rows)

fields=['task_id','L','Vf','rate','seed','burn','measurements_expected','output_dir','classification','exit_code','measurements_completed','S_pi','S_pi_err','average_sign','acceptance','max_sign_imag','max_observable_imag','sign_corrections','sign_recomputes','negative_signs','runtime_seconds','n_bins_used','min_update_denominator','bin_records_path','qc_flags','stderr_tail','parse_error']
with open(OUT+'/rigorous_153416_task_qc.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

reuse=[]
for r in rows:
    if r['L'] in (26,34) and r['rate'] in (0.5,1.0) and r['Vf'] in (3.8,3.9,4.0,4.1,4.2): reuse.append(r)
with open(OUT+'/rigorous_153416_reusable_240_qc.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(reuse)

def pooled_spi(rs):
    bins=[]
    for r in rs:
        if r['classification']!='success': continue
        try:
            with open(r['bin_records_path'],newline='') as f:
                for b in csv.DictReader(f): bins.append((float(b['sign_sum']),float(b['signed_S_pi_numerator'])))
        except Exception: pass
    if len(bins)<2: return None
    den=sum(x[0] for x in bins); num=sum(x[1] for x in bins)
    if not finite(den) or abs(den)<1e-12: return None
    val=num/den; leaves=[]
    for d,n in bins:
        if abs(den-d)>1e-12: leaves.append((num-n)/(den-d))
    if len(leaves)<2:return None
    avg=sum(leaves)/len(leaves); err=math.sqrt((len(leaves)-1)*sum((z-avg)**2 for z in leaves)/len(leaves))
    return val,err,len(bins)

pre=[]
for L in (26,34):
  for vf in (3.8,4.0,4.2):
    byrate={q:[r for r in rows if r['L']==L and abs(r['Vf']-vf)<1e-9 and r['rate']==q] for q in (0.5,1.0,2.0)}
    pools={q:pooled_spi(v) for q,v in byrate.items()}
    item={'L':L,'Vf':vf,**{f'n_R{q}':sum(r['classification']=='success' for r in byrate[q]) for q in byrate}}
    if all(pools[q] is not None and item[f'n_R{q}']==12 for q in pools):
        s05,e05,_=pools[0.5]; s1,e1,_=pools[1.0]; s2,e2,_=pools[2.0]
        q05=s1/s05; q1=s2/s1
        q05e=abs(q05)*math.hypot(e1/s1,e05/s05); q1e=abs(q1)*math.hypot(e2/s2,e1/s1)
        item.update({'S_R05':s05,'S_R05_err':e05,'S_R1':s1,'S_R1_err':e1,'S_R2':s2,'S_R2_err':e2,'Q05':q05,'Q05_err':q05e,'Q1':q1,'Q1_err':q1e,'D':q05-q1,'D_err':math.hypot(q05e,q1e),'status':'complete'})
    else: item['status']='incomplete'
    pre.append(item)
pfields=sorted({k for r in pre for k in r})
with open(OUT+'/rigorous_153416_preliminary_D.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=pfields);w.writeheader();w.writerows(pre)

reuse_success=sum(r['classification']=='success' for r in reuse)
reuse_warn=sum(bool(r.get('qc_flags')) for r in reuse if r['classification']=='success')
with open(OUT+'/rigorous_153416_summary.md','w') as f:
    f.write('# Original rigorous array 153416: collection and QC\n\n')
    f.write(f'Planned tasks: {len(rows)}. Success: {counts["success"]}. Expected incompatible-grid exclusion (exit_code=2): {counts["expected_incompatible_grid"]}. Other failures: {counts["other_failure"]}. Exit-code counts: {dict(sorted(exit_counts.items()))}.\n\n')
    f.write('All observed other failures have exit_code=3; representative `result.json.tmp` records identify `failure_reason=numerical_diagnostic_exceeds_tolerance`, rather than a scheduler failure.\n\n')
    f.write(f'## Reusable compatible subset\n\nExpected 240 tasks (L=26,34; R=0.5,1.0; Vf=3.8,3.9,4.0,4.1,4.2; 12 seeds): {reuse_success}/240 complete. QC-flagged successful reuse tasks: {reuse_warn}. Details: `rigorous_153416_reusable_240_qc.csv`.\n\n')
    f.write('## QC\n\nWarnings use flags for incomplete measurements, non-finite S_pi/sign, invalid acceptance, imaginary-part thresholds >1e-6, sign corrections, and missing bin records. Expected incompatible-grid tasks are excluded, not numerical failures.\n\n')
    f.write('## Preliminary D(V)\n\n`rigorous_153416_preliminary_D.csv` pools each rate separately from sign-reweighted bin numerators/denominators. Rates are treated as independent ensembles; Q and D errors use independent propagation. Rows are marked incomplete unless all 12 seeds exist at each of R=0.5,1,2. No Vc fit is performed here.\n')
print(dict(counts), 'reuse_success',reuse_success,'reuse_warn',reuse_warn)
