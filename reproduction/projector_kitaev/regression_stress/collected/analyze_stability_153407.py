#!/usr/bin/env python3
import csv, json, math, statistics
from pathlib import Path

BASE=Path('/home/sunxr/PfQMC-main/reproduction/projector_kitaev/regression_stress')
OUT=BASE/'collected'
MAN=BASE/'manifests/stability.csv'
SIGN_Z_MIN=3.0
METRICS=['S_pi','S_pi_dq','R_cdw','average_sign','acceptance','sign_corrections',
         'max_sign_imag','max_observable_imag','diagnostic_relative_frobenius_max',
         'diagnostic_S_pi_abs_diff_max','diagnostic_R_cdw_abs_diff_max',
         'diagnostic_sign_mismatch_count']

def ffloat(x):
    try:return float(x)
    except:return math.nan
def mean_se(x):
    x=[v for v in x if math.isfinite(v)]
    if not x:return math.nan,math.nan,0
    m=sum(x)/len(x)
    se=statistics.stdev(x)/math.sqrt(len(x)) if len(x)>1 else math.nan
    return m,se,len(x)
def write(name,rows,fields=None):
    p=OUT/name
    if fields is None: fields=list(rows[0]) if rows else []
    with p.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)

manifest=list(csv.DictReader(MAN.open()))
tasks=[]; byid={}
for r in manifest:
    d=BASE/r['output_dir']; status='missing'; rc=''
    sp=d/'status.txt'
    if sp.exists():
        txt=sp.read_text(); status='complete' if 'status=complete' in txt else 'failed'
        if 'exit_code=' in txt: rc=txt.split('exit_code=',1)[1].splitlines()[0]
    data={}
    if (d/'result.json').exists(): data=json.loads((d/'result.json').read_text())
    row={k:r[k] for k in r}
    row.update(status=status,exit_code=rc,result_present=int(bool(data)),measurement_rows=0)
    mp=d/'measurements.csv'
    if mp.exists(): row['measurement_rows']=max(0,sum(1 for _ in mp.open())-1)
    for k in METRICS+['average_sign_err','sign_corrections','sign_recomputes','sign_corrections',
                      'adaptive_rebuild_count','guard_trigger_fraction','min_update_denominator',
                      'diagnostic_comparisons']:
        row[k]=data.get(k,math.nan)
    ase=ffloat(row['average_sign_err']); av=ffloat(row['average_sign'])
    row['average_sign_z']=abs(av)/ase if ase>0 else (math.inf if av else math.nan)
    row['physics_valid']=int(row['average_sign_z']>=SIGN_Z_MIN)
    tasks.append(row);byid[int(r['task_id'])]=(r,data,d,row)
write('stability_153407_task_summary.csv',tasks)

# Build matched task pairs from manifest-defined groups.
pairs=[]
for seed in sorted({int(r['seed']) for r in manifest[:18]}):
    rr=[x for x in manifest if int(x['seed'])==seed and int(x['L'])==18]
    for aint,bint in [(5,10),(5,20),(10,20)]:
        a=next(x for x in rr if int(x['stabilization_interval'])==aint)
        b=next(x for x in rr if int(x['stabilization_interval'])==bint)
        pairs.append(('interval',f'{aint}_vs_{bint}',a,b))
for V,seeds in [(4,range(740007,740013)),(6,range(740013,740019))]:
    for seed in seeds:
        rr=[x for x in manifest if int(x['seed'])==seed and int(x['V'])==V]
        a=next(x for x in rr if int(x['guard'])==0); b=next(x for x in rr if int(x['guard'])==1)
        pairs.append(('guard',f'V{V}_off_vs_on',a,b))

diffs=[]
for kind,comparison,a,b in pairs:
    da=byid[int(a['task_id'])][1];db=byid[int(b['task_id'])][1]
    valid=byid[int(a['task_id'])][3]['physics_valid'] and byid[int(b['task_id'])][3]['physics_valid']
    for metric in METRICS:
        va=ffloat(da.get(metric));vb=ffloat(db.get(metric))
        diffs.append(dict(kind=kind,comparison=comparison,seed=a['seed'],metric=metric,
                          task_a=a['task_id'],task_b=b['task_id'],value_a=va,value_b=vb,
                          difference_b_minus_a=vb-va,physics_pair_valid=int(valid)))
write('stability_153407_matched_seed_differences.csv',diffs)

summ=[]
for key in sorted({(x['kind'],x['comparison'],x['metric']) for x in diffs}):
    kind,comp,metric=key; rows=[x for x in diffs if (x['kind'],x['comparison'],x['metric'])==key]
    physical=metric in ('S_pi','S_pi_dq','R_cdw')
    used=[x for x in rows if (not physical or x['physics_pair_valid'])]
    m,se,n=mean_se([x['difference_b_minus_a'] for x in used]); z=m/se if se>0 else math.nan
    if physical and not used: verdict='invalid_low_sign'
    elif math.isfinite(z) and abs(z)>=3: verdict='suspicious'
    else: verdict='pass'
    summ.append(dict(kind=kind,comparison=comp,metric=metric,n_pairs_total=len(rows),n_pairs_used=n,
                     mean_difference=m,se_difference=se,z_paired=z,verdict=verdict))
write('stability_153407_matched_pair_summary.csv',summ)

# Matched 20-bin ratios. A bin is usable only if its own sign mean is >=3 sigma from zero.
def bins(task):
    r,data,d,trow=byid[int(task['task_id'])]
    rows=list(csv.DictReader((d/'measurements.csv').open())); n=len(rows); out=[]
    for bi in range(20):
        lo=bi*n//20;hi=(bi+1)*n//20; q=rows[lo:hi]
        s=[float(x['sign']) for x in q]; sm=sum(s); mn=sm/len(s)
        se=math.sqrt(max(0,1-mn*mn)/(len(s)-1)) if len(s)>1 else math.nan
        z=abs(mn)/se if se>0 else (math.inf if mn else 0)
        np=sum(float(x['sign_S_pi_numerator']) for x in q)
        nd=sum(float(x['sign_S_pi_dq_numerator']) for x in q)
        p=np/sm if sm else math.nan; dq=nd/sm if sm else math.nan
        out.append(dict(bin=bi,sign_mean=mn,sign_z=z,valid=int(z>=SIGN_Z_MIN),S_pi=p,S_pi_dq=dq,
                        R_cdw=1-dq/p if p and math.isfinite(p) else math.nan))
    return out
bd=[]
for kind,comparison,a,b in pairs:
    ba,bb=bins(a),bins(b)
    for xa,xb in zip(ba,bb):
        for metric in ('S_pi','S_pi_dq','R_cdw'):
            valid=xa['valid'] and xb['valid']
            bd.append(dict(kind=kind,comparison=comparison,seed=a['seed'],bin=xa['bin'],metric=metric,
                           value_a=xa[metric],value_b=xb[metric],difference_b_minus_a=xb[metric]-xa[metric],
                           sign_z_a=xa['sign_z'],sign_z_b=xb['sign_z'],bin_pair_valid=int(valid)))
write('stability_153407_matched_bin_differences.csv',bd)
bs=[]
for key in sorted({(x['kind'],x['comparison'],x['metric']) for x in bd}):
    kind,comp,metric=key; rows=[x for x in bd if (x['kind'],x['comparison'],x['metric'])==key]
    used=[x for x in rows if x['bin_pair_valid']]
    m,se,n=mean_se([x['difference_b_minus_a'] for x in used]);z=m/se if se>0 else math.nan
    bs.append(dict(kind=kind,comparison=comp,metric=metric,n_bins_total=len(rows),n_bins_used=n,
                   mean_difference=m,se_difference=se,z_paired=z,
                   verdict=('invalid_low_sign' if not used else ('suspicious' if math.isfinite(z) and abs(z)>=3 else 'pass'))))
write('stability_153407_matched_bin_summary.csv',bs)
print('tasks',len(tasks),'complete',sum(x['status']=='complete' for x in tasks),'valid_physics',sum(x['physics_valid'] for x in tasks))
