#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function
import csv, glob, json, math, os, re, subprocess
from collections import defaultdict

BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAMPAIGNS=('tiny_enumeration','tiny_ed','gaussian_exact','gaussian_qmc','hs_crosscheck','delta_symmetry','delta_ed','stability','known_physics')
OBS=('S_pi','S_pi_dq','R_cdw')

def loadj(path):
    text=open(path).read();text=re.sub(r'(?<=:)inf(?=[,}])','1e999',text);text=re.sub(r'(?<=:)-inf(?=[,}])','-1e999',text);text=re.sub(r'(?<=:)nan(?=[,}])','null',text)
    return json.loads(text)
def mean(x): return sum(x)/len(x) if x else float('nan')
def sem(x):
    if len(x)<2:return float('nan')
    m=mean(x);return math.sqrt(sum((v-m)**2 for v in x)/(len(x)*(len(x)-1)))
def normalized_value(name,value):
    if name=='boundary' and str(value) in ('PBC','OBC'):return '0' if str(value)=='PBC' else '1'
    try:
        number=float(value)
        if number.is_integer():return str(int(number))
        return '%.12g'%number
    except:return str(value)
def key(r,names): return tuple(normalized_value(x,r[x]) for x in names)
def write(path,rows,fields):
    with open(path,'w') as h:w=csv.DictWriter(h,fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

def main():
    collected=os.path.join(BASE,'collected');plots=os.path.join(BASE,'plots')
    allrows=[];progress=[];bycamp={}
    for campaign in CAMPAIGNS:
        manifest=list(csv.DictReader(open(os.path.join(BASE,'manifests',campaign+'.csv'))));done=[]
        for row in manifest:
            path=os.path.join(BASE,row['output_dir'],'result.json')
            if not os.path.isfile(path) or not os.path.getsize(path):continue
            try:data=loadj(path)
            except Exception:continue
            merged=dict(row);merged.update(data);merged['campaign']=campaign;merged['result_path']=path;done.append(merged);allrows.append(merged)
        bycamp[campaign]=done;progress.append({'campaign':campaign,'completed':len(done),'planned':len(manifest),'missing':len(manifest)-len(done)})
    write(os.path.join(collected,'progress.csv'),progress,['campaign','completed','planned','missing'])
    if allrows:write(os.path.join(collected,'all_results.csv'),allrows,sorted(set(k for r in allrows for k in r)))

    comparison=[]
    for left_name,right_name,label in [('tiny_enumeration','tiny_ed','tiny_enum_vs_ed'),('gaussian_qmc','gaussian_exact','V0_qmc_vs_gaussian')]:
        names=('L','V','theta','beta_trial','dt','delta','mu','boundary','hs_scheme')
        right={key(r,names):r for r in bycamp[right_name]}
        for l in bycamp[left_name]:
            r=right.get(key(l,names));
            if not r:continue
            for ob in OBS:comparison.append({'comparison':label,'L':l['L'],'V':l['V'],'theta':l['theta'],'delta':l['delta'],'mu':l['mu'],'boundary':l['boundary'],'observable':ob,'left':l[ob],'right':r[ob],'difference':float(l[ob])-float(r[ob])})
    write(os.path.join(collected,'exact_comparisons.csv'),comparison,['comparison','L','V','theta','delta','mu','boundary','observable','left','right','difference'])

    grouped=[]
    group_fields=('campaign','L','V','theta','delta','mu','boundary','hs_scheme','stabilization_interval','guard')
    groups=defaultdict(list)
    for r in allrows:
        if r['campaign'] in ('hs_crosscheck','delta_symmetry','stability','known_physics','gaussian_qmc'):groups[key(r,group_fields)].append(r)
    for k,items in sorted(groups.items()):
        out=dict(zip(group_fields,k));out['n']=len(items)
        for ob in OBS+('average_sign',):
            vals=[float(x[ob]) for x in items if ob in x];out[ob]=mean(vals);out[ob+'_seed_sem']=sem(vals)
        for field in ('sign_corrections','diagnostic_sign_mismatch_count'):
            out[field]=sum(float(x.get(field,0) or 0) for x in items)
        out['max_sign_imag']=max(float(x.get('max_sign_imag',0) or 0) for x in items)
        out['max_fast_full_green_rel']=max(float(x.get('diagnostic_relative_frobenius_max',0) or 0) for x in items)
        grouped.append(out)
    if grouped:write(os.path.join(collected,'grouped_results.csv'),grouped,list(group_fields)+['n']+sum(([x,x+'_seed_sem'] for x in OBS+('average_sign',)),[])+['sign_corrections','diagnostic_sign_mismatch_count','max_sign_imag','max_fast_full_green_rel'])

    paired=[]
    def add_pairs(campaign,base_names,variant_name,left_value,right_value,label):
        index=defaultdict(dict)
        for row in bycamp[campaign]:index[key(row,base_names)][normalized_value(variant_name,row[variant_name])]=row
        for base,variants in index.items():
            if normalized_value(variant_name,left_value) not in variants or normalized_value(variant_name,right_value) not in variants:continue
            left=variants[normalized_value(variant_name,left_value)];right=variants[normalized_value(variant_name,right_value)]
            for ob in OBS+('average_sign',):paired.append({'comparison':label,'group':'|'.join(base[:-1]),'seed':base[-1],'observable':ob,'left':left[ob],'right':right[ob],'difference':float(left[ob])-float(right[ob])})
    add_pairs('hs_crosscheck',('L','V','theta','boundary','seed'),'hs_scheme',0,1,'hs0_minus_hs1')
    add_pairs('delta_symmetry',('L','V','boundary','seed'),'delta',-.6,.6,'delta_minus_minus_plus')
    add_pairs('stability',('L','V','theta','boundary','seed'),'stabilization_interval',5,10,'stb5_minus_stb10')
    add_pairs('stability',('L','V','theta','boundary','seed'),'stabilization_interval',20,10,'stb20_minus_stb10')
    add_pairs('stability',('L','V','theta','boundary','seed'),'guard',0,1,'guard_off_minus_on')
    write(os.path.join(collected,'paired_differences.csv'),paired,['comparison','group','seed','observable','left','right','difference'])
    paired_summary=[]
    pgroup=defaultdict(list)
    for row in paired:pgroup[(row['comparison'],row['group'],row['observable'])].append(float(row['difference']))
    for (comparison_name,group,ob),vals in sorted(pgroup.items()):
        err=sem(vals);avg=mean(vals);z=avg/err if err and not math.isnan(err) else float('nan')
        paired_summary.append({'comparison':comparison_name,'group':group,'observable':ob,'n_pairs':len(vals),'mean_difference':avg,'paired_sem':err,'z_score':z})
    write(os.path.join(collected,'paired_summary.csv'),paired_summary,['comparison','group','observable','n_pairs','mean_difference','paired_sem','z_score'])

    known=[r for r in grouped if r['campaign']=='known_physics']
    with open(os.path.join(plots,'known_physics_R_vs_V.dat'),'w') as h:
        h.write('L V R err theta\n')
        for r in known:h.write('%s %s %s %s %s\n'%(r['L'],r['V'],r['R_cdw'],r['R_cdw_seed_sem'],r['theta']))
    with open(os.path.join(plots,'known_physics_collapse_nu1.dat'),'w') as h:
        h.write('L x R err\n')
        for r in known:
            if abs(float(r['theta'])-float(r['L']))<1e-9:h.write('%s %.12g %s %s\n'%(r['L'],(float(r['V'])-4)*float(r['L']),r['R_cdw'],r['R_cdw_seed_sem']))
    gp=os.path.join(plots,'known_physics.gp')
    open(gp,'w').write("set terminal pngcairo size 1000,650\nset key outside\nset xlabel 'V'\nset ylabel 'R_cdw'\nset output '%s'\nplot for [L in '10 18 26 34'] '%s' using (strcol(1) eq L?$2:1/0):3:4 with yerrorlines title sprintf('L=%%s',L)\nset xlabel '(V-4)L (nu=1)'\nset output '%s'\nplot for [L in '10 18 26 34'] '%s' using (strcol(1) eq L?$2:1/0):3:4 with yerrorlines title sprintf('L=%%s',L)\n"%(os.path.join(plots,'R_cdw_vs_V.png'),os.path.join(plots,'known_physics_R_vs_V.dat'),os.path.join(plots,'R_cdw_collapse_nu1.png'),os.path.join(plots,'known_physics_collapse_nu1.dat')))
    try:subprocess.call(['gnuplot',gp])
    except OSError:pass

    lines=['# static projector regression/stress summary\n']
    for p in progress:lines.append('* %s: %s/%s completed%s.\n'%(p['campaign'],p['completed'],p['planned'],' — completed' if p['missing']==0 else ' — running/needs follow-up'))
    tiny=[abs(float(x['difference'])) for x in comparison if x['comparison']=='tiny_enum_vs_ed'];gauss=[abs(float(x['difference'])) for x in comparison if x['comparison']=='V0_qmc_vs_gaussian']
    lines+=['\n## tiny exact enumeration\n',('* PASS: max observable difference %.3e.\n'%max(tiny) if tiny and progress[0]['missing']==0 and progress[1]['missing']==0 else '* needs follow-up: incomplete.\n')]
    lines+=['\n## V=0 Gaussian exact\n',('* PASS: max observable difference %.3e.\n'%max(gauss) if gauss and progress[2]['missing']==0 and progress[3]['missing']==0 else '* needs follow-up: incomplete.\n')]
    comparison_for={'hs_crosscheck':'hs0_minus_hs1','delta_symmetry':'delta_minus_minus_plus','stability':None,'known_physics':None}
    for campaign,title in [('hs_crosscheck','HS scheme cross-check'),('delta_symmetry','delta-sign symmetry'),('stability','stability regression'),('known_physics','large-L known-physics sanity')]:
        p=next(x for x in progress if x['campaign']==campaign);subset=[x for x in grouped if x['campaign']==campaign]
        bad=[x for x in subset if float(x.get('diagnostic_sign_mismatch_count',0))>0 or float(x.get('max_fast_full_green_rel',0))>1e-6]
        relevant=[x for x in paired_summary if (comparison_for[campaign] is None or x['comparison']==comparison_for[campaign]) and x['observable'] in OBS]
        pair_bad=[x for x in relevant if x['n_pairs']>=4 and abs(float(x['z_score']))>3 and abs(float(x['mean_difference']))>1e-5]
        state='PASS' if p['missing']==0 and not bad and not pair_bad else ('suspicious' if bad or pair_bad else 'needs follow-up')
        lines+=['\n## %s\n'%title,'* %s; completed %d/%d; diagnostic anomalies %d; paired observable anomalies %d.\n'%(state,p['completed'],p['planned'],len(bad),len(pair_bad))]
    open(os.path.join(BASE,'summary.md'),'w').write(''.join(lines))
    print('completed %d/%d'%(sum(x['completed'] for x in progress),sum(x['planned'] for x in progress)))
if __name__=='__main__':main()
