#!/usr/bin/env python
from __future__ import division
import csv, glob, json, math, os, re, subprocess, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.dirname(os.path.abspath(__file__))

def load_json(path):
    raw = open(path).read()
    raw = re.sub(r'(?<=:)inf(?=[,}])', '1e999', raw)
    raw = re.sub(r'(?<=:)-inf(?=[,}])', '-1e999', raw)
    raw = re.sub(r'(?<=:)nan(?=[,}])', 'null', raw)
    return json.loads(raw)

def fnum(x):
    try: return float(x)
    except: return float('nan')

def mean(xs): return sum(xs) / len(xs) if xs else float('nan')
def sem(xs):
    if len(xs) < 2: return float('nan')
    m = mean(xs); return math.sqrt(sum((x-m)*(x-m) for x in xs)/(len(xs)-1)/len(xs))

def ratio_jk(seed_sums, key):
    # Independent seeds are the delete-one bins; preserves numerator/denominator.
    total_n = sum(x[key] for x in seed_sums); total_d = sum(x['den'] for x in seed_sums)
    val = total_n / total_d
    vals = [(total_n-x[key])/(total_d-x['den']) for x in seed_sums]
    jm = mean(vals); err = math.sqrt((len(vals)-1)/len(vals)*sum((x-jm)**2 for x in vals))
    return val, err

def read_measurements(path):
    rows=[]
    with open(path) as h:
        for r in csv.DictReader(h):
            try:
                sign=fnum(r['sign']); np=fnum(r['sign_S_pi_numerator']); nq=fnum(r['sign_S_pi_dq_numerator'])
            except KeyError: continue
            rows.append((sign,np,nq))
    return rows

def tau_int(values):
    n=len(values)
    if n < 4: return float('nan')
    m=mean(values); var=sum((x-m)*(x-m) for x in values)/n
    if var <= 0: return 0.5
    total=0.5
    # Initial-positive-sequence estimator; never crosses independent seed boundaries.
    for lag in range(1, min(n//2, 2000)):
        c=sum((values[i]-m)*(values[i+lag]-m) for i in range(n-lag))/(n-lag)/var
        if c <= 0: break
        total += c
    return total

def write_csv(path, rows, fields=None):
    if fields is None:
        fields=sorted(set(k for r in rows for k in r))
    with open(path,'w') as h:
        w=csv.DictWriter(h,fieldnames=fields,extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def png_plot(dat, out, title, xlabel, ylabel):
    gp=out+'.gp'
    open(gp,'w').write("set terminal pngcairo size 900,600\nset output '%s'\nset title '%s'\nset xlabel '%s'\nset ylabel '%s'\nset key outside\nplot '%s' using 1:2:3 with yerrorlines title columnhead(4)\n" % (out,title,xlabel,ylabel,dat))
    try: subprocess.call(['gnuplot',gp])
    except OSError: pass

def main():
    collected=os.path.join(BASE,'collected'); plots=os.path.join(BASE,'plots')
    for p in (collected,plots):
        if not os.path.isdir(p): os.makedirs(p)
    allrows=[]; failures=[]; group_rows=defaultdict(list)
    for campaign in ('trotter','autocorr','benchmark'):
        manifest=list(csv.DictReader(open(os.path.join(BASE,'manifests',campaign+'.csv'))))
        for t in manifest:
            out=os.path.join(BASE,t['output_dir']); rp=os.path.join(out,'result.json'); mp=os.path.join(out,'measurements.csv')
            try:
                r=load_json(rp)
                ms=read_measurements(mp)
                if len(ms) != int(t['measurements']): raise ValueError('measurement rows %d' % len(ms))
            except Exception as e:
                failures.append({'campaign':campaign,'task_id':t['task_id'],'output_dir':t['output_dir'],'reason':str(e)})
                continue
            row=dict(t); row.update(r); row['campaign']=campaign; row['result_path']=rp; row['measurement_path']=mp; row['measurement_rows']=len(ms)
            allrows.append(row); group_rows[campaign].append((row,ms))
    write_csv(os.path.join(collected,'all_task_results.csv'),allrows)
    write_csv(os.path.join(collected,'failed_tasks.csv'),failures,['campaign','task_id','output_dir','reason'])
    open(os.path.join(BASE,'failed_tasks.txt'),'w').write('\n'.join('%s,%s,%s' % (x['campaign'],x['task_id'],x['reason']) for x in failures)+'\n')
    warnings=[]
    for r in allrows:
        flags=[]
        if fnum(r.get('sign_corrections',0)) > 0: flags.append('SIGN_CORRECTION')
        if fnum(r.get('max_sign_imag',0)) > 1e-6: flags.append('SIGN_IMAG_WARN')
        if fnum(r.get('max_observable_imag',0)) > 1e-6: flags.append('OBSERVABLE_IMAG_WARN')
        if fnum(r.get('diagnostic_sign_mismatch_count',0)) > 0: flags.append('DIRECT_SIGN_MISMATCH')
        if flags:
            q={'campaign':r['campaign'],'task_id':r['task_id'],'seed':r['seed'],'flags':';'.join(flags),'sign_corrections':r.get('sign_corrections',0),'max_sign_imag':r.get('max_sign_imag',0),'max_observable_imag':r.get('max_observable_imag',0),'diagnostic_sign_mismatch_count':r.get('diagnostic_sign_mismatch_count',0)}
            warnings.append(q)
    write_csv(os.path.join(collected,'qc_warnings.csv'),warnings,['campaign','task_id','seed','flags','sign_corrections','max_sign_imag','max_observable_imag','diagnostic_sign_mismatch_count'])

    # Trotter: group matched seeds by dt and use leave-one-seed-out ratios.
    tgroups=defaultdict(list)
    for row,ms in group_rows['trotter']:
        tgroups[fnum(row['dt'])].append((row,ms))
    fitrows=[]; byobs=defaultdict(list); dt2_bad=[]
    for obs, key in [('S_pi','np'),('S_pi_dq','nq')]:
        pts=[]
        for dt in sorted(tgroups, reverse=True):
            sums=[]
            for row,ms in tgroups[dt]: sums.append({'den':sum(x[0] for x in ms),key:sum(x[1 if key=='np' else 2] for x in ms)})
            v,e=ratio_jk(sums,key); pts.append((dt,dt*dt,v,e))
        # weighted O0+a dt2
        sw=sum(1/(e*e) for _,_,_,e in pts); sx=sum(x/(e*e) for _,x,_,e in pts); sy=sum(y/(e*e) for _,x,y,e in pts); sxx=sum(x*x/(e*e) for _,x,_,e in pts); sxy=sum(x*y/(e*e) for _,x,y,e in pts)
        det=sw*sxx-sx*sx; o0=(sxx*sy-sx*sxy)/det; slope=(sw*sxy-sx*sy)/det; oerr=math.sqrt(sxx/det); aerr=math.sqrt(sw/det)
        chi=sum(((y-(o0+slope*x))/e)**2 for _,x,y,e in pts)
        for dt,x,y,e in pts:
            z=(y-o0)/math.sqrt(e*e+oerr*oerr)
            rr={'observable':obs,'dt':dt,'dt2':x,'value':y,'error':e,'O0':o0,'O0_err':oerr,'slope':slope,'slope_err':aerr,'chi2':chi,'dof':1,'dt01_minus_O0_sigma':z}
            fitrows.append(rr); byobs[obs].append(rr)
            if dt==0.2 and abs((y-(o0+slope*x))/e)>2: dt2_bad.append(obs)
    # R uses seed-level ratio transform under delete-one seed samples.
    pts=[]
    for dt in sorted(tgroups, reverse=True):
        ss=[]
        for row,ms in tgroups[dt]: ss.append({'den':sum(x[0] for x in ms),'np':sum(x[1] for x in ms),'nq':sum(x[2] for x in ms)})
        totalp=sum(x['np'] for x in ss); totalq=sum(x['nq'] for x in ss); val=1-totalq/totalp
        jk=[1-(totalq-x['nq'])/(totalp-x['np']) for x in ss]; jm=mean(jk); er=math.sqrt((len(jk)-1)/len(jk)*sum((x-jm)**2 for x in jk)); pts.append((dt,dt*dt,val,er))
    sw=sum(1/(e*e) for _,_,_,e in pts); sx=sum(x/(e*e) for _,x,_,e in pts); sy=sum(y/(e*e) for _,x,y,e in pts); sxx=sum(x*x/(e*e) for _,x,_,e in pts); sxy=sum(x*y/(e*e) for _,x,y,e in pts); det=sw*sxx-sx*sx; o0=(sxx*sy-sx*sxy)/det; slope=(sw*sxy-sx*sy)/det; oerr=math.sqrt(sxx/det); aerr=math.sqrt(sw/det); chi=sum(((y-(o0+slope*x))/e)**2 for _,x,y,e in pts)
    for dt,x,y,e in pts:
        fitrows.append({'observable':'R_cdw','dt':dt,'dt2':x,'value':y,'error':e,'O0':o0,'O0_err':oerr,'slope':slope,'slope_err':aerr,'chi2':chi,'dof':1,'dt01_minus_O0_sigma':(y-o0)/math.sqrt(e*e+oerr*oerr)})
        byobs['R_cdw'].append(fitrows[-1])
        if dt==0.2 and abs((y-(o0+slope*x))/e)>2: dt2_bad.append('R_cdw')
    write_csv(os.path.join(collected,'trotter_fit.csv'),fitrows,['observable','dt','dt2','value','error','O0','O0_err','slope','slope_err','chi2','dof','dt01_minus_O0_sigma'])
    for obs,rs in byobs.items():
        dat=os.path.join(plots,'trotter_%s.dat' % obs); h=open(dat,'w'); h.write('dt2 value error label\n'); [h.write('%.12g %.12g %.12g %s\n' % (r['dt2'],r['value'],r['error'],obs)) for r in sorted(rs,key=lambda x:x['dt2'])]; h.close(); png_plot(dat,os.path.join(plots,'trotter_%s.png'%obs),obs+' vs dt^2','dt^2',obs)

    # Autocorrelation: old V=2 raw data plus new V=4 campaign.
    acsources=[]
    old=sorted(glob.glob(os.path.join(ROOT,'results','trial_density_convergence_theta10','raw','task_*_L10_V2_th10_bt8_dt0.1_s*','measurements.csv')))
    acsources.append(('L10_V2_existing',old))
    acsources.append(('L10_V4_new',[row['measurement_path'] for row,ms in group_rows['autocorr']]))
    taurows=[]; blockrows=[]; recommendations=[]
    for label,paths in acsources:
        series=[]
        for p in paths:
            ms=read_measurements(p); series.append(ms)
            for name,col in [('sign',0),('sign_S_pi',1),('sign_S_pi_dq',2)]: taurows.append({'dataset':label,'seed_path':p,'observable':name,'tau_int':tau_int([x[col] for x in ms]),'measurements':len(ms)})
        for name,col in [('sign',0),('sign_S_pi',1),('sign_S_pi_dq',2)]:
            vals=[r['tau_int'] for r in taurows if r['dataset']==label and r['observable']==name]; taurows.append({'dataset':label,'seed_path':'POOLED_MEDIAN','observable':name,'tau_int':sorted(vals)[len(vals)//2],'measurements':sum(len(x) for x in series)})
        mt=max(r['tau_int'] for r in taurows if r['dataset']==label and r['seed_path']!='POOLED_MEDIAN'); recommendations.append({'dataset':label,'max_seed_tau_int':mt,'recommended_bin_size':max(25,int(math.ceil(10*mt)))})
        for bs in (25,50,100,250,500,1000):
            blocks=[]
            for s in series:
                for j in range(0,len(s)-bs+1,bs):
                    q=s[j:j+bs]; blocks.append((sum(x[0] for x in q),sum(x[1] for x in q),sum(x[2] for x in q)))
            for name,col in [('sign',0),('sign_S_pi',1),('sign_S_pi_dq',2)]:
                vals=[x[col]/bs for x in blocks]; blockrows.append({'dataset':label,'block_size':bs,'observable':name,'n_blocks':len(vals),'mean':mean(vals),'error':sem(vals)})
            # Ratio error stays explicitly as numerator/denominator ratio of block sums.
            for name,col in [('S_pi_ratio',1),('S_pi_dq_ratio',2)]:
                tn=sum(x[col] for x in blocks); td=sum(x[0] for x in blocks); jk=[(tn-x[col])/(td-x[0]) for x in blocks]; jm=mean(jk); er=math.sqrt((len(jk)-1)/len(jk)*sum((x-jm)**2 for x in jk)); blockrows.append({'dataset':label,'block_size':bs,'observable':name,'n_blocks':len(blocks),'mean':tn/td,'error':er})
    write_csv(os.path.join(collected,'autocorr_tau.csv'),taurows,['dataset','seed_path','observable','tau_int','measurements'])
    write_csv(os.path.join(collected,'blocking_errors.csv'),blockrows,['dataset','block_size','observable','n_blocks','mean','error'])
    write_csv(os.path.join(collected,'blocking_recommendations.csv'),recommendations)
    for label,_ in acsources:
        dat=os.path.join(plots,'blocking_%s.dat'%label); h=open(dat,'w'); h.write('block error unused label\n');
        for r in blockrows:
            if r['dataset']==label and r['observable'] in ('sign','sign_S_pi','sign_S_pi_dq'): h.write('%s %.12g 0 %s\n'%(r['block_size'],r['error'],r['observable']))
        h.close(); gp=dat+'.gp'; out=os.path.join(plots,'blocking_%s.png'%label); open(gp,'w').write("set terminal pngcairo size 900,600\nset output '%s'\nset title 'Blocking errors: %s'\nset xlabel 'block size'\nset ylabel 'standard error'\nset key outside\nplot '%s' using 1:2:3 with linespoints title 'sign', '' using 1:2:3 with linespoints title 'sign*S_pi', '' using 1:2:3 with linespoints title 'sign*S_pi_dq'\n"%(out,label,dat)); subprocess.call(['gnuplot',gp])

    # Benchmark pooled seed jackknife and ED z scores.
    bgroups=defaultdict(list)
    for row,ms in group_rows['benchmark']: bgroups[(int(row['boundary']),fnum(row['V']))].append((row,ms))
    brows=[]
    for (b,v),items in sorted(bgroups.items()):
        ss=[{'den':sum(x[0] for x in ms),'np':sum(x[1] for x in ms),'nq':sum(x[2] for x in ms)} for row,ms in items]
        sp,ep=ratio_jk(ss,'np'); sq,eq=ratio_jk(ss,'nq'); tp=sum(x['np'] for x in ss); tq=sum(x['nq'] for x in ss); rr=1-tq/tp; jk=[1-(tq-x['nq'])/(tp-x['np']) for x in ss]; jm=mean(jk); er=math.sqrt((len(jk)-1)/len(jk)*sum((x-jm)**2 for x in jk))
        bc='PBC' if b==0 else 'OBC'; ed=load_json(os.path.join(collected,'ed','%s_V%d.json'%(bc,int(v))))
        first=[row for row,ms in items]; signs=[fnum(x['average_sign']) for x in first]
        for name,val,err in [('S_pi',sp,ep),('S_pi_dq',sq,eq),('R_cdw',rr,er)]: brows.append({'boundary':bc,'V':v,'observable':name,'QMC':val,'QMC_err':err,'ED':ed[name],'QMC_minus_ED':val-ed[name],'z_score':((val-ed[name])/err if err > 1e-12 else float('nan')),'average_sign':mean(signs),'average_sign_seed_SEM':sem(signs),'acceptance':mean([fnum(x['acceptance']) for x in first]),'sign_corrections':sum(fnum(x.get('sign_corrections',0)) for x in first),'max_sign_imag':max(fnum(x.get('max_sign_imag',0)) for x in first),'max_observable_imag':max(fnum(x.get('max_observable_imag',0)) for x in first),'n_seeds':len(items)})
    write_csv(os.path.join(collected,'benchmark_ed_qmc.csv'),brows,['boundary','V','observable','QMC','QMC_err','ED','QMC_minus_ED','z_score','average_sign','average_sign_seed_SEM','acceptance','sign_corrections','max_sign_imag','max_observable_imag','n_seeds'])
    # Optional dt=.15 manifest/script only when dt=.20 has a >2 sigma fit residual.
    if dt2_bad:
        src=list(csv.DictReader(open(os.path.join(BASE,'manifests','trotter.csv')))); opt=[]
        for t in src[:12]:
            q=dict(t); q['task_id']=str(len(opt)); q['dt']='0.15'; q['output_dir']='raw/trotter_dt015/task_%03d_s%s'%(len(opt),q['seed']); opt.append(q)
        write_csv(os.path.join(BASE,'manifests','trotter_dt015_optional.csv'),opt,list(src[0].keys()))
    summary=[]
    summary.append('# validation_hpc_234 summary\n')
    summary.append('Completed QMC: %d/138; failed: %d. Same-contour ED: 8/8.\n' % (len(allrows),len(failures)))
    summary.append('## Trotter O(dt^2)\n')
    for o in ('S_pi','S_pi_dq','R_cdw'):
        r=[x for x in fitrows if x['observable']==o][0]; d=[x for x in fitrows if x['observable']==o and abs(x['dt']-.1)<1e-10][0]; summary.append('* %s: O0=%.8g +/- %.3g, a=%.8g +/- %.3g, chi2/dof=%.3g/1; dt=0.1 minus O0 = %.3g sigma.\n'%(o,r['O0'],r['O0_err'],r['slope'],r['slope_err'],r['chi2'],d['dt01_minus_O0_sigma']))
    summary.append('dt=0.2 >2sigma fit-residual observables: %s.\n'%(', '.join(sorted(set(dt2_bad))) if dt2_bad else 'none'))
    summary.append('## Autocorrelation / blocking\n')
    for r in recommendations: summary.append('* %s: max seed tau_int=%.3g; recommended production bin >= %d.\n'%(r['dataset'],r['max_seed_tau_int'],r['recommended_bin_size']))
    summary.append('## L=6 ED-QMC benchmark\n')
    for r in brows: summary.append('* %s V=%g %s: QMC %.8g +/- %.3g; ED %.8g; z=%.3g.\n'%(r['boundary'],r['V'],r['observable'],r['QMC'],r['QMC_err'],r['ED'],r['z_score']))
    summary.append('## Failures / anomalies\n')
    summary.append('* Failed tasks: none.\n' if not failures else ''.join('* failed %s %s: %s\n'%(x['campaign'],x['task_id'],x['reason']) for x in failures))
    if warnings:
        summary.append('* Diagnostic WARN tasks (%d; full detail: collected/qc_warnings.csv): %s.\n' % (len(warnings), ', '.join('%s:%s[%s]'%(x['campaign'],x['task_id'],x['flags']) for x in warnings)))
    else: summary.append('* Diagnostic WARN tasks: none.\n')
    open(os.path.join(BASE,'summary.md'),'w').write(''.join(summary))
    print('completed=%d failures=%d' % (len(allrows),len(failures)))

if __name__=='__main__':
    main()
    # The legacy block above writes diagonal-weight Trotter fits.  Always
    # replace them with the authoritative matched-(seed,bin) correlated
    # jackknife products before returning from a full collector rerun.
    subprocess.check_call([sys.executable, os.path.join(BASE, 'correlated_trotter_jackknife.py')])
