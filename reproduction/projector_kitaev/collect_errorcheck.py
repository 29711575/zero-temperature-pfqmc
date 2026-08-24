#!/usr/bin/env python3
import csv, glob, json, math, os

base="results/hpc_errorcheck"
manifest={}
with open("manifest_errorcheck.csv",encoding="utf-8") as f:
    for row in csv.DictReader(f): manifest[os.path.basename(row["output_dir"])]=row

fields=["campaign","L","V","theta","beta_trial","seed","S_pi_mean","S_pi_err",
        "S_pi_dq_mean","S_pi_dq_err","R_cdw_mean","R_cdw_err",
        "average_sign_mean","average_sign_err","acceptance","runtime_seconds",
        "n_bins","sign_corrections","max_sign_imag","max_observable_imag","result_path"]
rows=[]
for path in sorted(glob.glob(base+"/runs/task_*/result.json")):
    with open(path,encoding="utf-8") as f: data=json.load(f)
    meta=manifest[os.path.basename(os.path.dirname(path))]
    row={k:data.get(k,meta.get(k,"")) for k in fields}; row["campaign"]=meta["campaign"]; row["result_path"]=path
    rows.append(row)
with open(base+"/per_seed.csv","w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def mean_stderr(values):
    n=len(values); mean=sum(values)/n
    err=math.sqrt(sum((x-mean)**2 for x in values)/(n*(n-1))) if n>1 else 0.0
    return mean,err
def weighted(values,errors):
    pairs=[(x,e) for x,e in zip(values,errors) if e>0 and math.isfinite(e)]
    if not pairs: return mean_stderr(values)
    weights=[1/(e*e) for _,e in pairs]; sw=sum(weights)
    return sum(w*x for w,(x,_) in zip(weights,pairs))/sw,math.sqrt(1/sw)

groups={}
for r in rows:
    key=(r["campaign"],int(r["L"]),float(r["V"]),float(r["theta"]),float(r["beta_trial"]))
    groups.setdefault(key,[]).append(r)
out=[]
for key,group in sorted(groups.items()):
    row={"campaign":key[0],"L":key[1],"V":key[2],"theta":key[3],"beta_trial":key[4],"n_seeds":len(group)}
    for obs,errkey in [("S_pi_mean","S_pi_err"),("S_pi_dq_mean","S_pi_dq_err"),
                       ("R_cdw_mean","R_cdw_err"),("average_sign_mean","average_sign_err")]:
        vals=[float(x[obs]) for x in group]; errs=[float(x[errkey]) for x in group]
        m,e=mean_stderr(vals); wm,we=weighted(vals,errs); stem=obs.removesuffix("_mean")
        row[stem+"_seed_mean"]=m; row[stem+"_seed_stderr"]=e
        row[stem+"_ivw_mean"]=wm; row[stem+"_ivw_err"]=we
    row["acceptance_mean"]=sum(float(x["acceptance"]) for x in group)/len(group)
    row["runtime_mean"]=sum(float(x["runtime_seconds"]) for x in group)/len(group)
    out.append(row)
summary_fields=list(out[0]) if out else []
with open(base+"/seed_summary.csv","w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=summary_fields); w.writeheader(); w.writerows(out)
print(f"per_seed={len(rows)} groups={len(out)}")
