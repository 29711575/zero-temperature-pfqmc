#!/usr/bin/env python3
import csv, glob, json, math, os

base="results/hpc_contactcheck"
manifest={}
with open("manifest_contactcheck.csv",encoding="utf-8") as f:
    for row in csv.DictReader(f): manifest[os.path.basename(row["output_dir"])]=row

fields=["campaign","L","V","theta","beta_trial","seed","onsite_contact",
        "S_pi_offsite_mean","S_pi_offsite_err","S_pi_dq_offsite_mean","S_pi_dq_offsite_err",
        "S_pi_mean","S_pi_err","S_pi_dq_mean","S_pi_dq_err","R_cdw_mean","R_cdw_err",
        "average_sign_mean","average_sign_err","acceptance","runtime_seconds","n_bins",
        "sign_corrections","max_sign_imag","max_observable_imag","result_path"]
rows=[]
for path in sorted(glob.glob(base+"/runs/task_*/result.json")):
    with open(path,encoding="utf-8") as f: data=json.load(f)
    meta=manifest[os.path.basename(os.path.dirname(path))]
    row={k:data.get(k,meta.get(k,"")) for k in fields}
    row["campaign"]=meta["campaign"]; row["result_path"]=path
    rows.append(row)
os.makedirs(base,exist_ok=True)
with open(base+"/per_seed.csv","w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def mean_stderr(values):
    n=len(values); mean=sum(values)/n
    err=math.sqrt(sum((x-mean)**2 for x in values)/(n*(n-1))) if n>1 else 0.0
    return mean,err

groups={}
for row in rows:
    key=(int(row["L"]),float(row["V"]),float(row["theta"]),float(row["beta_trial"]))
    groups.setdefault(key,[]).append(row)
out=[]
for key,group in sorted(groups.items()):
    row={"L":key[0],"V":key[1],"theta":key[2],"beta_trial":key[3],"n_seeds":len(group),
         "onsite_contact":float(group[0]["onsite_contact"])}
    for obs in ["S_pi_offsite_mean","S_pi_dq_offsite_mean","S_pi_mean","S_pi_dq_mean",
                "R_cdw_mean","average_sign_mean"]:
        mean,err=mean_stderr([float(x[obs]) for x in group]); stem=obs[:-5] if obs.endswith("_mean") else obs
        row[stem+"_seed_mean"]=mean; row[stem+"_seed_stderr"]=err
    row["acceptance_mean"]=sum(float(x["acceptance"]) for x in group)/len(group)
    row["max_sign_imag"]=max(float(x["max_sign_imag"]) for x in group)
    row["max_observable_imag"]=max(float(x["max_observable_imag"]) for x in group)
    row["sign_corrections"]=sum(int(x["sign_corrections"]) for x in group)
    out.append(row)
with open(base+"/seed_summary.csv","w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]) if out else []); w.writeheader(); w.writerows(out)
print(f"per_seed={len(rows)} groups={len(out)}")
