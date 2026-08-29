#!/usr/bin/env python3
import csv, json, os, subprocess, sys
from pathlib import Path

base=Path(__file__).resolve().parent
task=int(sys.argv[1])
rows=list(csv.DictReader((base/'replay_manifest.csv').open()))
row=rows[task]
result_root=Path(os.environ.get('MPZ2_RESULT_ROOT',str(base/'results')))
out=result_root/'replay'/row['label']
out.mkdir(parents=True,exist_ok=False)
cmd=[str(base/'bin'/'projector_real_z2_driver'),row['L'],row['theta'],row['beta_trial'],row['dt'],row['V'],row['delta'],row['mu'],row['boundary'],row['hs_scheme'],row['seed'],row['burn'],row['measurements'],'1',str(out/'measurements.csv'),row['diagnostic_stride'],row['sign_stride'],str(out/'legacy.csv'),row['mp_spot_stride']]
p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
(out/'command.txt').write_text(' '.join(cmd)+'\n')
(out/'stdout.log').write_text(p.stdout)
(out/'stderr.log').write_text(p.stderr)
if p.returncode: raise SystemExit(p.returncode)
result=json.loads(p.stdout.strip().splitlines()[-1])
result['expected_trajectory_hash']=int(row['expected_trajectory_hash'])
result['trajectory_gate_pass']=result['shadow_trajectory_match'] and result['trajectory_hash']==int(row['expected_trajectory_hash'])
result['expected_average_sign']=float(row['expected_average_sign'])
(out/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
if not result['trajectory_gate_pass'] or result.get('mp_correction_count')!=0: raise SystemExit(3)
