#!/usr/bin/env python3
import csv, subprocess, sys
from pathlib import Path

base=Path(__file__).resolve().parent
task=int(sys.argv[1]); row=list(csv.DictReader((base/'precision_manifest.csv').open()))[task]
out=base/'results'/'precision'/f'task_{task:02d}'
out.mkdir(parents=True,exist_ok=False)
common=['18','18','6','8','0.1','1','1',row['seed'],row['half_step']]
if row['mode']=='checkpoint': cmd=[str(base/'bin'/'validation_oracle'),'checkpoint',*common,row['precision'],row['expected_hs_hash'],str(out/'precision.csv')]
else: cmd=[str(base/'bin'/'validation_oracle'),'adaptive_checkpoint',*common,row['expected_hs_hash'],str(out/'adaptive.csv')]
p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
(out/'command.txt').write_text(' '.join(cmd)+'\n');(out/'stdout.log').write_text(p.stdout);(out/'stderr.log').write_text(p.stderr)
if p.returncode: raise SystemExit(p.returncode)

