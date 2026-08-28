#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path

base=Path(__file__).resolve().parent; task=int(sys.argv[1]); out=base/'results'/'controls'/f'task_{task:02d}';out.mkdir(parents=True,exist_ok=False)
exe=str(base/'bin'/'validation_oracle')
if task==0: cmd=[exe,'control','6','6','6','8','0.1','1','1','2126143300','4','1',str(out/'control.csv')]
elif task==1: cmd=[exe,'control','6','6','6','8','0.1','1','1','2126143301','4','1',str(out/'control.csv')]
elif task==2: cmd=[exe,'control','12','12','6','8','0.1','1','1','2126203300','3','0',str(out/'control.csv')]
elif task==3: cmd=[str(base/'run_regressions.sh'),str(out)]
else: raise SystemExit('invalid task')
p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
(out/'command.txt').write_text(' '.join(cmd)+'\n');(out/'stdout.log').write_text(p.stdout);(out/'stderr.log').write_text(p.stderr)
if p.returncode: raise SystemExit(p.returncode)

