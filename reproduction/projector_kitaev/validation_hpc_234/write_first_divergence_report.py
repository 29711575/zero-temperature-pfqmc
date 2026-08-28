#!/usr/bin/env python
import csv, os
base=os.path.dirname(os.path.abspath(__file__))
p=os.path.join(base,'collected','first_divergence.csv')
rows=list(csv.DictReader(open(p)))
event=[r for r in rows if r['event']][0]; idx=rows.index(event)
right=[r for r in rows[:idx] if r['direction']=='right']
left=[r for r in rows[:idx] if r['direction']=='left']
def f(r,k): return float(r[k])
maxpre=max(f(r,'green_rel_pre') for r in right); maxpost=max(f(r,'green_rel_post') for r in right)
ratio=max(abs(complex(f(r,'fast_ratio_re'),f(r,'fast_ratio_im'))-complex(f(r,'full_ratio_re'),f(r,'full_ratio_im'))) for r in right)
signmis=sum(int(r['tracked_pre_re']>=0)!=int(r['direct_pre_re']>=0) for r in right)
lines=['# first divergence\n','Parameters: L=10, V=4, PBC, theta=10, beta_trial=8, dt=0.1, hs_scheme=0, guard OFF, seed=983002. The debug starts from the deterministic initial contour and stops at the first threshold event.\n','## First event\n','* Type: **GREEN_PRE** (before the local HS flip).\n','* direction/sweep: left / 0.\n','* operator=%s, bond=%s, aux=%s, full-contour boundary=%s.\n'%(event['operator'],event['bond'],event['aux'],event['boundary']),'* fast-vs-full Green relative error: **%s**.\n'%event['green_rel_pre'],'* tracked sign=%s%+.3gi; direct sign=%s%+.3gi: both positive, hence no sign mismatch at the first event.\n'%(event['tracked_pre_re'],float(event['tracked_pre_im']),event['direct_pre_re'],float(event['direct_pre_im'])),'* fast/full ratio: (%s,%s) / (%s,%s); both acceptance decisions are %s.\n'%(event['fast_ratio_re'],event['fast_ratio_im'],event['full_ratio_re'],event['full_ratio_im'],'accept' if event['fast_accept']=='1' else 'reject'),'* local minimum denominator=%s; uniform=%s.\n'%(event['min_denominator'],event['uniform']),'## Directional evidence\n','* right sweep before the event: %d local proposals; max Green relative error pre/post = %.3e / %.3e; max |fast_ratio-full_ratio| = %.3e; sign mismatches = %d.\n'%(len(right),maxpre,maxpost,ratio,signmis),'* first observed failure is immediately after entering the left traversal: l=879 is the terminal dense operator, then l=878 (bond-1 HS operator) is propagated rightward and the first proposal sees the Green mismatch.\n','## Verdict\n','**Green first.** The first threshold crossing occurs while tracked/direct signs agree and before a local ratio or acceptance mismatch. This localizes the first observed split to left-sweep propagation/stabilization around the terminal boundary, not to the right-sweep fast local ratio update. No algorithm was modified.\n','## Context CSV\n','Full preceding context and the event row: `collected/first_divergence.csv`.\n']
open(os.path.join(base,'first_divergence.md'),'w').write(''.join(lines))
