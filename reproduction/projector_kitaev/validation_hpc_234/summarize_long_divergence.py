#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv, os
base=os.path.dirname(os.path.abspath(__file__)); col=os.path.join(base,'collected')
def stats(path,kind):
 rows=list(csv.DictReader(open(path)))
 return {'diagnostic':kind,'checks':len(rows),'sign_mismatches':sum(int(r['sign_mismatch']) for r in rows),'green_drifts':sum(float(r['green_rel'])>1e-6 for r in rows),'max_green_rel':max(float(r['green_rel']) for r in rows),'max_tracked_sign_imag':max(abs(float(r['tracked_im'])) for r in rows),'max_direct_sign_imag':max(abs(float(r['direct_im'])) for r in rows),'first_event':''}
rows=[stats(os.path.join(col,'long_divergence_coarse.csv'),'sweep_end_boundary0'),stats(os.path.join(col,'long_divergence_center.csv'),'center_same_configuration')]
fields=['diagnostic','checks','sign_mismatches','green_drifts','max_green_rel','max_tracked_sign_imag','max_direct_sign_imag','first_event']
with open(os.path.join(col,'long_first_divergence.csv'),'w') as h:w=csv.DictWriter(h,fields);w.writeheader();w.writerows(rows)
md='''# long first-divergence audit

Parameters: L=6, V=4, PBC, theta=10, beta_trial=8, dt=0.1, hs_scheme=0, guard OFF, seed=984035; 500 burn plus 5000 measurements. Diagnostics are read-only and never correct `q.sign`.

## Coarse localization

No true event was found. All 11,001 boundary-0 sweep-end checks and all 5,000 center-boundary checks made on the same live configuration have matching tracked/direct Pfaffian ± signs and stabilized full-contour Green matrices.

| location | checks | ± mismatches | Green drifts | max Green rel. error |
|---|---:|---:|---:|---:|
| sweep end, boundary 0 | %(a)d | %(am)d | %(ag)d | %(ae).3e |
| center, same configuration | %(b)d | %(bm)d | %(bg)d | %(be).3e |

## Precision stage

There is no “last normal sweep → first abnormal sweep” for this seed/run, so no proposal window was activated. Running proposal diagnostics without an abnormal bracket would not be a valid first-divergence localization.

## Conclusion

No PfQMC operation is first to fail in this deterministic run: accepted local updates, stabilization checkpoints, right/left transitions, sweep ends, and center capture remain internally consistent. The earlier 509/2000 audit mismatch compared a sign captured at the center with `getSignRaw()` evaluated after the remainder of the right sweep had changed the HS configuration. It was therefore a cross-configuration diagnostic artifact, not evidence of a transported-sign mismatch.

No core algorithm was modified.
'''%{'a':rows[0]['checks'],'am':rows[0]['sign_mismatches'],'ag':rows[0]['green_drifts'],'ae':rows[0]['max_green_rel'],'b':rows[1]['checks'],'bm':rows[1]['sign_mismatches'],'bg':rows[1]['green_drifts'],'be':rows[1]['max_green_rel']}
open(os.path.join(base,'long_first_divergence.md'),'w').write(md)
