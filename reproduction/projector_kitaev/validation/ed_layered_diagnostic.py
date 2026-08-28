#!/usr/bin/env python
import json, os
import numpy as np
from scipy.linalg import expm
L=6; V=2.; N=1<<L; I=np.eye(N,dtype=complex)
def ann(site):
 a=np.zeros((N,N),complex)
 for s in range(N):
  if (s>>site)&1: a[s^(1<<site),s]=(-1)**(bin(s&((1<<site)-1)).count('1'))
 return a
def info(x):
 p=np.unravel_index(np.abs(x).argmax(),x.shape)
 return {'max_abs':float(abs(x[p])),'index':[int(p[0]),int(p[1])],'real':float(x[p].real),'imag':float(x[p].imag)}
c=[ann(i) for i in range(L)]; cd=[x.conj().T for x in c]
g=[c[i]+cd[i] for i in range(L)]+[-1j*(c[i]-cd[i]) for i in range(L)]
z=[cd[i].dot(c[i])-.5*I for i in range(L)]
A=np.zeros((2*L,2*L),complex)
for i in range(L-1):
 for k in range(2):
  x,y=k*L+i,k*L+i+1; q=1j if i%2==0 else -1j; A[x,y]=q; A[y,x]=-q
 x,y=i,i+1; A[x,y]+=1j; A[y,x]-=1j
 x,y=L+i,L+i+1; A[x,y]-=1j; A[y,x]+=1j
A[2*L-1,L]=-2j; A[L,2*L-1]=2j
errs=[]
for a in range(2*L):
 for b in range(2*L): errs.append((info(g[a].dot(g[b])+g[b].dot(g[a])-(2 if a==b else 0)*I),a,b))
worst=max(errs,key=lambda x:x[0]['max_abs']); alg=dict(worst[0],gamma_indices=[worst[1],worst[2]])
Hu=sum((.5*A[a,b]*g[a].dot(g[b]) for a in range(2*L) for b in range(a+1,2*L)),np.zeros((N,N),complex))
Ha=sum((.25*A[a,b]*g[a].dot(g[b]) for a in range(2*L) for b in range(2*L)),np.zeros((N,N),complex))
eig=np.linalg.eigvalsh(Hu)
kin={'majorana_convention':'{gamma_a,gamma_b}=2 delta_ab I','algebra':alg,'A_hermiticity':info(A-A.conj().T),'A_eigenvalues':[float(x) for x in np.linalg.eigvalsh(A)],'bilinear_all_vs_upper':info(Ha-Hu),'hermiticity':info(Hu-Hu.conj().T),'ground_energy':float(eig[0]),'first_eigenvalues':[float(x) for x in eig[:8]]}
bonds={}
for i,j in ((0,1),(5,0)):
 D=V*z[i].dot(z[j]); M=-(V/4)*g[i].dot(g[L+i]).dot(g[j]).dot(g[L+j])
 bonds['%d,%d'%(i,j)]={'D_minus_M':info(D-M),'D_hermiticity':info(D-D.conj().T),'M_hermiticity':info(M-M.conj().T)}
out=os.environ.get('ED_DIAG_OUT','.')
if not os.path.isdir(out): os.makedirs(out)
with open(os.path.join(out,'ed_majorana_algebra.json'),'w') as f: json.dump({'L':L,'max_anticommutator_error':alg,'all_pair_max_abs':max(x[0]['max_abs'] for x in errs)},f,indent=2)
with open(os.path.join(out,'ed_kinetic_only.json'),'w') as f: json.dump(kin,f,indent=2)
with open(os.path.join(out,'ed_single_bond_check.json'),'w') as f: json.dump({'L':L,'V':V,'bonds':bonds},f,indent=2)
ok=alg['max_abs']<1e-12 and kin['bilinear_all_vs_upper']['max_abs']<1e-12 and kin['hermiticity']['max_abs']<1e-12 and abs(kin['ground_energy']+6)<1e-10 and all(x['D_minus_M']['max_abs']<1e-12 for x in bonds.values())
if ok:
 HM=Hu+sum((V*z[i].dot(z[(i+1)%L]) for i in range(L)),np.zeros((N,N),complex))
 HC=sum((V*z[i].dot(z[(i+1)%L]) for i in range(L)),np.zeros((N,N),complex))
 for i in range(L):
  j=(i+1)%L
  HC += -cd[i].dot(c[j])-cd[j].dot(c[i])+cd[j].dot(cd[i])+c[i].dot(c[j])
 em=np.linalg.eigvalsh(HM); ec=np.linalg.eigvalsh(HC); sd=em-ec; first=int(np.where(np.abs(sd)>1e-10)[0][0]) if np.any(np.abs(sd)>1e-10) else None
 def obs(H,e):
  w,v=np.linalg.eigh(H); rho=np.outer(v[:,0],v[:,0].conj()); q=[]
  for x in (np.pi,np.pi+2*np.pi/L):
   Q=sum((np.exp(1j*x*i)*z[i] for i in range(L)),np.zeros((N,N),complex)); q.append(float(np.trace(rho.dot(Q.conj().T.dot(Q))).real/L**2))
  return {'E0':float(w[0]),'S_pi':q[0],'S_pi_dq':q[1],'R_cdw':1-q[1]/q[0]}
 full={'majorana':obs(HM,em),'complex_standard':obs(HC,ec),'spectrum_max_abs_diff':float(np.max(np.abs(sd))),'first_spectrum_mismatch':None if first is None else {'index':first,'majorana':float(em[first]),'complex_standard':float(ec[first]),'difference':float(sd[first])},'raw_matrix_difference':info(HM-HC),'kinetic_raw_matrix_difference':info(Hu-(HC-sum((V*z[i].dot(z[(i+1)%L]) for i in range(L)),np.zeros((N,N),complex))))}
 with open(os.path.join(out,'ed_full_hamiltonian_check.json'),'w') as f: json.dump(full,f,indent=2)
 E=sum((z[i].dot(z[i+1]) for i in range(0,L-1,2)),np.zeros((N,N),complex))
 O=sum((z[i].dot(z[i+1]) for i in range(1,L-1,2)),z[L-1].dot(z[0]))
 dt=.1; kh=expm(-.5*dt*Hu); U=kh.dot(expm(-dt*V*E)).dot(expm(-dt*V*O)).dot(kh)
 rho=np.eye(N,dtype=complex)
 for unused in range(80): rho=expm(-dt*Hu).dot(rho); rho/=np.trace(rho).real
 P=np.eye(N,dtype=complex)
 for unused in range(100): P=U.dot(P); P/=np.linalg.norm(P)
 rho=P.dot(rho).dot(P.conj().T); rho/=np.trace(rho).real
 q=[]
 for x in (np.pi,np.pi+2*np.pi/L):
  Q=sum((np.exp(1j*x*i)*z[i] for i in range(L)),np.zeros((N,N),complex)); q.append(float(np.trace(rho.dot(Q.conj().T.dot(Q))).real/L**2))
 qmc={'S_pi':(.08276162,.00036981),'S_pi_dq':(.04637077,.00010956),'R_cdw':(.43970683,.00363175)}
 project={'method':'static_projector_exact_PBC','L':L,'V':V,'delta':1.,'mu':0.,'theta':10.,'beta_trial':8.,'dt':dt,'symmetric_trotter':'exp(-dt K/2) exp(-dt V E) exp(-dt V O) exp(-dt K/2)','interaction_groups':{'even_bonds':[[0,1],[2,3],[4,5]],'odd_bonds':[[1,2],[3,4],[5,0]]},'S_pi':q[0],'S_pi_dq':q[1],'R_cdw':1-q[1]/q[0]}
 project['z_scores']={name:(project[name]-mean)/err for name,(mean,err) in qmc.items()}
 with open(os.path.join(out,'ed_projector_contour_fixed.json'),'w') as f: json.dump(project,f,indent=2)
print(json.dumps({'passed_prerequisites':ok,'algebra_error':alg['max_abs'],'kinetic_ground_energy':kin['ground_energy'],'bond_errors':{k:v['D_minus_M']['max_abs'] for k,v in bonds.items()}}))
