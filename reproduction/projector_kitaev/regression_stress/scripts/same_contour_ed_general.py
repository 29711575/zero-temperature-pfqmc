#!/usr/bin/env python
import argparse
import json
import numpy as np
from scipy.linalg import expm

def ann(L,i):
 n=1<<L; a=np.zeros((n,n),complex)
 for s in range(n):
  if (s>>i)&1: a[s^(1<<i),s]=(-1)**(bin(s&((1<<i)-1)).count('1'))
 return a

def main():
 p=argparse.ArgumentParser(); p.add_argument('--L',type=int,required=True); p.add_argument('--V',type=float,required=True); p.add_argument('--theta',type=float,required=True); p.add_argument('--beta-trial',type=float,required=True); p.add_argument('--dt',type=float,required=True); p.add_argument('--boundary',type=int,choices=(0,1),default=0); p.add_argument('--delta',type=float,default=1.0); p.add_argument('--mu',type=float,default=0.0); a=p.parse_args()
 L=a.L; n=1<<L; c=[ann(L,i) for i in range(L)]; cd=[x.conj().T for x in c]; g=[x+y for x,y in zip(c,cd)]+[-1j*(x-y) for x,y in zip(c,cd)]
 A=np.zeros((2*L,2*L),complex); E=np.zeros((n,n),complex); O=np.zeros((n,n),complex); z=[cd[i].dot(c[i])-0.5*np.eye(n) for i in range(L)]
 for i in range(L-1):
  j=(i+1)%L; phase=1j if i%2==0 else -1j
  for k in range(2):
   x,y=k*L+i,k*L+j; A[x,y]+=phase; A[y,x]+=-phase
   if k==0: A[x,y]+=1j*a.delta; A[y,x]+=-1j*a.delta
   else: A[x,y]+=-1j*a.delta; A[y,x]+=1j*a.delta
  term=z[i].dot(z[j])
  if i%2==0: E+=term
  else: O+=term
 if a.boundary==0 and L%2==0:
  O+=z[L-1].dot(z[0])
  a0,b0=L-1,0; a1,b1=2*L-1,L
  A[a0,b0]+=1j*(a.delta-1); A[b0,a0]+=-1j*(a.delta-1)
  A[a1,b1]+=-1j*(a.delta+1); A[b1,a1]+=+1j*(a.delta+1)
 for i in range(L):
  A[i,L+i]+=-1j*a.mu; A[L+i,i]+=1j*a.mu
 if np.max(np.abs(A+A.T)) >= 1e-12:
  raise AssertionError('Majorana kinetic matrix is not antisymmetric')
 K=np.zeros((n,n),complex)
 for i in range(2*L):
  for j in range(i+1,2*L):
   K+=0.5*A[i,j]*g[i].dot(g[j])
 kh=expm(-0.5*a.dt*K); U=kh.dot(expm(-a.dt*a.V*E)).dot(expm(-a.dt*a.V*O)).dot(kh)
 H=K+a.V*(E+O); spectrum=np.linalg.eigvalsh(H).real
 rho=np.eye(n,dtype=complex)
 remaining=a.beta_trial
 while remaining>1e-12:
  step=min(a.dt,remaining); rho=expm(-step*K).dot(rho); rho/=np.trace(rho).real; remaining-=step
 P=np.eye(n,dtype=complex)
 for _ in range(int(round(a.theta/a.dt))): P=U.dot(P); P/=np.linalg.norm(P)
 rho=P.dot(rho).dot(P.conj().T); rho/=np.trace(rho).real
 def values(state):
  vals=[]
  for q in (np.pi,np.pi+2*np.pi/L):
   Q=sum((np.exp(1j*q*i)*z[i] for i in range(L)),np.zeros((n,n),complex)); vals.append(float((np.trace(state.dot(Q.conj().T.dot(Q)))/np.trace(state)).real/L**2))
  return vals
 full=values(rho); parity=np.array([1 if bin(s).count('1')%2==0 else -1 for s in range(n)])
 even=rho*(((parity[:,None]+1)/2)*((parity[None,:]+1)/2)); odd=rho*(((1-parity[:,None])/2)*((1-parity[None,:])/2))
 ev=values(even); od=values(odd)
 boundary_name='PBC' if a.boundary==0 else 'OBC'
 interaction_bonds=list(range(L)) if a.boundary==0 else list(range(L-1))
 print(json.dumps({'method':'static_projector_same_contour_exact_general','L':L,'V':a.V,'theta':a.theta,'beta_trial':a.beta_trial,'dt':a.dt,'delta':a.delta,'mu':a.mu,'boundary':boundary_name,'S_pi':full[0],'S_pi_dq':full[1],'R_cdw':1-full[1]/full[0],'even_parity':[ev[0],ev[1],1-ev[1]/ev[0]],'odd_parity':[od[0],od[1],1-od[1]/od[0]],'ground_energy':float(spectrum[0]),'first_eigenvalues':spectrum[:6].tolist(),'kinetic_boundary_A':[complex(A[L-1,0]).real,complex(A[2*L-1,L]).real,complex(A[L-1,0]).imag,complex(A[2*L-1,L]).imag],'interaction_bonds':interaction_bonds},separators=(',',':')))
if __name__=='__main__': main()
