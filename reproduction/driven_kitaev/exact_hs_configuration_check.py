#!/usr/bin/env python3
"""Independent 2^L Fock-space referee for frozen driven HS contours."""
import argparse,csv,json,math
from pathlib import Path
import numpy as np
from scipy.linalg import expm

def annihilation(L,i):
 n=1<<L;a=np.zeros((n,n),complex)
 for s in range(n):
  if s>>i&1:a[s^(1<<i),s]=(-1)**(s&((1<<i)-1)).bit_count()
 return a

def kinetic(L):
 c=[annihilation(L,i) for i in range(L)];cd=[x.conj().T for x in c];g=[x+y for x,y in zip(c,cd)]+[-1j*(x-y) for x,y in zip(c,cd)];A=np.zeros((2*L,2*L),complex)
 for i in range(L-1):
  for k in range(2):
   a,b=k*L+i,k*L+i+1;z=1j if i%2==0 else -1j;A[a,b]=z;A[b,a]=-z
  A[i,i+1]+=1j;A[i+1,i]+=-1j;A[L+i,L+i+1]+=-1j;A[L+i+1,L+i]+=1j
 K=np.zeros((1<<L,1<<L),complex)
 for a in range(2*L):
  for b in range(a+1,2*L):K+=.5*A[a,b]*(g[a]@g[b])
 return K,g

def indices(L,bond,q,maj):
 ix=2*q+bond
 if bond==0:return maj*L+ix,maj*L+(ix+1)%L
 return maj*L+(ix+1)%L,maj*L+ix

def hs_matrix(L,g,op):
 n=1<<L;U=np.eye(n,dtype=complex);V=op['local_V']
 if V==0:return U
 lam=np.arccosh(np.exp(.5*V*DT));ch=np.cosh(.5*lam);sh=np.sinh(.5*lam)
 for q,s in enumerate(op['s']):
  for maj in (0,1):
   a,b=indices(L,op['bond'],q,maj);U=(ch*np.eye(n)+1j*s*sh*(g[a]@g[b]))@U
 return U

def schedule(snapshot,L=6):
 by={x['index']:x for x in snapshot['operators']};last=max(by)+1;K,g=kinetic(L);trial=80;kh=expm(-.5*DT*K);kt=expm(-DT*K);ops=[]
 for j in range(last+1):
  if j in by:ops.append(hs_matrix(L,g,by[j]))
  else:ops.append(kt if j<trial else kh)
 return ops,g

def product(ops,start=0):
 dtype=np.clongdouble if ARGS.longdouble else complex;X=np.eye(ops[0].shape[0],dtype=dtype);logscale=0.;n=len(ops)
 for z in range(n):
  X=ops[(start+z)%n].astype(dtype)@X
  if z%8==7:
   sc=np.max(np.abs(X));X/=sc;logscale+=math.log(sc)
 return X,logscale

def weight(ops):
 X,ls=product(ops);return np.trace(X),ls

def green(ops,g,start):
 X,_=product(ops,start);w=np.trace(X);G=np.zeros((len(g),len(g)),complex)
 for i in range(len(g)):
  for j in range(len(g)):
   if i!=j:G[i,j]=np.trace((g[i]@g[j])@X)/w
 return G,G+np.eye(len(g))

def cplx(x):return complex(x[0],x[1])
def fast_after(G,op,ia,L=6):
 G=G.copy();s=op['s'][ia];th=np.tanh(np.arccosh(np.exp(.5*op['local_V']*DT)))
 for maj in (0,1):
  a,b=indices(L,op['bond'],ia,maj);tmp=1-1j*th*s*G[a,b];x1=-G[:,a].copy();x2=-G[:,b].copy();x1[a]+=2;x2[b]+=2;alpha=1j*s*th/tmp;G+=alpha*np.outer(x1,x2)-alpha*np.outer(x2,x1)
 return G

def check(path):
 s=json.load(open(path));global DT;DT=ARGS.dt;ops,g=schedule(s);w0,l0=weight(ops);target=next(x for x in s['operators'] if x['index']==s['operator_index']);Gex,gex=green(ops,g,s['operator_index']);gf=np.array([cplx(x) for x in s['G_fast']]).reshape(12,12);gfull=np.array([cplx(x) for x in s['G_full']]).reshape(12,12);gfa=fast_after(gf,target,s['aux_index']);target['s'][s['aux_index']]*=-1;ops1,_=schedule(s);w1,l1=weight(ops1);rex=w1/w0*np.exp(l1-l0);Gexa,gexa=green(ops1,g,s['operator_index']);target['s'][s['aux_index']]*=-1;rf,rfull=cplx(s['R_fast']),cplx(s['R_full'])
 return dict(file=str(path),flip=s['flip'],operator_index=s['operator_index'],region=s['region'],local_V=target['local_V'],bond=target['bond'],aux_index=s['aux_index'],trace_scaled_abs=abs(w0),R_fast_real=rf.real,R_fast_imag=rf.imag,R_full_real=rfull.real,R_full_imag=rfull.imag,R_exact_real=rex.real,R_exact_imag=rex.imag,R_fast_exact_abs_error=abs(rf-rex),R_full_exact_abs_error=abs(rfull-rex),R_fast_exact_magnitude_error=abs(abs(rf)-abs(rex)),R_full_exact_magnitude_error=abs(abs(rfull)-abs(rex)),G_fast_exact_max_error=float(np.max(np.abs(gf-gex))),G_full_exact_max_error=float(np.max(np.abs(gfull-gex))),G_fast_after_exact_max_error=float(np.max(np.abs(gfa-gexa))))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--dt',type=float,required=True);p.add_argument('--output',required=True);p.add_argument('--longdouble',action='store_true');p.add_argument('snapshots',nargs='+');ARGS=p.parse_args();rows=[check(Path(x)) for x in ARGS.snapshots]
 with open(ARGS.output,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print(json.dumps(rows,indent=2))
