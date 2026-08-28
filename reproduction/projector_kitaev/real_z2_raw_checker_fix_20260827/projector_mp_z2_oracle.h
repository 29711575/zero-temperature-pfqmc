#ifndef PROJECTOR_MP_Z2_ORACLE_H
#define PROJECTOR_MP_Z2_ORACLE_H

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include "operator.h"

// Read-only 160-digit full-contour oracle.  It consumes no RNG and constructs
// every one-body factor through Operator::left_multiply, so its contour order
// is exactly the production order.  Intended for one-time real-Z2
// initialization and rare trusted-mismatch adjudication only.
namespace projector_mp_z2 {
using Real=boost::multiprecision::number<boost::multiprecision::cpp_dec_float<160>>;
struct C { Real r=0,i=0; C(){} C(Real x):r(x){} C(Real x,Real y):r(x),i(y){} };
struct M {
    int n=0; std::vector<C> a;
    explicit M(int d=0):n(d),a(size_t(d)*d){}
    C& operator()(int x,int y){return a[size_t(x)*n+y];}
    const C& operator()(int x,int y)const{return a[size_t(x)*n+y];}
    static M eye(int n){M x(n);for(int i=0;i<n;++i)x(i,i)=C(1);return x;}
};
inline C add(C a,const C&b){a.r+=b.r;a.i+=b.i;return a;}
inline C sub(C a,const C&b){a.r-=b.r;a.i-=b.i;return a;}
inline C neg(C a){a.r=-a.r;a.i=-a.i;return a;}
inline C mul(const C&a,const C&b){return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};}
inline C div(const C&a,const C&b){Real d=b.r*b.r+b.i*b.i;return {(a.r*b.r+a.i*b.i)/d,(a.i*b.r-a.r*b.i)/d};}
inline Real abs2(const C&a){return a.r*a.r+a.i*a.i;}
inline Real abs(const C&a){using boost::multiprecision::sqrt;return sqrt(abs2(a));}
inline M multiply(const M&a,const M&b){M c(a.n);for(int i=0;i<a.n;++i)for(int k=0;k<a.n;++k)for(int j=0;j<a.n;++j)c(i,j)=add(c(i,j),mul(a(i,k),b(k,j)));return c;}
inline M inverse(M a){const int n=a.n;M b=M::eye(n);for(int k=0;k<n;++k){int p=k;Real best=abs(a(k,k));for(int i=k+1;i<n;++i)if(abs(a(i,k))>best){best=abs(a(i,k));p=i;}if(best==0)throw std::runtime_error("MP Z2 oracle singular matrix");if(p!=k)for(int j=0;j<n;++j){std::swap(a(k,j),a(p,j));std::swap(b(k,j),b(p,j));}C d=a(k,k);for(int j=0;j<n;++j){a(k,j)=div(a(k,j),d);b(k,j)=div(b(k,j),d);}for(int i=0;i<n;++i)if(i!=k){C f=a(i,k);if(abs(f)==0)continue;for(int j=0;j<n;++j){a(i,j)=sub(a(i,j),mul(f,a(k,j)));b(i,j)=sub(b(i,j),mul(f,b(k,j)));}}}return b;}
inline M green(const M&b){M x=b;for(int i=0;i<x.n;++i)x(i,i)=add(x(i,i),C(1));x=inverse(x);for(C&z:x.a){z.r*=2;z.i*=2;}for(int i=0;i<x.n;++i)x(i,i)=sub(x(i,i),C(1));return x;}
inline C pfaffian(M a){C pf(1);for(int k=0;k<a.n;k+=2){int p=k+1;Real best=abs(a(k,p));for(int j=k+2;j<a.n;++j)if(abs(a(k,j))>best){best=abs(a(k,j));p=j;}if(best==0)throw std::runtime_error("MP Z2 oracle zero Pfaffian pivot");if(p!=k+1){for(int j=0;j<a.n;++j)std::swap(a(k+1,j),a(p,j));for(int i=0;i<a.n;++i)std::swap(a(i,k+1),a(i,p));pf=neg(pf);}C pivot=a(k,k+1);pf=mul(pf,pivot);for(int i=k+2;i<a.n;++i)for(int j=i+1;j<a.n;++j){C correction=div(sub(mul(a(k,i),a(k+1,j)),mul(a(k,j),a(k+1,i))),pivot);a(i,j)=sub(a(i,j),correction);a(j,i)=neg(a(i,j));}}return pf;}
inline C phase(C z){Real m=abs(z);if(m==0)throw std::runtime_error("MP Z2 oracle zero phase");return div(z,C(m));}
inline M factor(Operator*op,int n){MatType id=MatType::Identity(n,n),b;op->left_multiply(id,b);M x(n);for(int i=0;i<n;++i)for(int j=0;j<n;++j)x(i,j)=C(Real(b(i,j).real()),Real(b(i,j).imag()));return x;}
inline C weightPhase(Operator*op){DataType z=op->getSignOfWeight();return phase(C(Real(z.real()),Real(z.imag())));}
inline int fullContourZ2(const std::vector<Operator*>&ops){
    if(ops.empty())throw std::runtime_error("MP Z2 oracle empty contour");
    MatType probe=MatType::Identity(1,1),unused; (void)probe; (void)unused;
    MatType id; // n is recovered from the first factor by probing its auxiliary Green.
    ops[0]->getGreensMat(id); const int n=int(id.rows());
    M product=factor(ops[0],n),gcur=green(product);C sign=weightPhase(ops[0]);
    const C extra=((n/2)%2==0)?C(1):C(-1);
    for(int i=1;i<int(ops.size());++i){M next=factor(ops[i],n),gnext=green(next),block(2*n);for(int r=0;r<n;++r)for(int c=0;c<n;++c){block(r,c)=gnext(r,c);block(n+r,n+c)=gcur(r,c);}for(int r=0;r<n;++r){block(r,n+r)=C(-1);block(n+r,r)=C(1);}sign=mul(sign,mul(weightPhase(ops[i]),mul(phase(pfaffian(block)),extra)));if(i+1<int(ops.size())){product=multiply(next,product);gcur=green(product);}}
    sign=phase(sign);const Real tol("1e-40");if(abs(sign.i)>tol||abs(sign.r)<Real("0.5"))throw std::runtime_error("MP Z2 oracle did not resolve a real sign");return sign.r<0?-1:1;
}
}
#endif
