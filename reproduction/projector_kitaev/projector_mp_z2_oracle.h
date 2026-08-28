#ifndef PROJECTOR_MP_Z2_ORACLE_H
#define PROJECTOR_MP_Z2_ORACLE_H

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "pfqmc.h"

// The physical full-contour convention is the production operator vector in
// its fixed index-zero order. Cyclic rotations are diagnostic views only.
namespace projector_mp_z2 {

template <unsigned Digits> struct Arithmetic {
    using Real=boost::multiprecision::number<boost::multiprecision::cpp_dec_float<Digits>>;
    struct C { Real r=0,i=0; C(){} C(Real x):r(x){} C(Real x,Real y):r(x),i(y){} };
    struct M { int n=0; std::vector<C>a; explicit M(int d=0):n(d),a(size_t(d)*d){} C&operator()(int x,int y){return a[size_t(x)*n+y];}const C&operator()(int x,int y)const{return a[size_t(x)*n+y];}static M eye(int n){M x(n);for(int i=0;i<n;++i)x(i,i)=C(Real(1));return x;} };
    struct Diagnostics { Real min_pivot=0,max_pivot=0;bool have_pivot=false;void pivot(const Real&v){if(!have_pivot){min_pivot=max_pivot=v;have_pivot=true;}else{min_pivot=std::min(min_pivot,v);max_pivot=std::max(max_pivot,v);}} };
    static C add(C a,const C&b){a.r+=b.r;a.i+=b.i;return a;}static C sub(C a,const C&b){a.r-=b.r;a.i-=b.i;return a;}static C neg(C a){a.r=-a.r;a.i=-a.i;return a;}static C mul(const C&a,const C&b){return{a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};}static C div(const C&a,const C&b){Real d=b.r*b.r+b.i*b.i;return{(a.r*b.r+a.i*b.i)/d,(a.i*b.r-a.r*b.i)/d};}static Real abs2(const C&a){return a.r*a.r+a.i*a.i;}static Real abs(const C&a){using boost::multiprecision::sqrt;return sqrt(abs2(a));}
    static M multiply(const M&a,const M&b){M c(a.n);for(int i=0;i<a.n;++i)for(int k=0;k<a.n;++k)for(int j=0;j<a.n;++j)c(i,j)=add(c(i,j),mul(a(i,k),b(k,j)));return c;}
    static M inverse(M a,Diagnostics&d){int n=a.n;M b=M::eye(n);for(int k=0;k<n;++k){int p=k;Real best=abs(a(k,k));for(int i=k+1;i<n;++i)if(abs(a(i,k))>best){best=abs(a(i,k));p=i;}if(best==0)throw std::runtime_error("MP Z2 oracle singular matrix");d.pivot(best);if(p!=k)for(int j=0;j<n;++j){std::swap(a(k,j),a(p,j));std::swap(b(k,j),b(p,j));}C q=a(k,k);for(int j=0;j<n;++j){a(k,j)=div(a(k,j),q);b(k,j)=div(b(k,j),q);}for(int i=0;i<n;++i)if(i!=k){C f=a(i,k);if(abs(f)==0)continue;for(int j=0;j<n;++j){a(i,j)=sub(a(i,j),mul(f,a(k,j)));b(i,j)=sub(b(i,j),mul(f,b(k,j)));}}}return b;}
    static M green(const M&b,Diagnostics&d){M x=b;for(int i=0;i<x.n;++i)x(i,i)=add(x(i,i),C(Real(1)));x=inverse(x,d);for(C&z:x.a){z.r*=2;z.i*=2;}for(int i=0;i<x.n;++i)x(i,i)=sub(x(i,i),C(Real(1)));return x;}
    static C pfaffian(M a,Diagnostics&d){C pf(Real(1));for(int k=0;k<a.n;k+=2){int p=k+1;Real best=abs(a(k,p));for(int j=k+2;j<a.n;++j)if(abs(a(k,j))>best){best=abs(a(k,j));p=j;}if(best==0)throw std::runtime_error("MP Z2 oracle zero Pfaffian pivot");d.pivot(best);if(p!=k+1){for(int j=0;j<a.n;++j)std::swap(a(k+1,j),a(p,j));for(int i=0;i<a.n;++i)std::swap(a(i,k+1),a(i,p));pf=neg(pf);}C pivot=a(k,k+1);pf=mul(pf,pivot);for(int i=k+2;i<a.n;++i)for(int j=i+1;j<a.n;++j){C correction=div(sub(mul(a(k,i),a(k+1,j)),mul(a(k,j),a(k+1,i))),pivot);a(i,j)=sub(a(i,j),correction);a(j,i)=neg(a(i,j));}}return pf;}
    static C phase(C z){Real m=abs(z);if(m==0)throw std::runtime_error("MP Z2 oracle zero phase");return div(z,C(m));}
    static M factor(Operator*op,int n){MatType id=MatType::Identity(n,n),b;op->left_multiply(id,b);M x(n);for(int i=0;i<n;++i)for(int j=0;j<n;++j)x(i,j)=C(Real(b(i,j).real()),Real(b(i,j).imag()));return x;}
    static C weightPhase(Operator*op){DataType z=op->getSignOfWeight();return phase(C(Real(z.real()),Real(z.imag())));}
};

template <unsigned Digits> MpZ2Result fullContourZ2AtPrecision(const std::vector<Operator*>&ops){
    using A=Arithmetic<Digits>;using Real=typename A::Real;using C=typename A::C;using M=typename A::M;typename A::Diagnostics diag;MpZ2Result result;result.precision_digits=int(Digits);result.canonical_order=true;
    try{if(ops.empty())throw std::runtime_error("MP Z2 oracle empty contour");MatType id;ops[0]->getGreensMat(id);int n=int(id.rows());M product=A::factor(ops[0],n),gcur=A::green(product,diag);C sign=A::weightPhase(ops[0]);C extra=((n/2)%2==0)?C(Real(1)):C(Real(-1));for(int i=1;i<int(ops.size());++i){M next=A::factor(ops[i],n),gnext=A::green(next,diag),block(2*n);for(int r=0;r<n;++r)for(int c=0;c<n;++c){block(r,c)=gnext(r,c);block(n+r,n+c)=gcur(r,c);}for(int r=0;r<n;++r){block(r,n+r)=C(Real(-1));block(n+r,r)=C(Real(1));}sign=A::mul(sign,A::mul(A::weightPhase(ops[i]),A::mul(A::phase(A::pfaffian(block,diag)),extra)));if(i+1<int(ops.size())){product=A::multiply(next,product);gcur=A::green(product,diag);}}sign=A::phase(sign);Real reality=A::abs(sign.i)/std::max(A::abs(sign.r),Real("1e-100000"));result.reality_error=reality.template convert_to<double>();result.z2=sign.r<0?-1:1;Real logCondition=0;if(diag.have_pivot&&diag.min_pivot>0){using boost::multiprecision::log10;logCondition=log10(diag.max_pivot/diag.min_pivot);}result.residual_or_condition=logCondition.template convert_to<double>();if(!std::isfinite(result.reality_error)||result.reality_error>1e-30){result.status=MpZ2Status::untrusted_reality;result.message="phase is not resolved as real";return result;}if(!std::isfinite(result.residual_or_condition)||result.residual_or_condition>1e12){result.status=MpZ2Status::untrusted_condition;result.message="nonfinite/catastrophic pivot conditioning diagnostic";return result;}result.status=MpZ2Status::untrusted_precision;result.message="single precision result requires agreement at a second precision";return result;}catch(const std::exception&e){result.status=MpZ2Status::unavailable;result.z2=0;result.message=e.what();return result;}}

inline bool numericallyResolved(const MpZ2Result&r){return r.z2!=0&&r.status==MpZ2Status::untrusted_precision;}
inline MpZ2Result fullContourZ2(const std::vector<Operator*>&ops){const MpZ2Result r160=fullContourZ2AtPrecision<160>(ops);const MpZ2Result r320=fullContourZ2AtPrecision<320>(ops);if(numericallyResolved(r160)&&numericallyResolved(r320)&&r160.z2==r320.z2){MpZ2Result out=r320;out.status=MpZ2Status::trusted;out.converged=true;out.precision_escalations=1;out.message="canonical MP160/MP320 agreement";return out;}const MpZ2Result r640=fullContourZ2AtPrecision<640>(ops);if(numericallyResolved(r320)&&numericallyResolved(r640)&&r320.z2==r640.z2){MpZ2Result out=r640;out.status=MpZ2Status::trusted;out.converged=true;out.precision_escalations=2;out.message="canonical MP320/MP640 agreement after escalation";return out;}MpZ2Result out=r640;out.status=MpZ2Status::untrusted_precision;out.converged=false;out.precision_escalations=2;std::ostringstream msg;msg<<"canonical MP Z2 failed consecutive-precision agreement: 160="<<r160.z2<<", 320="<<r320.z2<<", 640="<<r640.z2;out.message=msg.str();return out;}

}
#endif
