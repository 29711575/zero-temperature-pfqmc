#ifndef DRIVEN_MULTIPRECISION_REBUILD_H
#define DRIVEN_MULTIPRECISION_REBUILD_H

// Debug-only 50-decimal rebuild for the driven contour.  Factors are rebuilt
// from the live HS fields and chain parameters; casting existing double B
// matrices is intentionally not used.
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include "kitaevChain.h"

namespace driven_multiprecision {
using Real = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<80>>;
struct Complex { Real r=0,i=0; Complex(){} Complex(Real a):r(a){} Complex(Real a,Real b):r(a),i(b){} };
inline Complex operator+(Complex a,const Complex&b){a.r+=b.r;a.i+=b.i;return a;}
inline Complex operator-(Complex a,const Complex&b){a.r-=b.r;a.i-=b.i;return a;}
inline Complex operator-(Complex a){a.r=-a.r;a.i=-a.i;return a;}
inline Complex operator*(const Complex&a,const Complex&b){return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};}
inline Complex operator/(const Complex&a,const Complex&b){Real d=b.r*b.r+b.i*b.i;return {(a.r*b.r+a.i*b.i)/d,(a.i*b.r-a.r*b.i)/d};}
inline Real abs(const Complex&a){using boost::multiprecision::sqrt;return sqrt(a.r*a.r+a.i*a.i);}
struct Matrix {int n;std::vector<Complex>a;explicit Matrix(int d=0):n(d),a(size_t(d)*d){}Complex&operator()(int x,int y){return a[size_t(x)*n+y];}const Complex&operator()(int x,int y)const{return a[size_t(x)*n+y];}static Matrix eye(int n){Matrix x(n);for(int i=0;i<n;++i)x(i,i)=Complex(Real(1));return x;}};
inline Matrix mul(const Matrix&a,const Matrix&b){Matrix c(a.n);for(int i=0;i<a.n;++i)for(int k=0;k<a.n;++k)for(int j=0;j<a.n;++j)c(i,j)=c(i,j)+a(i,k)*b(k,j);return c;}
inline Real maxabs(const Matrix&a){Real x=0;for(auto&z:a.a)x=std::max(x,abs(z));return x;}
inline Real norminf(const Matrix&a){Real x=0;for(int i=0;i<a.n;++i){Real y=0;for(int j=0;j<a.n;++j)y+=abs(a(i,j));x=std::max(x,y);}return x;}
inline void scale(Matrix&a,const Real&s){for(auto&z:a.a){z.r*=s;z.i*=s;}}
inline Matrix expm(Matrix a){using boost::multiprecision::pow;Real n=norminf(a);int q=0;while(n>Real("0.25")){scale(a,Real("0.5"));n*=Real("0.5");++q;}Matrix sum=Matrix::eye(a.n),term=sum;Real tol=pow(Real(10),Real(-70));for(int k=1;k<640;++k){term=mul(a,term);scale(term,Real(1)/Real(k));for(size_t z=0;z<sum.a.size();++z)sum.a[z]=sum.a[z]+term.a[z];if(maxabs(term)<tol)break;}for(int k=0;k<q;++k)sum=mul(sum,sum);return sum;}
inline Matrix inverse(Matrix a){int n=a.n;Matrix b=Matrix::eye(n);for(int k=0;k<n;++k){int p=k;Real best=abs(a(k,k));for(int i=k+1;i<n;++i)if(abs(a(i,k))>best){best=abs(a(i,k));p=i;}if(best==0)throw std::runtime_error("singular MP rebuild");for(int j=0;j<n;++j){std::swap(a(k,j),a(p,j));std::swap(b(k,j),b(p,j));}Complex d=a(k,k);for(int j=0;j<n;++j){a(k,j)=a(k,j)/d;b(k,j)=b(k,j)/d;}for(int i=0;i<n;++i)if(i!=k){Complex f=a(i,k);for(int j=0;j<n;++j){a(i,j)=a(i,j)-f*a(k,j);b(i,j)=b(i,j)-f*b(k,j);}}}return b;}
inline std::string compact_double(double x){char b[64];for(int p=1;p<=17;++p){std::snprintf(b,sizeof(b),"%.*g",p,x);double y=std::strtod(b,nullptr);if(std::memcmp(&x,&y,sizeof(double))==0)return b;}return "0";}
inline Matrix kinetic(const SpinlessTvChainUtils&c){int L=c.Lx;Matrix h(2*L);Complex I(0,1),id(0,c.delta),im(0,c.mu);for(int x=0;x<L-1;++x){Complex z=x%2? -I:I;for(int k=0;k<2;++k){int a=k*L+x,b=a+1;h(a,b)=z;h(b,a)=-z;}int a=x,b=x+1;h(a,b)=h(a,b)+id;h(b,a)=h(b,a)-id;a=L+x;b=a+1;h(a,b)=h(a,b)-id;h(b,a)=h(b,a)+id;}for(int x=0;x<L;++x){int a=x,b=L+x;h(a,b)=h(a,b)-im;h(b,a)=h(b,a)+im;}return h;}
inline Matrix hs(const SpinlessTvChainUtils&c,const SpinlessVOperator&v){using boost::multiprecision::exp;using boost::multiprecision::acosh;using boost::multiprecision::cosh;using boost::multiprecision::sinh;int n=c.nDim;Matrix b=Matrix::eye(n);Real V(compact_double(v.localV)),dt(compact_double(c.dt)),lam=acosh(exp(Real("0.5")*V*dt)),ch=cosh(lam),sh=sinh(lam);for(int q=0;q<v.s->size();++q)for(int m=0;m<2;++m){int a,d;c.aux2MajoranaIdx(q,m,v.bondType,a,d);b(a,a)=ch;b(d,d)=ch;b(a,d)=Complex(0,sh*Real((*v.s)(q)));b(d,a)=Complex(0,-sh*Real((*v.s)(q)));}return b;}
inline bool rebuild(const SpinlessTvChainUtils&c,const std::vector<Operator*>&ops,int boundary,int trial_dense,MatType&out){try{int n=c.nDim,N=ops.size();if(!N)return false;Real dt(compact_double(c.dt));Matrix h=kinetic(c),full=h,half=h;scale(full,-dt);scale(half,-dt*Real("0.5"));Matrix kt=expm(full),kh=expm(half),p=Matrix::eye(n);Real ls=0;for(int off=0;off<N;++off){int j=(boundary+off)%N;Matrix b;if(auto*v=dynamic_cast<SpinlessVOperator*>(ops[j]))b=hs(c,*v);else b=(j<trial_dense?kt:kh);p=mul(b,p);Real s=maxabs(p);scale(p,Real(1)/s);using boost::multiprecision::log;ls+=log(s);}using boost::multiprecision::exp;Real e=exp(-ls);for(int i=0;i<n;++i)p(i,i)=p(i,i)+e;Matrix g=inverse(p);scale(g,Real(2)*e);out.resize(n,n);for(int i=0;i<n;++i)for(int j=0;j<n;++j)out(i,j)=DataType(g(i,j).r.convert_to<double>(),g(i,j).i.convert_to<double>());return true;}catch(...){return false;}}
}
#endif
