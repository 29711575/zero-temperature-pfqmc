// Independent arbitrary-boundary Majorana Green rebuild for frozen snapshots.
// Debug/validation helper only: no production PfQMC code is used or modified.
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

using boost::property_tree::ptree;

template <class R> struct C {
    R re{}, im{};
    C() = default;
    C(const R &r) : re(r), im(0) {}
    C(const R &r, const R &i) : re(r), im(i) {}
    C operator-() const { return C(-re, -im); }
    C &operator+=(const C &z) { re += z.re; im += z.im; return *this; }
    C &operator-=(const C &z) { re -= z.re; im -= z.im; return *this; }
    C &operator*=(const C &z) {
        R r = re * z.re - im * z.im;
        im = re * z.im + im * z.re; re = r; return *this;
    }
    C &operator/=(const C &z) {
        R d = z.re * z.re + z.im * z.im;
        R r = (re * z.re + im * z.im) / d;
        im = (im * z.re - re * z.im) / d; re = r; return *this;
    }
};
template <class R> C<R> operator+(C<R> a, const C<R> &b) { return a += b; }
template <class R> C<R> operator-(C<R> a, const C<R> &b) { return a -= b; }
template <class R> C<R> operator*(C<R> a, const C<R> &b) { return a *= b; }
template <class R> C<R> operator/(C<R> a, const C<R> &b) { return a /= b; }
template <class R> R abs2(const C<R> &z) { return z.re*z.re + z.im*z.im; }
template <class R> R absv(const C<R> &z) { using boost::multiprecision::sqrt; return sqrt(abs2(z)); }

template <class R> struct Matrix {
    int n = 0;
    std::vector<C<R>> a;
    Matrix() = default;
    explicit Matrix(int d) : n(d), a(std::size_t(d)*d) {}
    C<R> &operator()(int i, int j) { return a[std::size_t(i)*n+j]; }
    const C<R> &operator()(int i, int j) const { return a[std::size_t(i)*n+j]; }
    static Matrix identity(int n) { Matrix x(n); for (int i=0;i<n;++i) x(i,i)=C<R>(R(1)); return x; }
};

template <class R> Matrix<R> multiply(const Matrix<R> &x, const Matrix<R> &y) {
    Matrix<R> z(x.n);
    for (int i=0;i<x.n;++i) for (int k=0;k<x.n;++k) {
        const C<R> xik=x(i,k);
        for (int j=0;j<x.n;++j) z(i,j) += xik*y(k,j);
    }
    return z;
}
template <class R> R max_abs(const Matrix<R> &x) {
    R m=0; for (const auto &z:x.a) m=std::max(m,absv(z)); return m;
}
template <class R> R inf_norm(const Matrix<R> &x) {
    R m=0; for(int i=0;i<x.n;++i){R s=0;for(int j=0;j<x.n;++j)s+=absv(x(i,j));m=std::max(m,s);}return m;
}
template <class R> void scale_inplace(Matrix<R> &x,const R &s){for(auto &z:x.a){z.re*=s;z.im*=s;}}

// Scaling/squaring Taylor exponential.  Contour kinetic steps have small norm,
// so this converges rapidly and avoids any double-precision linear algebra.
template <class R> Matrix<R> expm(Matrix<R> x, int digits) {
    int squarings=0; R norm=inf_norm(x);
    while(norm>R("0.25")){scale_inplace(x,R("0.5"));norm*=R("0.5");++squarings;}
    Matrix<R> sum=Matrix<R>::identity(x.n), term=sum;
    const R tol=pow(R(10),R(-(digits-8)));
    for(int k=1;k<8*digits;++k){
        term=multiply(x,term); scale_inplace(term,R(1)/R(k));
        for(std::size_t q=0;q<sum.a.size();++q)sum.a[q]+=term.a[q];
        if(max_abs(term)<tol)break;
        if(k+1==8*digits)throw std::runtime_error("matrix exponential did not converge");
    }
    for(int k=0;k<squarings;++k)sum=multiply(sum,sum);
    return sum;
}

template <class R> Matrix<R> inverse_partial_pivot(Matrix<R> a, R *pivot_spread=nullptr) {
    const int n=a.n; Matrix<R> b=Matrix<R>::identity(n); R pmax=0,pmin=std::numeric_limits<R>::max();
    for(int k=0;k<n;++k){
        int p=k;R best=abs2(a(k,k));
        for(int i=k+1;i<n;++i)if(abs2(a(i,k))>best){best=abs2(a(i,k));p=i;}
        if(best==0)throw std::runtime_error("singular multiprecision denominator");
        if(p!=k)for(int j=0;j<n;++j){std::swap(a(k,j),a(p,j));std::swap(b(k,j),b(p,j));}
        R pv=absv(a(k,k));pmax=std::max(pmax,pv);pmin=std::min(pmin,pv);
        C<R> d=a(k,k);for(int j=0;j<n;++j){a(k,j)/=d;b(k,j)/=d;}
        for(int i=0;i<n;++i)if(i!=k){C<R> f=a(i,k);if(abs2(f)==0)continue;
            for(int j=0;j<n;++j){a(i,j)-=f*a(k,j);b(i,j)-=f*b(k,j);}}
    }
    if(pivot_spread)*pivot_spread=pmax/pmin;return b;
}

struct Op { int index=0,bond=0; std::string local_v; std::vector<int> s; };
struct Snapshot { int flip=0,boundary=0,aux=0,L=6; std::vector<Op> ops; std::vector<std::pair<std::string,std::string>> production; };

static std::string shortest_roundtrip(double x) {
    char buf[64];
    for (int p=1;p<=17;++p) {
        std::snprintf(buf,sizeof(buf),"%.*g",p,x);
        const double y=std::strtod(buf,nullptr);
        if (std::memcmp(&x,&y,sizeof(double))==0) return buf;
    }
    throw std::runtime_error("could not format binary64 value");
}

Snapshot read_snapshot(const std::string &path, bool accepted) {
    ptree root; boost::property_tree::read_json(path,root); Snapshot s;
    s.flip=root.get<int>("flip");s.boundary=root.get<int>("operator_index");s.aux=root.get<int>("aux_index");
    for(const auto &node:root.get_child("operators")){
        Op o;o.index=node.second.get<int>("index");o.bond=node.second.get<int>("bond");
        // Match the frozen Python referee: JSON is first decoded as binary64,
        // then converted to its compact physical decimal schedule value.
        o.local_v=shortest_roundtrip(node.second.get<double>("local_V"));
        for(const auto &v:node.second.get_child("s"))o.s.push_back(v.second.get_value<int>());
        s.ops.push_back(o);
    }
    if(accepted)for(auto &o:s.ops)if(o.index==s.boundary){o.s.at(s.aux)*=-1;break;}
    if(auto child=root.get_child_optional("G_accepted"))for(const auto &v:*child){
        auto it=v.second.begin();std::string re=(it++)->second.get_value<std::string>();std::string im=it->second.get_value<std::string>();s.production.emplace_back(re,im);
    }
    return s;
}

template <class R> Matrix<R> kinetic_generator(int L) {
    Matrix<R> a(2*L); const C<R> ip(R(0),R(1)), im(R(0),R(-1));
    for(int i=0;i<L-1;++i){
        C<R> z=(i%2==0)?ip:im;
        for(int k=0;k<2;++k){int x=k*L+i,y=x+1;a(x,y)=z;a(y,x)=-z;}
        a(i,i+1)+=ip;a(i+1,i)+=im;a(L+i,L+i+1)+=im;a(L+i+1,L+i)+=ip;
    }
    return a;
}
static std::pair<int,int> aux_indices(int L,int bond,int q,int maj){int ix=2*q+bond;return bond==0?std::make_pair(maj*L+ix,maj*L+(ix+1)%L):std::make_pair(maj*L+(ix+1)%L,maj*L+ix);}

template <class R> Matrix<R> hs_matrix(int L,const Op &o,const R &dt){
    Matrix<R> b=Matrix<R>::identity(2*L);R V(o.local_v);if(V==0)return b;
    using boost::multiprecision::exp;using boost::multiprecision::acosh;using boost::multiprecision::cosh;using boost::multiprecision::sinh;
    R lam=acosh(exp(R("0.5")*V*dt)),ch=cosh(lam),sh=sinh(lam);
    for(int q=0;q<(int)o.s.size();++q)for(int maj=0;maj<2;++maj){auto p=aux_indices(L,o.bond,q,maj);b(p.first,p.first)=C<R>(ch);b(p.second,p.second)=C<R>(ch);b(p.first,p.second)=C<R>(0,sh*R(o.s[q]));b(p.second,p.first)=C<R>(0,-sh*R(o.s[q]));}
    return b;
}

template <class R> struct Rebuild { Matrix<R> g; R pivot_spread, log_scale, residual; double seconds=0; };
template <class R> Rebuild<R> rebuild(const Snapshot &s,int digits) {
    using boost::multiprecision::exp;using boost::multiprecision::log;
    auto t0=std::chrono::steady_clock::now();const int n=2*s.L;const R dt("0.1");
    Matrix<R> K=kinetic_generator<R>(s.L);scale_inplace(K,R(-dt));Matrix<R> kt=expm(K,digits);
    Matrix<R> H=kinetic_generator<R>(s.L);scale_inplace(H,-dt*R("0.5"));Matrix<R> kh=expm(H,digits);
    int last=0;for(const auto&o:s.ops)last=std::max(last,o.index+1);const int N=last+1;
    std::vector<const Op*> by(N,nullptr);for(const auto&o:s.ops)by[o.index]=&o;
    Matrix<R> p=Matrix<R>::identity(n);R logscale=0;
    for(int off=0;off<N;++off){int idx=(s.boundary+off)%N;Matrix<R> b=by[idx]?hs_matrix<R>(s.L,*by[idx],dt):(idx<80?kt:kh);p=multiply(b,p);R sc=max_abs(p);scale_inplace(p,R(1)/sc);logscale+=log(sc);}
    R eps=exp(-logscale);Matrix<R> den=p;for(int i=0;i<n;++i)den(i,i)+=C<R>(eps);
    R piv;Matrix<R> inv=inverse_partial_pivot(den,&piv);Matrix<R> g=inv;scale_inplace(g,R(2)*eps);
    Matrix<R> check=multiply(den,inv);for(int i=0;i<n;++i)check(i,i)-=C<R>(R(1));R residual=inf_norm(check);
    auto t1=std::chrono::steady_clock::now();return {g,piv,logscale,residual,std::chrono::duration<double>(t1-t0).count()};
}

template <class R> R production_error(const Matrix<R>&g,const Snapshot&s){if(s.production.size()!=g.a.size())return R(-1);R m=0;for(std::size_t i=0;i<g.a.size();++i){C<R> z(R(s.production[i].first),R(s.production[i].second));m=std::max(m,absv(g.a[i]-z));}return m;}
template <class R> void emit_matrix(std::ostream&os,const Matrix<R>&g){os<<"[";for(std::size_t i=0;i<g.a.size();++i){if(i)os<<",";os<<"[\""<<g.a[i].re<<"\",\""<<g.a[i].im<<"\"]";}os<<"]";}

template <unsigned Digits> void run_one(const Snapshot&s,const std::string&label){
    using R=boost::multiprecision::number<boost::multiprecision::cpp_dec_float<Digits>>;
    auto x=rebuild<R>(s,Digits);std::cout<<std::setprecision(Digits);
    std::cout<<"{\"label\":\""<<label<<"\",\"flip\":"<<s.flip<<",\"L\":"<<s.L<<",\"digits\":"<<Digits
             <<",\"seconds\":"<<std::setprecision(9)<<x.seconds<<std::setprecision(Digits)
             <<",\"production_max_error\":\""<<production_error(x.g,s)<<"\",\"pivot_spread\":\""<<x.pivot_spread
             <<"\",\"log_scale\":\""<<x.log_scale<<"\",\"inverse_residual_inf\":\""<<x.residual<<"\",\"green\":";emit_matrix(std::cout,x.g);std::cout<<"}\n";
}

int main(int argc,char**argv){
    try{if(argc<3)throw std::runtime_error("usage: multiprecision_green_check snapshot.json digits [L] [accepted=1] [label]");
        int digits=std::stoi(argv[2]);Snapshot s=read_snapshot(argv[1],argc<5||std::stoi(argv[4])!=0);if(argc>=4)s.L=std::stoi(argv[3]);
        // Runtime-only representative scaling case: deterministically extend
        // each frozen half-chain HS field vector when L exceeds the snapshot L.
        // L=6 frozen validation is unchanged.
        for(auto &o:s.ops)for(int q=(int)o.s.size();q<s.L/2;++q)o.s.push_back(((o.index+3*q+o.bond)&1)?-1:1);
        std::string label=argc>=6?argv[5]:argv[1];
        if(digits==50)run_one<50>(s,label);else if(digits==80)run_one<80>(s,label);else if(digits==100)run_one<100>(s,label);else throw std::runtime_error("digits must be 50, 80, or 100");
    }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}return 0;
}
