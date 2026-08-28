#include <algorithm>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"
#include "projector_mp_z2_oracle.h"

namespace {

class Walker : public Spinless_tV {
  public:
    int center=-1,ntrial=0,nphys=0;
    Walker(const SpinlessTvChainUtils*c,rdGenerator*r,double theta,double beta){
        build_projector_static_contour(*this,c,r,theta,beta,center,ntrial,nphys);
    }
};

void burnAdvance(PfQMC&q,int burn,int target){
    for(int i=0;i<burn;++i){q.rightSweep();q.leftSweep();}
    for(int k=0;k<target;++k){q.rightSweep();q.leftSweep();}
}

std::string compact(double x){
    char b[64];
    for(int p=1;p<=17;++p){
        std::snprintf(b,sizeof(b),"%.*g",p,x);
        double y=std::strtod(b,nullptr);
        if(std::memcmp(&x,&y,sizeof(double))==0)return b;
    }
    throw std::runtime_error("cannot round-trip double");
}

template<unsigned Digits> struct MP {
    using Real=boost::multiprecision::number<boost::multiprecision::cpp_dec_float<Digits>>;
    struct C{Real r=0,i=0;C(){}C(Real x):r(x){}C(Real x,Real y):r(x),i(y){}};
    struct M{
        int n=0;std::vector<C>a;
        explicit M(int d=0):n(d),a(size_t(d)*d){}
        C&operator()(int x,int y){return a[size_t(x)*n+y];}
        const C&operator()(int x,int y)const{return a[size_t(x)*n+y];}
        static M eye(int n){M x(n);for(int i=0;i<n;++i)x(i,i)=C(1);return x;}
    };
    struct Solve{M inverse;C determinant;Real residual=0,condition=0;};
    static C add(C a,const C&b){a.r+=b.r;a.i+=b.i;return a;}
    static C sub(C a,const C&b){a.r-=b.r;a.i-=b.i;return a;}
    static C neg(C a){a.r=-a.r;a.i=-a.i;return a;}
    static C mulc(const C&a,const C&b){return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};}
    static C divc(const C&a,const C&b){Real d=b.r*b.r+b.i*b.i;return {(a.r*b.r+a.i*b.i)/d,(a.i*b.r-a.r*b.i)/d};}
    static C conjc(C a){a.i=-a.i;return a;}
    static Real abs2(const C&a){return a.r*a.r+a.i*a.i;}
    static Real absc(const C&a){using boost::multiprecision::sqrt;return sqrt(abs2(a));}
    static M mulm(const M&a,const M&b){M c(a.n);for(int i=0;i<a.n;++i)for(int k=0;k<a.n;++k){C z=a(i,k);for(int j=0;j<a.n;++j)c(i,j)=add(c(i,j),mulc(z,b(k,j)));}return c;}
    static M subm(const M&a,const M&b){M c(a.n);for(size_t k=0;k<c.a.size();++k)c.a[k]=sub(a.a[k],b.a[k]);return c;}
    static void scale(M&a,const Real&s){for(C&z:a.a){z.r*=s;z.i*=s;}}
    static Real maxabs(const M&a){Real x=0;for(const C&z:a.a)x=std::max(x,absc(z));return x;}
    static Real norminf(const M&a){Real x=0;for(int i=0;i<a.n;++i){Real y=0;for(int j=0;j<a.n;++j)y+=absc(a(i,j));x=std::max(x,y);}return x;}
    static Real frob(const M&a){using boost::multiprecision::sqrt;Real x=0;for(const C&z:a.a)x+=abs2(z);return sqrt(x);}
    static M expm(M a){
        using boost::multiprecision::pow;
        Real n=norminf(a);int q=0;
        while(n>Real("0.25")){scale(a,Real("0.5"));n*=Real("0.5");++q;}
        M sum=M::eye(a.n),term=sum;
        const Real tol=pow(Real(10),-int(Digits-12));
        for(int k=1;k<4096;++k){
            term=mulm(a,term);scale(term,Real(1)/Real(k));
            for(size_t z=0;z<sum.a.size();++z)sum.a[z]=add(sum.a[z],term.a[z]);
            if(maxabs(term)<tol)break;
            if(k==4095)throw std::runtime_error("MP expm did not converge");
        }
        for(int k=0;k<q;++k)sum=mulm(sum,sum);
        return sum;
    }
    static Solve solve(M a){
        const M original=a;int n=a.n;M b=M::eye(n);C det(1);int parity=1;
        for(int k=0;k<n;++k){
            int p=k;Real best=absc(a(k,k));
            for(int i=k+1;i<n;++i)if(absc(a(i,k))>best){best=absc(a(i,k));p=i;}
            if(best==0)throw std::runtime_error("singular MP matrix");
            if(p!=k){for(int j=0;j<n;++j){std::swap(a(k,j),a(p,j));std::swap(b(k,j),b(p,j));}parity=-parity;}
            C d=a(k,k);det=mulc(det,d);
            for(int j=0;j<n;++j){a(k,j)=divc(a(k,j),d);b(k,j)=divc(b(k,j),d);}
            for(int i=0;i<n;++i)if(i!=k){C f=a(i,k);if(absc(f)==0)continue;for(int j=0;j<n;++j){a(i,j)=sub(a(i,j),mulc(f,a(k,j)));b(i,j)=sub(b(i,j),mulc(f,b(k,j)));}}
        }
        if(parity<0)det=neg(det);
        M residual=subm(mulm(original,b),M::eye(n));
        Solve out;out.inverse=std::move(b);out.determinant=det;
        out.condition=norminf(original)*norminf(out.inverse);
        out.residual=norminf(residual)/std::max(out.condition,Real(1));
        return out;
    }
    static M green(const M&product,Solve*solveOut=nullptr){
        M a=product;for(int i=0;i<a.n;++i)a(i,i)=add(a(i,i),C(1));
        Solve s=solve(a);M g=s.inverse;scale(g,Real(2));for(int i=0;i<g.n;++i)g(i,i)=sub(g(i,i),C(1));if(solveOut)*solveOut=std::move(s);return g;
    }
    static Real rel(const M&a,const M&b){return frob(subm(a,b))/std::max(frob(b),Real(1e-100));}
    // The oracle's Cayley Green is Gc=2(I+B)^-1-I, while production stores
    // Gp=Gc+I.  Thus Gp+Gp^T-2I = Gc+Gc^T.  Report the production diagnostic
    // normalization requested by the audit without changing ratio evaluation.
    static Real structure(const M&gc){M numerator(gc.n),gp=gc;for(int i=0;i<gc.n;++i){gp(i,i)=add(gp(i,i),C(1));for(int j=0;j<gc.n;++j)numerator(i,j)=add(gc(i,j),gc(j,i));}return frob(numerator)/std::max(frob(gp),Real("1e-100"));}
    static std::string str(const Real&x){std::ostringstream s;s<<std::setprecision(std::min<unsigned>(Digits-5,80))<<std::scientific<<x;return s.str();}
    static std::string strc(const C&x){return str(x.r)+";"+str(x.i);}
};

template<unsigned D> typename MP<D>::M kineticGenerator(const SpinlessTvChainUtils&c){
    using X=MP<D>;using R=typename X::Real;using C=typename X::C;using M=typename X::M;
    int L=c.Lx;M h(2*L);C I(0,1),id(0,R(compact(c.delta))),im(0,R(compact(c.mu)));
    for(int x=0;x<L-1;++x){C z=x%2?X::neg(I):I;for(int k=0;k<2;++k){int a=k*L+x,b=a+1;h(a,b)=z;h(b,a)=X::neg(z);}int a=x,b=x+1;h(a,b)=X::add(h(a,b),id);h(b,a)=X::sub(h(b,a),id);a=L+x;b=a+1;h(a,b)=X::sub(h(a,b),id);h(b,a)=X::add(h(b,a),id);}
    if(c.boundaryType==0){int a=L-1,b=0,aa=2*L-1,bb=L;if(L%2==0){h(a,b)=X::add(h(a,b),X::sub(id,I));h(b,a)=X::sub(h(b,a),X::sub(id,I));h(aa,bb)=X::sub(h(aa,bb),X::add(id,I));h(bb,aa)=X::add(h(bb,aa),X::add(id,I));}else{h(a,bb)=X::sub(h(a,bb),X::add(I,id));h(bb,a)=X::add(h(bb,a),X::add(I,id));h(aa,b)=X::add(h(aa,b),X::sub(I,id));h(b,aa)=X::sub(h(b,aa),X::sub(I,id));}}
    for(int x=0;x<L;++x){int a=x,b=L+x;h(a,b)=X::sub(h(a,b),im);h(b,a)=X::add(h(b,a),im);}return h;
}

template<unsigned D> typename MP<D>::M hsFactor(const SpinlessTvChainUtils&cfg,const SpinlessVOperator&v,int flipAux){
    using X=MP<D>;using R=typename X::Real;using C=typename X::C;using M=typename X::M;
    using boost::multiprecision::acosh;using boost::multiprecision::cosh;using boost::multiprecision::exp;using boost::multiprecision::sinh;
    M b=M::eye(cfg.nDim);R lam=acosh(exp(R("0.5")*R(compact(v.localV))*R(compact(cfg.dt)))),ch=cosh(lam),sh=sinh(lam);
    for(int q=0;q<v.s->size();++q){int sigma=(*v.s)(q);if(q==flipAux)sigma=-sigma;int a,d,c,e;cfg.aux2MajoranaIdx(q,0,v.bondType,a,d);cfg.aux2MajoranaIdx(q,1,v.bondType,c,e);if(v.hsScheme==0){for(auto p:{std::pair<int,int>{a,d},std::pair<int,int>{c,e}}){b(p.first,p.first)=C(ch);b(p.second,p.second)=C(ch);b(p.first,p.second)=C(0,sh*R(sigma));b(p.second,p.first)=C(0,-sh*R(sigma));}}else{b(a,a)=C(ch);b(c,c)=C(ch);b(a,c)=C(0,sh*R(sigma));b(c,a)=C(0,-sh*R(sigma));b(d,d)=C(ch);b(e,e)=C(ch);b(d,e)=C(0,-sh*R(sigma));b(e,d)=C(0,sh*R(sigma));}}
    return b;
}

template<unsigned D> struct OracleResult{
    using X=MP<D>;typename X::C direct,sequential,q;typename X::Real orderProduct=0,orderGreen=0,orderRatio=0,r2q=0,structure=0,solveResidual=0,condition=0,symKinetic=0,symHS=0,symProduct=0;
};

template<unsigned D> OracleResult<D> evaluateOracle(const SpinlessTvChainUtils&cfg,const std::vector<Operator*>&ops,int trialSlices,int boundary,int aux){
    using X=MP<D>;using R=typename X::Real;using C=typename X::C;using M=typename X::M;
    M h=kineticGenerator<D>(cfg),full=h,half=h;X::scale(full,-R(compact(cfg.dt)));X::scale(half,-R(compact(cfg.dt))/R(2));full=X::expm(full);half=X::expm(half);
    auto factor=[&](int index,bool flipped){if(auto*v=dynamic_cast<SpinlessVOperator*>(ops[index]))return hsFactor<D>(cfg,*v,(flipped&&index==boundary)?aux:-1);if(dynamic_cast<DenseOperator*>(ops[index]))return index<trialSlices?full:half;throw std::runtime_error("unsupported contour factor");};
    auto forward=[&](bool flipped){M p=M::eye(cfg.nDim);for(int off=0;off<int(ops.size());++off){int i=(boundary+off)%int(ops.size());p=X::mulm(factor(i,flipped),p);}return p;};
    auto reverse=[&](bool flipped){M p=M::eye(cfg.nDim);for(int step=1;step<=int(ops.size());++step){int i=(boundary-step+int(ops.size()))%int(ops.size());p=X::mulm(p,factor(i,flipped));}return p;};
    M b=forward(false),bp=forward(true),bind=reverse(false),bpind=reverse(true);typename X::Solve sold,soln,soldind;M g=X::green(b,&sold),gp=X::green(bp,&soln),gind=X::green(bind,&soldind);X::green(bpind,nullptr);
    auto ratio=[&](M gin){auto*v=dynamic_cast<SpinlessVOperator*>(ops[boundary]);if(!v)throw std::runtime_error("target is not HS operator");int a,d,c,e;cfg.aux2MajoranaIdx(aux,0,v->bondType,a,d);cfg.aux2MajoranaIdx(aux,1,v->bondType,c,e);C I(0,1);R lam=boost::multiprecision::acosh(boost::multiprecision::exp(R("0.5")*R(compact(v->localV))*R(compact(cfg.dt)))),th=boost::multiprecision::tanh(lam),ch=boost::multiprecision::cosh(lam),eta=ch*ch,sigma=R((*v->s)(aux));C d1,d2,r;if(v->hsScheme==0){d1=X::sub(C(1),X::mulc(X::mulc(I,C(th*sigma)),gin(a,d)));d2=X::sub(C(1),X::mulc(X::mulc(I,C(th*sigma)),gin(c,e)));r=X::add(X::mulc(d1,d2),X::mulc(C(th*th),X::sub(X::mulc(gin(a,c),gin(d,e)),X::mulc(gin(d,c),gin(a,e)))));}else{d1=X::sub(C(1),X::mulc(X::mulc(I,C(th*sigma)),gin(a,c)));d2=X::add(C(1),X::mulc(X::mulc(I,C(th*sigma)),gin(d,e)));r=X::sub(X::mulc(d1,d2),X::mulc(C(th*th),X::sub(X::mulc(gin(a,d),gin(c,e)),X::mulc(gin(c,d),gin(a,e)))));}r=X::mulc(r,C(eta));
        int p1=v->hsScheme==0?a:a,q1=v->hsScheme==0?d:c,p2=v->hsScheme==0?c:d,q2=v->hsScheme==0?e:e;R second=v->hsScheme==0?R(1):R(-1);C r1=X::mulc(d1,C(ch));std::vector<C>x1(gin.n),x2(gin.n);for(int i=0;i<gin.n;++i){x1[i]=X::neg(gin(i,p1));x2[i]=X::neg(gin(i,q1));}x1[p1]=X::add(x1[p1],C(2));x2[q1]=X::add(x2[q1],C(2));C alpha=X::divc(X::mulc(I,C(sigma*th)),d1);for(int i=0;i<gin.n;++i)for(int j=0;j<gin.n;++j)gin(i,j)=X::add(gin(i,j),X::mulc(alpha,X::sub(X::mulc(x1[i],x2[j]),X::mulc(x2[i],x1[j]))));C d2actual=X::sub(C(1),X::mulc(X::mulc(I,C(second*th*sigma)),gin(p2,q2)));C seq=X::mulc(r1,X::mulc(d2actual,C(ch)));return std::pair<C,C>{r,seq};};
    auto rr=ratio(g),rri=ratio(gind);C q=X::divc(soln.determinant,sold.determinant);C r2=X::mulc(rr.first,rr.first);OracleResult<D> out;out.direct=rr.first;out.sequential=rr.second;out.q=q;out.orderProduct=std::max(X::rel(b,bind),X::rel(bp,bpind));out.orderGreen=X::rel(g,gind);out.orderRatio=X::absc(X::sub(rr.first,rri.first))/std::max(X::absc(rr.first),R(1e-100));out.r2q=X::absc(X::sub(r2,q))/std::max(X::absc(q),R(1e-100));out.structure=X::structure(g);out.solveResidual=std::max(sold.residual,soln.residual);out.condition=std::max(sold.condition,soln.condition);
    auto sym=[&](const M&m){M x(m.n);for(int i=0;i<m.n;++i){R di=(i%cfg.Lx)%2?R(-1):R(1);for(int j=0;j<m.n;++j){R dj=(j%cfg.Lx)%2?R(-1):R(1);x(i,j)=X::sub(X::mulc(C(di*dj),X::conjc(m(i,j))),m(i,j));}}R numerator=X::frob(x),denominator=std::max(X::frob(m),R("1e-100"));return R(numerator/denominator);};out.symKinetic=std::max(sym(full),sym(half));out.symHS=0;for(Operator*o:ops)if(auto*v=dynamic_cast<SpinlessVOperator*>(o))out.symHS=std::max(out.symHS,sym(hsFactor<D>(cfg,*v,-1)));out.symProduct=sym(b);return out;
}

const SpinlessTvChainUtils* currentCfg=nullptr;const std::vector<Operator*>*currentOps=nullptr;int currentTrial=0;

template<unsigned D> void emit(std::ostream&out,const std::string&label,int L,double V,int seed,int measurement,const std::string&sweep,int boundary,int aux,DataType doubleRatio,unsigned long long fieldHash){
    using X=MP<D>;auto z=evaluateOracle<D>(*currentCfg,*currentOps,currentTrial,boundary,aux);using boost::multiprecision::atan2;auto phase=atan2(z.direct.i,z.direct.r);out<<label<<','<<L<<','<<V<<','<<seed<<','<<measurement<<','<<sweep<<','<<boundary<<','<<aux<<','<<fieldHash<<','<<D<<','<<std::setprecision(17)<<doubleRatio.real()<<','<<doubleRatio.imag()<<','<<X::str(z.direct.r)<<','<<X::str(z.direct.i)<<','<<X::str(phase)<<','<<X::str(z.sequential.r)<<','<<X::str(z.sequential.i)<<','<<X::str(z.q.r)<<','<<X::str(z.q.i)<<','<<X::str(z.r2q)<<','<<X::str(z.orderProduct)<<','<<X::str(z.orderGreen)<<','<<X::str(z.orderRatio)<<','<<X::str(z.structure)<<','<<X::str(z.solveResidual)<<','<<X::str(z.condition)<<','<<X::str(z.symKinetic)<<','<<X::str(z.symHS)<<','<<X::str(z.symProduct)<<'\n';}

unsigned long long fieldHash(const std::vector<Operator*>&ops){unsigned long long h=1469598103934665603ULL;for(Operator*o:ops)if(auto*s=o->getAuxField())for(int i=0;i<s->size();++i){h^=static_cast<unsigned long long>(static_cast<long long>((*s)(i)));h*=1099511628211ULL;}return h;}

struct Captured{};
#ifndef PFQMC_MP_ORACLE_SMALL_ONLY
void runEvent(const std::string&label,int L,double V,int seed,int burn,int measurement,const std::string&sweep,int boundary,int aux,const std::string&trajectory,const std::string&path){
    constexpr double dt=.1,beta=8.;SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(L/dt)),0,1,0,0);rdGenerator rd(seed);Walker w(&cfg,&rd,L,beta);PfQMC q(&w,10);q.configureReadOnlyUpdateAudit(true);burnAdvance(q,burn,measurement);if(sweep=="left")q.rightSweep();std::ofstream out(path);out<<"label,L,V,seed,measurement,sweep,boundary,aux,field_hash,precision_digits,double_ratio_real,double_ratio_imag,rank4_real,rank4_imag,rank4_phase,sequential_real,sequential_imag,Q_real,Q_imag,relative_r2_minus_Q,product_order_relative,green_order_relative,ratio_order_relative,green_structure_residual,solve_residual,condition_estimate,kinetic_symmetry_residual,hs_symmetry_residual,product_symmetry_residual\n";
    if(trajectory=="full_reference"&&sweep=="right"){MatType fullStart;q.rebuildGreenFromFullContourAtBoundary(0,fullStart);q.g=fullStart;}
    bool captured=false;auto evaluate=[&](DataType ratio){if(captured)return;captured=true;currentCfg=&cfg;currentOps=&q.op_array;currentTrial=w.ntrial;auto h=fieldHash(q.op_array);emit<80>(out,label,L,V,seed,measurement,sweep,boundary,aux,ratio,h);emit<160>(out,label,L,V,seed,measurement,sweep,boundary,aux,ratio,h);emit<320>(out,label,L,V,seed,measurement,sweep,boundary,aux,ratio,h);throw Captured();};
    q.configureFlipReplayHook([&](const FlipReplayRecord&rec){if(rec.stage==0&&rec.sweep_direction==(sweep=="right"?0:1)&&rec.boundary==boundary&&rec.aux==aux)evaluate(rec.ratio);});q.configureGreenRecoveryEventHook([&](const GreenRecoveryEvent&e){if(e.source=="COMPLEX_RATIO"&&e.boundary==boundary&&e.aux==aux)evaluate(e.rebuilt_ratio);});
    if(trajectory=="full_reference")q.configureGreenRecoveryPrototype(true,std::numeric_limits<double>::infinity(),-1,-1,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());else if(trajectory=="prototype_old")q.configureGreenRecoveryPrototype(true,.2,1e-10,-1,.8,100);
    try{if(sweep=="right")q.rightSweep();else q.leftSweep();}catch(const Captured&){}if(!captured)throw std::runtime_error("target proposal not captured");
}
#endif

MatType gammaMatrix(int L,int which){int dim=1<<L,site=which%L;bool second=which>=L;MatType g=MatType::Zero(dim,dim);for(int state=0;state<dim;++state){int flipped=state^(1<<site);int parity=__builtin_popcount(unsigned(state&((1<<site)-1)))&1;double sign=parity?-1:1;if(!second)g(flipped,state)=sign;else g(flipped,state)=DataType(0,(state&(1<<site))?sign:-sign);}return g;}
MatType fockExpm(const MatType&a){MatType x=a;return expm(x,1.0);}
MatType fockFactor(const SpinlessTvChainUtils&cfg,Operator*op,int index,int trial,const std::vector<MatType>&gamma,int flipAux){MatType A=MatType::Zero(cfg.nDim,cfg.nDim);if(dynamic_cast<DenseOperator*>(op)){cfg.KineticGenerator(A);A*=-(index<trial?cfg.dt:cfg.dt/2);}else if(auto*v=dynamic_cast<SpinlessVOperator*>(op)){double lam=std::acosh(std::exp(.5*v->localV*cfg.dt));for(int q=0;q<v->s->size();++q){int sigma=(*v->s)(q);if(q==flipAux)sigma=-sigma;int a,b,c,d;cfg.aux2MajoranaIdx(q,0,v->bondType,a,b);cfg.aux2MajoranaIdx(q,1,v->bondType,c,d);if(v->hsScheme==0){A(a,b)=DataType(0,lam*sigma);A(b,a)=-A(a,b);A(c,d)=DataType(0,lam*sigma);A(d,c)=-A(c,d);}else{A(a,c)=DataType(0,lam*sigma);A(c,a)=-A(a,c);A(b,d)=DataType(0,-lam*sigma);A(d,b)=-A(b,d);}}}else throw std::runtime_error("unsupported Fock factor");MatType gen=MatType::Zero(gamma[0].rows(),gamma[0].cols());for(int i=0;i<cfg.nDim;++i)for(int j=i+1;j<cfg.nDim;++j)gen+=.5*A(i,j)*(gamma[i]*gamma[j]);return fockExpm(gen);}
void runSmall(int L,int seed,int samples,const std::string&path){constexpr double dt=.1,beta=8.;double V=5;SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(L/dt)),0,1,0,0);rdGenerator rd(seed);Walker w(&cfg,&rd,L,beta);std::vector<int>hs;for(int i=0;i<int(w.op_array.size());++i)if(dynamic_cast<SpinlessVOperator*>(w.op_array[i]))hs.push_back(i);std::vector<MatType>gamma;for(int i=0;i<2*L;++i)gamma.push_back(gammaMatrix(L,i));std::ofstream out(path);out<<"L,seed,sample,boundary,aux,field_hash,fock_ratio_real,fock_ratio_imag,mp_rank4_real,mp_rank4_imag,mp_sequential_real,mp_sequential_imag,Q_real,Q_imag,relative_fock_minus_rank4,relative_rank4_squared_minus_Q\n";currentCfg=&cfg;currentOps=&w.op_array;currentTrial=w.ntrial;for(int s=0;s<samples;++s){int boundary=hs[(size_t(s)*hs.size())/samples],aux=s%dynamic_cast<SpinlessVOperator*>(w.op_array[boundary])->s->size();MatType p=MatType::Identity(1<<L,1<<L),pn=p;for(int off=0;off<int(w.op_array.size());++off){int i=(boundary+off)%int(w.op_array.size());p=fockFactor(cfg,w.op_array[i],i,w.ntrial,gamma,-1)*p;pn=fockFactor(cfg,w.op_array[i],i,w.ntrial,gamma,(i==boundary)?aux:-1)*pn;}DataType fr=pn.trace()/p.trace();auto z=evaluateOracle<160>(cfg,w.op_array,w.ntrial,boundary,aux);using X=MP<160>;using R=typename X::Real;using C=typename X::C;C f(R(compact(fr.real())),R(compact(fr.imag())));auto ferr=X::absc(X::sub(f,z.direct))/std::max(X::absc(z.direct),R(1e-100));out<<L<<','<<seed<<','<<s<<','<<boundary<<','<<aux<<','<<fieldHash(w.op_array)<<','<<std::setprecision(17)<<fr.real()<<','<<fr.imag()<<','<<X::str(z.direct.r)<<','<<X::str(z.direct.i)<<','<<X::str(z.sequential.r)<<','<<X::str(z.sequential.i)<<','<<X::str(z.q.r)<<','<<X::str(z.q.i)<<','<<X::str(ferr)<<','<<X::str(z.r2q)<<'\n';}}

template<unsigned D> typename MP<D>::C mpPhase(typename MP<D>::C z){using X=MP<D>;auto a=X::absc(z);if(a==0)throw std::runtime_error("zero MP phase");return X::divc(z,typename X::C(a));}
template<unsigned D> typename MP<D>::C mpPfaffian(typename MP<D>::M a){using X=MP<D>;using C=typename X::C;using R=typename X::Real;C pf(1);for(int k=0;k<a.n;k+=2){int p=k+1;R best=X::absc(a(k,p));for(int j=k+2;j<a.n;++j)if(X::absc(a(k,j))>best){best=X::absc(a(k,j));p=j;}if(best==0)throw std::runtime_error("zero MP Pfaffian pivot");if(p!=k+1){for(int j=0;j<a.n;++j)std::swap(a(k+1,j),a(p,j));for(int i=0;i<a.n;++i)std::swap(a(i,k+1),a(i,p));pf=X::neg(pf);}C pivot=a(k,k+1);pf=X::mulc(pf,pivot);for(int i=k+2;i<a.n;++i)for(int j=i+1;j<a.n;++j){C correction=X::divc(X::sub(X::mulc(a(k,i),a(k+1,j)),X::mulc(a(k,j),a(k+1,i))),pivot);a(i,j)=X::sub(a(i,j),correction);a(j,i)=X::neg(a(i,j));}}return pf;}

template<unsigned D> typename MP<D>::C mpRawSign(const SpinlessTvChainUtils&cfg,const std::vector<Operator*>&ops,int trial,int boundary,int aux,bool flipped,int rotation=0){using X=MP<D>;using R=typename X::Real;using C=typename X::C;using M=typename X::M;M h=kineticGenerator<D>(cfg),full=h,half=h;X::scale(full,-R(compact(cfg.dt)));X::scale(half,-R(compact(cfg.dt))/R(2));full=X::expm(full);half=X::expm(half);auto factor=[&](int i){if(auto*v=dynamic_cast<SpinlessVOperator*>(ops[i]))return hsFactor<D>(cfg,*v,(flipped&&i==boundary)?aux:-1);if(dynamic_cast<DenseOperator*>(ops[i]))return i<trial?full:half;throw std::runtime_error("unsupported MP raw factor");};const int n=int(ops.size());rotation=(rotation%n+n)%n;auto at=[&](int step){return(rotation+step)%n;};M product=factor(at(0)),gcur=X::green(product);C sign(1);const C extra=((cfg.nDim/2)%2==0)?C(1):C(-1);for(int step=1;step<n;++step){int i=at(step);M gnext=X::green(factor(i)),block(2*cfg.nDim);for(int r=0;r<cfg.nDim;++r)for(int c=0;c<cfg.nDim;++c){block(r,c)=gnext(r,c);block(cfg.nDim+r,cfg.nDim+c)=gcur(r,c);}for(int r=0;r<cfg.nDim;++r){block(r,cfg.nDim+r)=C(-1);block(cfg.nDim+r,r)=C(1);}sign=X::mulc(sign,X::mulc(mpPhase<D>(mpPfaffian<D>(block)),extra));if(step+1<n){product=X::mulm(factor(i),product);gcur=X::green(product);}}return mpPhase<D>(sign);}

struct TargetCaptured{};
template<class F> void captureProposal(PfQMC&q,int half,int slice,int aux,const std::string&direction,F evaluate){int cycle=half/2;for(int i=0;i<cycle;++i){q.rightSweep();q.leftSweep();}if(direction=="left")q.rightSweep();MatType tmp=MatType::Identity(q.nDim,q.nDim),A=tmp;try{if(direction=="right"){int seg=0;for(int l=0;l<q.op_length;++l){Operator*op=q.op_array[l];if(l==slice){for(int a=0;a<=aux;++a){double u=0;op->prepareSingleFlip(q.g,a,&u);DataType ratio=op->preparedRatio();if(a==aux){evaluate(ratio,u,l);throw TargetCaptured();}bool accept=u<std::abs(ratio);q.sign*=op->finishSingleFlip(q.g,accept,true);}}else q.sign*=op->update(q.g);op->left_multiply(A,tmp);std::swap(A,tmp);if(q.need_stabilization[(l+1)%q.op_length]){if(seg==0)q.udtR[seg]=UDT(A);else q.udtR[seg]=A*q.udtR[seg-1];A.setIdentity();if(seg==q.checkpoints-1)q.udtR[seg].onePlusInv(q.g);else q.g=onePlusInv(q.udtL[seg+1],q.udtR[seg]);++seg;}else op->left_propagate(q.g,tmp);}}else{int seg=q.checkpoints-1;for(int l=q.op_length-1;l>=0;--l){Operator*op=q.op_array[l];op->right_propagate(q.g,tmp);if(l==slice){for(int a=0;a<=aux;++a){double u=0;op->prepareSingleFlip(q.g,a,&u);DataType ratio=op->preparedRatio();if(a==aux){evaluate(ratio,u,l);throw TargetCaptured();}bool accept=u<std::abs(ratio);q.sign*=op->finishSingleFlip(q.g,accept,true);}}else q.sign*=op->update(q.g);op->right_multiply(A,tmp);std::swap(A,tmp);if(q.need_stabilization[l]){A.adjointInPlace();if(seg==q.checkpoints-1)q.udtL[seg]=UDT(A);else q.udtL[seg]=A*q.udtL[seg+1];A.setIdentity();if(seg==0){q.udtL[seg].onePlusInv(q.g);q.g.adjointInPlace();}else q.g=onePlusInv(q.udtL[seg],q.udtR[seg-1]);--seg;}}}}catch(const TargetCaptured&){return;}throw std::runtime_error("target proposal not reached");}

template<unsigned D> void runTarget(int L,double theta,double V,double beta,double dt,int boundaryType,int hs,int seed,int half,int slice,int aux,const std::string&direction,bool dense,const std::string&path){SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(theta/dt)),boundaryType,1,0,hs);rdGenerator rd(seed);Walker w(&cfg,&rd,theta,beta);PfQMC q(&w,10);q.configureAdaptiveGuard(false,.1,100);std::ofstream out(path);out<<"precision,L,theta,V,boundary_type,hs_scheme,seed,half_step,direction,slice,aux,field_hash,double_fast_re,double_fast_im,uniform,mp_ratio_re,mp_ratio_im,mp_ratio_phase,mp_Q_re,mp_Q_im,mp_raw_pre_re,mp_raw_pre_im,mp_raw_post_re,mp_raw_post_im,dense_enabled,dense_weight_pre_re,dense_weight_pre_im,dense_weight_post_re,dense_weight_post_im,dense_ratio_re,dense_ratio_im,dense_sign_pre_re,dense_sign_pre_im,dense_sign_post_re,dense_sign_post_im,mp_condition,mp_residual,relative_r2_minus_Q\n";captureProposal(q,half,slice,aux,direction,[&](DataType fast,double u,int boundary){using X=MP<D>;auto z=evaluateOracle<D>(cfg,q.op_array,w.ntrial,boundary,aux);auto sr0=mpRawSign<D>(cfg,q.op_array,w.ntrial,boundary,aux,false),sr1=mpRawSign<D>(cfg,q.op_array,w.ntrial,boundary,aux,true);DataType wp(0),wn(0),fr(0),sp(0),sn(0);if(dense){std::vector<MatType>gamma;for(int i=0;i<2*L;++i)gamma.push_back(gammaMatrix(L,i));MatType p=MatType::Identity(1<<L,1<<L),pn=p;for(int off=0;off<int(q.op_array.size());++off){int i=(boundary+off)%int(q.op_array.size());p=fockFactor(cfg,q.op_array[i],i,w.ntrial,gamma,-1)*p;pn=fockFactor(cfg,q.op_array[i],i,w.ntrial,gamma,(i==boundary)?aux:-1)*pn;}wp=p.trace();wn=pn.trace();fr=wn/wp;sp=wp/std::abs(wp);sn=wn/std::abs(wn);}using boost::multiprecision::atan2;out<<D<<','<<L<<','<<theta<<','<<V<<','<<boundaryType<<','<<hs<<','<<seed<<','<<half<<','<<direction<<','<<boundary<<','<<aux<<','<<fieldHash(q.op_array)<<','<<std::setprecision(17)<<fast.real()<<','<<fast.imag()<<','<<u<<','<<X::str(z.direct.r)<<','<<X::str(z.direct.i)<<','<<X::str(atan2(z.direct.i,z.direct.r))<<','<<X::str(z.q.r)<<','<<X::str(z.q.i)<<','<<X::str(sr0.r)<<','<<X::str(sr0.i)<<','<<X::str(sr1.r)<<','<<X::str(sr1.i)<<','<<dense<<','<<wp.real()<<','<<wp.imag()<<','<<wn.real()<<','<<wn.imag()<<','<<fr.real()<<','<<fr.imag()<<','<<sp.real()<<','<<sp.imag()<<','<<sn.real()<<','<<sn.imag()<<','<<X::str(z.condition)<<','<<X::str(z.solveResidual)<<','<<X::str(z.r2q)<<'\n';});}

template<unsigned D> void runBoundary(int L,double theta,double V,double beta,double dt,int boundaryType,int hs,int seed,int half,int slice,int aux,const std::string&direction,const std::string&path){SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(theta/dt)),boundaryType,1,0,hs);rdGenerator rd(seed);Walker w(&cfg,&rd,theta,beta);PfQMC q(&w,10);std::ofstream out(path);out<<"precision,L,theta,V,boundary_type,hs_scheme,seed,half_step,direction,target_slice,target_aux,field_hash,rotation,rotation_role,pre_re,pre_im,post_re,post_im,ratio_re,ratio_im,pre_z2,post_z2\n";captureProposal(q,half,slice,aux,direction,[&](DataType,double,int boundary){using X=MP<D>;std::vector<std::pair<int,std::string>> rotations={{0,"canonical_zero"},{w.ntrial,"trial_physical"},{w.center,"measurement_center"},{boundary,"proposal_boundary"},{(boundary+1)%q.op_length,"proposal_after"},{q.op_length/4,"quarter"},{q.op_length/2,"half"},{3*q.op_length/4,"three_quarter"}};std::vector<int>seen;for(const auto&item:rotations){int rot=(item.first%q.op_length+q.op_length)%q.op_length;if(std::find(seen.begin(),seen.end(),rot)!=seen.end())continue;seen.push_back(rot);auto a=mpRawSign<D>(cfg,q.op_array,w.ntrial,boundary,aux,false,rot),b=mpRawSign<D>(cfg,q.op_array,w.ntrial,boundary,aux,true,rot),ratio=X::divc(b,a);out<<D<<','<<L<<','<<theta<<','<<V<<','<<boundaryType<<','<<hs<<','<<seed<<','<<half<<','<<direction<<','<<boundary<<','<<aux<<','<<fieldHash(q.op_array)<<','<<rot<<','<<item.second<<','<<X::str(a.r)<<','<<X::str(a.i)<<','<<X::str(b.r)<<','<<X::str(b.i)<<','<<X::str(ratio.r)<<','<<X::str(ratio.i)<<','<<(a.r<0?-1:1)<<','<<(b.r<0?-1:1)<<'\n';}});}

template<unsigned D> void runCheckpoint(int L,double theta,double V,double beta,double dt,int boundaryType,int hs,int seed,int half,unsigned long long expectedHash,const std::string&path){
    SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(theta/dt)),boundaryType,1,0,hs);rdGenerator rd(seed);Walker w(&cfg,&rd,theta,beta);
    PfQMC q(&w,10);q.configureAdaptiveGuard(false,.1,100);q.configureLeftSweepGreenRecovery(false);
    int cycle=half/2;bool right=(half%2)==0;for(int c=0;c<cycle;++c){q.rightSweep();q.leftSweep();}q.rightSweep();if(!right)q.leftSweep();
    auto hash=fieldHash(q.op_array);if(hash!=expectedHash)throw std::runtime_error("checkpoint field hash mismatch");
    using X=MP<D>;std::ofstream out(path);out<<"precision,L,theta,V,boundary_type,hs_scheme,seed,half_step,direction,field_hash,rotation,rotation_role,mp_raw_re,mp_raw_im,mp_raw_z2,expected_z2\n";
    std::vector<std::pair<int,std::string>> rotations={{0,"canonical_zero"},{w.ntrial,"trial_physical"},{w.center,"measurement_center"},{q.op_length/4,"quarter"},{q.op_length/2,"half"},{3*q.op_length/4,"three_quarter"}};std::vector<int>seen;
    for(const auto&item:rotations){int rot=(item.first%q.op_length+q.op_length)%q.op_length;if(std::find(seen.begin(),seen.end(),rot)!=seen.end())continue;seen.push_back(rot);auto a=mpRawSign<D>(cfg,q.op_array,w.ntrial,-1,-1,false,rot);out<<D<<','<<L<<','<<theta<<','<<V<<','<<boundaryType<<','<<hs<<','<<seed<<','<<half<<','<<(right?"right":"left")<<','<<hash<<','<<rot<<','<<item.second<<','<<X::str(a.r)<<','<<X::str(a.i)<<','<<(a.r<0?-1:1)<<",1\n";}
}

void runAdaptiveCheckpoint(int L,double theta,double V,double beta,double dt,int boundaryType,int hs,int seed,int half,unsigned long long expectedHash,const std::string&path){SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(theta/dt)),boundaryType,1,0,hs);rdGenerator rd(seed);Walker w(&cfg,&rd,theta,beta);PfQMC q(&w,10);int cycle=half/2;bool right=(half%2)==0;for(int c=0;c<cycle;++c){q.rightSweep();q.leftSweep();}q.rightSweep();if(!right)q.leftSweep();auto hash=fieldHash(q.op_array);if(hash!=expectedHash)throw std::runtime_error("adaptive checkpoint field hash mismatch");MpZ2Result r=projector_mp_z2::fullContourZ2(q.op_array);std::ofstream out(path);out<<"L,theta,V,boundary_type,hs_scheme,seed,half_step,direction,field_hash,z2,status,precision_digits,canonical_order,converged,reality_error,residual_or_condition,precision_escalations,message,expected_z2\n";out<<L<<','<<theta<<','<<V<<','<<boundaryType<<','<<hs<<','<<seed<<','<<half<<','<<(right?"right":"left")<<','<<hash<<','<<r.z2<<','<<mpZ2StatusName(r.status)<<','<<r.precision_digits<<','<<r.canonical_order<<','<<r.converged<<','<<std::setprecision(17)<<r.reality_error<<','<<r.residual_or_condition<<','<<r.precision_escalations<<",\""<<r.message<<"\",1\n";}

DataType denseContourPhase(const SpinlessTvChainUtils&cfg,const Walker&w){int L=cfg.Lx,dim=1<<L;std::vector<MatType>gamma;for(int i=0;i<2*L;++i)gamma.push_back(gammaMatrix(L,i));MatType p=MatType::Identity(dim,dim);for(int i=0;i<int(w.op_array.size());++i){p=fockFactor(cfg,w.op_array[i],i,w.ntrial,gamma,-1)*p;double n=p.norm();if(!std::isfinite(n)||n==0)throw std::runtime_error("dense Fock contour scaling failed");p/=n;}DataType z=p.trace();if(std::abs(z)==0)throw std::runtime_error("dense Fock zero trace");return z/std::abs(z);}
void runControl(int L,double theta,double V,double beta,double dt,int boundaryType,int hs,int seed,int samples,bool dense,const std::string&path){SpinlessTvChainUtils c0(L,dt,V,2*int(std::llround(theta/dt)),boundaryType,1,0,hs),c1(L,dt,V,2*int(std::llround(theta/dt)),boundaryType,1,0,hs);rdGenerator r0(seed),r1(seed);Walker w0(&c0,&r0,theta,beta),w1(&c1,&r1,theta,beta);PfQMC legacy(&w0,10);auto prod=[](const std::vector<Operator*>&ops){return projector_mp_z2::fullContourZ2(ops);};PfQMC real(&w1,10,PfQMCSignMode::real_z2,prod);std::ofstream out(path);out<<"L,theta,V,boundary_type,hs_scheme,seed,configuration,field_hash,dense_enabled,dense_phase_re,dense_phase_im,mp160_re,mp160_im,mp320_re,mp320_im,adaptive_z2,adaptive_status,adaptive_precision,transported_z2,legacy_z2,hs_rng_shadow\n";for(int s=0;s<samples;++s){DataType densePhase(0);if(dense)densePhase=denseContourPhase(c1,w1);auto a=mpRawSign<160>(c1,w1.op_array,w1.ntrial,-1,-1,false,0);auto b=mpRawSign<320>(c1,w1.op_array,w1.ntrial,-1,-1,false,0);MpZ2Result pr=prod(w1.op_array);bool shadow=fieldHash(w0.op_array)==fieldHash(w1.op_array)&&r0.diagnosticStateHash()==r1.diagnosticStateHash();out<<L<<','<<theta<<','<<V<<','<<boundaryType<<','<<hs<<','<<seed<<','<<s<<','<<fieldHash(w1.op_array)<<','<<dense<<','<<std::setprecision(17)<<densePhase.real()<<','<<densePhase.imag()<<','<<MP<160>::str(a.r)<<','<<MP<160>::str(a.i)<<','<<MP<320>::str(b.r)<<','<<MP<320>::str(b.i)<<','<<pr.z2<<','<<mpZ2StatusName(pr.status)<<','<<pr.precision_digits<<','<<real.physicalZ2Sign()<<','<<(legacy.sign.real()>=0?1:-1)<<','<<shadow<<'\n';legacy.rightSweep();real.rightSweep();legacy.leftSweep();real.leftSweep();}}

}

int main(int argc,char**argv)try{
    mkl_set_num_threads(1);if(argc<2)throw std::runtime_error("usage: oracle target|small ...");std::string mode=argv[1];
#ifndef PFQMC_MP_ORACLE_SMALL_ONLY
    if(mode=="event"){if(argc!=13)throw std::runtime_error("event label L V seed burn measurement right|left boundary aux trajectory out.csv");runEvent(argv[2],std::stoi(argv[3]),std::stod(argv[4]),std::stoi(argv[5]),std::stoi(argv[6]),std::stoi(argv[7]),argv[8],std::stoi(argv[9]),std::stoi(argv[10]),argv[11],argv[12]);}
    else
#endif
    if(mode=="small"){if(argc!=6)throw std::runtime_error("small L seed samples out.csv");runSmall(std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]),argv[5]);}
    else if(mode=="target"){if(argc!=17)throw std::runtime_error("target L theta V beta dt boundary hs seed half slice aux direction precision dense out.csv");int D=std::stoi(argv[14]);bool dense=std::stoi(argv[15]);if(D==80)runTarget<80>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),argv[13],dense,argv[16]);else if(D==160)runTarget<160>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),argv[13],dense,argv[16]);else if(D==320)runTarget<320>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),argv[13],dense,argv[16]);else throw std::runtime_error("precision must be 80,160,320");}
    else if(mode=="boundary"){if(argc!=16)throw std::runtime_error("boundary L theta V beta dt boundary hs seed half slice aux direction precision out.csv");int D=std::stoi(argv[14]);if(D==160)runBoundary<160>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),argv[13],argv[15]);else if(D==320)runBoundary<320>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),argv[13],argv[15]);else throw std::runtime_error("precision must be 160 or 320");}
    else if(mode=="checkpoint"){if(argc!=14)throw std::runtime_error("checkpoint L theta V beta dt boundary hs seed half precision expected_hash out.csv");int D=std::stoi(argv[11]);auto h=std::stoull(argv[12]);if(D==160)runCheckpoint<160>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),h,argv[13]);else if(D==320)runCheckpoint<320>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),h,argv[13]);else if(D==640)runCheckpoint<640>(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),h,argv[13]);else throw std::runtime_error("precision must be 160, 320 or 640");}
    else if(mode=="adaptive_checkpoint"){if(argc!=13)throw std::runtime_error("adaptive_checkpoint L theta V beta dt boundary hs seed half expected_hash out.csv");runAdaptiveCheckpoint(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),std::stoull(argv[11]),argv[12]);}
    else if(mode=="control"){if(argc!=13)throw std::runtime_error("control L theta V beta dt boundary hs seed samples dense out.csv");runControl(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10]),std::stoi(argv[11]),argv[12]);}
    else throw std::runtime_error("invalid mode");std::cout<<"{\"status\":\"complete\"}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
