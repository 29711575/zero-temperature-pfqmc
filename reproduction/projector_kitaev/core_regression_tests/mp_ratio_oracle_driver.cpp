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

unsigned long long fieldHash(const std::vector<Operator*>&ops){unsigned long long h=1469598103934665603ULL;for(Operator*o:ops)if(auto*s=o->getAuxField())for(int i=0;i<s->size();++i){h^=unsigned((*s)(i)+2);h*=1099511628211ULL;}return h;}

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

}

int main(int argc,char**argv)try{
    mkl_set_num_threads(1);if(argc<2)throw std::runtime_error("usage: mp_ratio_oracle event|small ...");std::string mode=argv[1];
#ifndef PFQMC_MP_ORACLE_SMALL_ONLY
    if(mode=="event"){if(argc!=13)throw std::runtime_error("event label L V seed burn measurement right|left boundary aux trajectory out.csv");runEvent(argv[2],std::stoi(argv[3]),std::stod(argv[4]),std::stoi(argv[5]),std::stoi(argv[6]),std::stoi(argv[7]),argv[8],std::stoi(argv[9]),std::stoi(argv[10]),argv[11],argv[12]);}
    else
#endif
    if(mode=="small"){if(argc!=6)throw std::runtime_error("small L seed samples out.csv");runSmall(std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]),argv[5]);}else throw std::runtime_error("invalid mode");std::cout<<"{\"status\":\"complete\"}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
