#include "test_common.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using MPReal=boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1800>>;
struct MPC { MPReal r=0,i=0; MPC(){} MPC(MPReal x):r(x){} MPC(MPReal x,MPReal y):r(x),i(y){} };
struct MPM {
    int n=0; std::vector<MPC> a;
    explicit MPM(int d=0):n(d),a(std::size_t(d)*d){}
    MPC& operator()(int i,int j){return a[std::size_t(i)*n+j];}
    const MPC& operator()(int i,int j)const{return a[std::size_t(i)*n+j];}
    static MPM eye(int n){MPM x(n);for(int i=0;i<n;++i)x(i,i)=MPC(1);return x;}
};
MPC mpAdd(MPC a,const MPC&b){a.r+=b.r;a.i+=b.i;return a;}
MPC mpSub(MPC a,const MPC&b){a.r-=b.r;a.i-=b.i;return a;}
MPC mpMul(const MPC&a,const MPC&b){return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};}
MPC mpConj(MPC a){a.i=-a.i;return a;}
MPReal mpAbs2(const MPC&a){return a.r*a.r+a.i*a.i;}
MPReal mpAbs(const MPC&a){using boost::multiprecision::sqrt;return sqrt(mpAbs2(a));}
MPC mpScale(MPC a,const MPReal&s){a.r*=s;a.i*=s;return a;}
MPM mpAddM(const MPM&a,const MPM&b){MPM c(a.n);for(std::size_t k=0;k<c.a.size();++k)c.a[k]=mpAdd(a.a[k],b.a[k]);return c;}
MPM mpSubM(const MPM&a,const MPM&b){MPM c(a.n);for(std::size_t k=0;k<c.a.size();++k)c.a[k]=mpSub(a.a[k],b.a[k]);return c;}
MPM mpAdj(const MPM&a){MPM c(a.n);for(int i=0;i<a.n;++i)for(int j=0;j<a.n;++j)c(i,j)=mpConj(a(j,i));return c;}
MPM mpMulM(const MPM&a,const MPM&b){MPM c(a.n);for(int i=0;i<a.n;++i)for(int k=0;k<a.n;++k){const MPC z=a(i,k);for(int j=0;j<a.n;++j)c(i,j)=mpAdd(c(i,j),mpMul(z,b(k,j)));}return c;}
MPReal mpNorm(const MPM&a){using boost::multiprecision::sqrt;MPReal x=0;for(const MPC&z:a.a)x+=mpAbs2(z);return sqrt(x);}
MPReal mpRel(const MPM&a,const MPM&b){return mpNorm(mpSubM(a,b))/std::max(mpNorm(b),MPReal("1e-150"));}
MPM mpFromDouble(const MatType&a){MPM x(a.rows());for(int i=0;i<a.rows();++i)for(int j=0;j<a.cols();++j)x(i,j)=MPC(MPReal(a(i,j).real()),MPReal(a(i,j).imag()));return x;}
MPReal mpPow2(int e){using boost::multiprecision::pow;return pow(MPReal(2),e);}
int mpExponent2(const MPReal&z){using boost::multiprecision::floor;using boost::multiprecision::log;return floor(log(z)/log(MPReal(2))).convert_to<int>()+1;}
std::string mpString(const MPReal&x){std::ostringstream s;s<<std::setprecision(32)<<std::scientific<<x;return s.str();}

MPM mpReconstruct(const UDT&u){
    MPM d(u.nDim);for(int i=0;i<u.nDim;++i)d(i,i)=MPC(MPReal(u.D(i))*mpPow2(u.Dexp(i)));
    return mpMulM(mpMulM(mpFromDouble(u.U),d),mpFromDouble(u.T));
}
MPReal mpBackwardResidual(const MPM&A,const MatType&g){
    MPM system=mpAddM(MPM::eye(A.n),A),gm=mpFromDouble(g),lhs=mpMulM(system,gm),identity=MPM::eye(A.n);
    for(MPC&z:lhs.a)z=mpScale(z,MPReal("0.5"));
    return mpNorm(mpSubM(lhs,identity))/(mpNorm(system)*mpNorm(gm)*MPReal("0.5")+mpNorm(identity));
}
struct MPSolve {MPM inverse;MPReal residual=0,condition=0;};
MPSolve mpInverse(const MPM&input){
    const int n=input.n;MPM a=input,b=MPM::eye(n),original=input;
    for(int k=0;k<n;++k){int pivot=k;MPReal best=mpAbs(a(k,k));for(int i=k+1;i<n;++i)if(mpAbs(a(i,k))>best){best=mpAbs(a(i,k));pivot=i;}if(best==0)throw std::runtime_error("MP singular matrix");if(pivot!=k)for(int j=0;j<n;++j){std::swap(a(k,j),a(pivot,j));std::swap(b(k,j),b(pivot,j));}MPC diag=a(k,k);MPReal den=mpAbs2(diag);MPC invdiag(diag.r/den,-diag.i/den);for(int j=0;j<n;++j){a(k,j)=mpMul(a(k,j),invdiag);b(k,j)=mpMul(b(k,j),invdiag);}for(int i=0;i<n;++i)if(i!=k){MPC f=a(i,k);if(mpAbs(f)==0)continue;for(int j=0;j<n;++j){a(i,j)=mpSub(a(i,j),mpMul(f,a(k,j)));b(i,j)=mpSub(b(i,j),mpMul(f,b(k,j)));}}}
    MPSolve s;s.inverse=b;MPM identity=MPM::eye(n);s.residual=mpNorm(mpSubM(mpMulM(original,b),identity))/(mpNorm(original)*mpNorm(b)+mpNorm(identity));s.condition=mpNorm(original)*mpNorm(b);return s;
}

MatType randomUnitary(int n,std::mt19937&rng){
    std::normal_distribution<double>d(0,1);MatType a(n,n);
    for(int j=0;j<n;++j)for(int i=0;i<n;++i)a(i,j)=DataType(d(rng),d(rng));
    Eigen::HouseholderQR<MatType>qr(a);return (qr.householderQ()*MatType::Identity(n,n)).eval();
}
MatType randomDense(int n,std::mt19937&rng){
    std::normal_distribution<double>d(0,.15/std::sqrt(double(n)));MatType a=MatType::Identity(n,n);
    for(int j=0;j<n;++j)for(int i=0;i<n;++i)a(i,j)+=DataType(d(rng),d(rng));return a;
}
void assignScales(UDT&u,int range,std::mt19937&rng,bool cancelling,const UDT*other){
    std::uniform_real_distribution<double>m(.5,.999999);
    for(int i=0;i<u.nDim;++i){int e=-range+int((2LL*range*i)/std::max(1,u.nDim-1));if(cancelling&&other)e=-other->Dexp(i)+(i%3-1);u.D(i)=m(rng);u.Dexp(i)=e;}
}
UDT makeLegacy(int n,int range,std::mt19937&rng,bool cancelling,const UDT*other=nullptr){
    UDT u(n);u.U=randomUnitary(n,rng);u.T=randomUnitary(n,rng);assignScales(u,range,rng,cancelling,other);return u;
}
UDT makeCanonical(int n,int range,std::mt19937&rng,bool cancelling,const UDT*other=nullptr){
    MatType a=randomDense(n,rng);UDT u(a);assignScales(u,range,rng,cancelling,other);return u;
}

struct Metric {double adj=0,right=0,trans=0,maxOff=0,minCol=0,maxCol=0,minSv=0,maxSv=0;bool finite=false;};
Metric metric(const MatType&a){
    Metric m;m.finite=scaleSafeFinite(a);const int rows=a.rows(),cols=a.cols();
    MatType gc=a.adjoint()*a-MatType::Identity(cols,cols);m.adj=gc.norm()/std::sqrt(double(cols));
    MatType gr=a*a.adjoint()-MatType::Identity(rows,rows);m.right=gr.norm()/std::sqrt(double(rows));
    if(rows==cols)m.trans=(a.transpose()*a-MatType::Identity(cols,cols)).norm()/std::sqrt(double(cols));else m.trans=std::numeric_limits<double>::quiet_NaN();
    m.maxOff=0;for(int i=0;i<cols;++i)for(int j=0;j<cols;++j)if(i!=j)m.maxOff=std::max(m.maxOff,std::abs((a.col(i).dot(a.col(j)))));
    m.minCol=std::numeric_limits<double>::infinity();m.maxCol=0;for(int j=0;j<cols;++j){double z=a.col(j).norm();m.minCol=std::min(m.minCol,z);m.maxCol=std::max(m.maxCol,z);}
    Eigen::JacobiSVD<MatType>svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);m.minSv=svd.singularValues().minCoeff();m.maxSv=svd.singularValues().maxCoeff();return m;
}
struct TContract {int uniqueTerminal=0;double maxPivotMagnitudeError=0;double maxBelow=0;};
TContract tContract(const MatType&t){
    const int n=t.rows();std::set<int> terminal;TContract c;
    for(int j=0;j<n;++j){int p=-1;for(int i=0;i<n;++i)if(std::abs(t(i,j))>1e-13)p=i;if(p>=0){terminal.insert(p);c.maxPivotMagnitudeError=std::max(c.maxPivotMagnitudeError,std::abs(std::abs(t(p,j))-1.0));for(int i=p+1;i<n;++i)c.maxBelow=std::max(c.maxBelow,std::abs(t(i,j)));}}
    c.uniqueTerminal=terminal.size();return c;
}

void writeStage(std::ofstream&out,const std::string&generator,int n,int range,int seed,const std::string&stage,int k,int pivot,int exponent,double selectedNorm,const Metric&m,double p1=0,double p2=0,double minRemain=0,double maxRemain=0,const std::string&recon="nan"){
    out<<generator<<','<<n<<','<<range<<','<<seed<<','<<stage<<','<<k<<','<<pivot<<','<<exponent<<','<<std::setprecision(17)<<selectedNorm<<','<<m.adj<<','<<m.right<<','<<m.trans<<','<<m.maxOff<<','<<m.minCol<<','<<m.maxCol<<','<<m.minSv<<','<<m.maxSv<<','<<p1<<','<<p2<<','<<minRemain<<','<<maxRemain<<','<<(m.finite?1:0)<<','<<recon<<'\n';
}

struct Prepared {MatType work; iVecType exponent;};
Prepared prepareBi(const UDT&left,const MatType&x,const UDT&right){
    const int n=left.nDim;Prepared p{MatType::Zero(n,n),iVecType(n)};
    for(int j=0;j<n;++j){int ce=std::numeric_limits<int>::min();for(int i=0;i<n;++i){DataType z=x(i,j)*left.D(i)*right.D(j);if(std::abs(z)>0){int e=0;std::frexp(std::abs(z),&e);ce=std::max(ce,left.Dexp(i)+right.Dexp(j)+e);}}p.exponent(j)=ce;for(int i=0;i<n;++i)p.work(i,j)=scaleSafeLdexp(x(i,j)*left.D(i)*right.D(j),left.Dexp(i)+right.Dexp(j)-ce);}
    return p;
}
struct Instrumented {UDT local;std::vector<int>pivotPositions;std::vector<int>pivotOriginal;Prepared prepared;};
Instrumented instrument(const Prepared&input,std::ofstream&stages,const std::string&generator,int n,int range,int seed){
    MatType work=input.work,q=MatType::Zero(n,n),t=MatType::Zero(n,n);iVecType exponent=input.exponent;std::vector<int>perm(n);for(int j=0;j<n;++j)perm[j]=j;
    auto renorm=[&](int j){double z=work.col(j).stableNorm();int e=0;std::frexp(z,&e);work.col(j)*=scaleSafePow2(-e);exponent(j)+=e;return z;};for(int j=0;j<n;++j)renorm(j);
    Instrumented result;result.prepared=input;result.local.nDim=n;result.local.D.resize(n);result.local.Dexp.resize(n);
    for(int k=0;k<n;++k){
        int pivot=k;for(int j=k+1;j<n;++j)if(exponent(j)>exponent(pivot)||(exponent(j)==exponent(pivot)&&work.col(j).squaredNorm()>work.col(pivot).squaredNorm()))pivot=j;
        result.pivotPositions.push_back(pivot);if(pivot!=k){work.col(k).swap(work.col(pivot));std::swap(exponent(k),exponent(pivot));std::swap(perm[k],perm[pivot]);}result.pivotOriginal.push_back(perm[k]);
        double norm=work.col(k).stableNorm();int de=0;double dm=std::frexp(norm,&de);result.local.D(k)=dm;result.local.Dexp(k)=exponent(k)+de;q.col(k)=work.col(k)/norm;t(k,perm[k])=1.0;
        double max1=0,max2=0,minRemain=std::numeric_limits<double>::infinity(),maxRemain=0;
        for(int j=k+1;j<n;++j){DataType c1=q.col(k).dot(work.col(j));work.col(j)-=q.col(k)*c1;max1=std::max(max1,std::abs(q.col(k).dot(work.col(j))));DataType c2=q.col(k).dot(work.col(j));work.col(j)-=q.col(k)*c2;max2=std::max(max2,std::abs(q.col(k).dot(work.col(j))));int shift=exponent(j)-result.local.Dexp(k);t(k,perm[j])=scaleSafeLdexp((c1+c2)/dm,shift);double rn=work.col(j).stableNorm();minRemain=std::min(minRemain,rn);maxRemain=std::max(maxRemain,rn);renorm(j);}
        MatType partial=q.leftCols(k+1);Metric pm=metric(partial);writeStage(stages,generator,n,range,seed,"mgs_step",k,pivot,result.local.Dexp(k),norm,pm,max1,max2,k+1<n?minRemain:0,k+1<n?maxRemain:0);
    }
    result.local.U=q;result.local.T=t;return result;
}

struct MPRef {MatType q;MatType t;dVecType d;iVecType dexp;MPReal orth;MPReal recon;double minOverlap=1,maxDdiff=0,tRelative=0;};
MPRef mpReference(const Prepared&input,const Instrumented&prod){
    const int n=input.work.rows();std::vector<std::vector<MPC>>v(n,std::vector<MPC>(n)),q(n,std::vector<MPC>(n));std::vector<int>exp(n),perm(n);for(int j=0;j<n;++j){exp[j]=input.exponent(j);perm[j]=j;for(int i=0;i<n;++i)v[j][i]=MPC(MPReal(input.work(i,j).real()),MPReal(input.work(i,j).imag()));}
    auto dot=[&](const std::vector<MPC>&a,const std::vector<MPC>&b){MPC z;for(int i=0;i<n;++i)z=mpAdd(z,mpMul(mpConj(a[i]),b[i]));return z;};
    auto norm=[&](const std::vector<MPC>&a){MPReal z=0;for(const MPC&x:a)z+=mpAbs2(x);using boost::multiprecision::sqrt;return sqrt(z);};
    auto renorm=[&](int j){MPReal z=norm(v[j]);int e=mpExponent2(z);MPReal scale=mpPow2(-e);for(MPC&x:v[j])x=mpScale(x,scale);exp[j]+=e;};for(int j=0;j<n;++j)renorm(j);
    std::vector<std::vector<MPC>>t(n,std::vector<MPC>(n));std::vector<MPReal>d(n);std::vector<int>de(n);
    for(int k=0;k<n;++k){int pivot=prod.pivotPositions[k];if(pivot!=k){std::swap(v[k],v[pivot]);std::swap(exp[k],exp[pivot]);std::swap(perm[k],perm[pivot]);}
        for(int pass=0;pass<2;++pass)for(int h=0;h<k;++h){MPC c=dot(q[h],v[k]);for(int i=0;i<n;++i)v[k][i]=mpSub(v[k][i],mpMul(q[h][i],c));}
        MPReal z=norm(v[k]);int de0=mpExponent2(z);MPReal dm=z/mpPow2(de0);d[k]=dm;de[k]=exp[k]+de0;for(int i=0;i<n;++i)q[k][i]=mpScale(v[k][i],MPReal(1)/z);t[k][perm[k]]=MPC(1);
        for(int j=k+1;j<n;++j){MPC c1=dot(q[k],v[j]);for(int i=0;i<n;++i)v[j][i]=mpSub(v[j][i],mpMul(q[k][i],c1));MPC c2=dot(q[k],v[j]);for(int i=0;i<n;++i)v[j][i]=mpSub(v[j][i],mpMul(q[k][i],c2));t[k][perm[j]]=mpScale(mpAdd(c1,c2),mpPow2(exp[j]-de[k])/d[k]);renorm(j);}
    }
    MPRef r;r.q=MatType(n,n);r.t=MatType(n,n);r.d.resize(n);r.dexp.resize(n);for(int k=0;k<n;++k){r.d(k)=d[k].convert_to<double>();r.dexp(k)=de[k];for(int i=0;i<n;++i)r.q(i,k)=DataType(q[k][i].r.convert_to<double>(),q[k][i].i.convert_to<double>());for(int j=0;j<n;++j)r.t(k,j)=DataType(t[k][j].r.convert_to<double>(),t[k][j].i.convert_to<double>());}
    MPReal os=0;for(int a=0;a<n;++a)for(int b=0;b<n;++b){MPC z=dot(q[a],q[b]);if(a==b)z.r-=1;os+=mpAbs2(z);}using boost::multiprecision::sqrt;r.orth=sqrt(os/MPReal(n));
    UDT u(n);u.U=r.q;u.T=r.t;u.D=r.d;u.Dexp=r.dexp;MPM target=mpFromDouble(input.work);for(int j=0;j<n;++j){MPReal z=mpPow2(input.exponent(j));for(int i=0;i<n;++i)target(i,j)=mpScale(target(i,j),z);}r.recon=mpRel(mpReconstruct(u),target);
    r.minOverlap=1;for(int k=0;k<n;++k)r.minOverlap=std::min(r.minOverlap,std::abs(prod.local.U.col(k).dot(r.q.col(k))));r.maxDdiff=0;for(int k=0;k<n;++k)r.maxDdiff=std::max(r.maxDdiff,double(std::abs(prod.local.Dexp(k)-r.dexp(k))));r.tRelative=(prod.local.T-r.t).norm()/std::max(r.t.norm(),1e-300);return r;
}

double greenRel(const MatType&a,const MatType&b){return (a-b).norm()/std::max(b.norm(),1e-300);}
void process(const std::string&generator,int n,int range,int seed,std::ofstream&stages,std::ofstream&refs){
    std::mt19937 rng(seed);UDT right,left;if(generator=="legacy_invalid_T"){right=makeLegacy(n,range,rng,false);left=makeLegacy(n,range,rng,true,&right);}else{right=makeCanonical(n,range,rng,false);left=makeCanonical(n,range,rng,true,&right);}MatType bridge=left.T*right.U;Prepared prepared=prepareBi(left,bridge,right);
    Metric ml=metric(left.U),mr=metric(right.U),mb=metric(bridge),mw=metric(prepared.work);writeStage(stages,generator,n,range,seed,"input_left_U",-1,-1,0,0,ml);writeStage(stages,generator,n,range,seed,"input_right_U",-1,-1,0,0,mr);writeStage(stages,generator,n,range,seed,"bridge",-1,-1,0,0,mb);writeStage(stages,generator,n,range,seed,"prepared_work",-1,-1,0,0,mw);
    Instrumented ins=instrument(prepared,stages,generator,n,range,seed);Metric mq=metric(ins.local.U);writeStage(stages,generator,n,range,seed,"local_q",n,-1,0,0,mq);
    UDT product=ins.local;product.U=left.U*product.U;product.T=product.T*right.T;Metric mo=metric(product.U);MPM target=mpMulM(mpReconstruct(left),mpReconstruct(right));MPReal prodRecon=mpRel(mpReconstruct(product),target);writeStage(stages,generator,n,range,seed,"output_U",n,-1,0,0,mo,0,0,0,0,mpString(prodRecon));
    MPRef ref=mpReference(prepared,ins);UDT refProduct(n);refProduct.U=left.U*ref.q;refProduct.D=ref.d;refProduct.Dexp=ref.dexp;refProduct.T=ref.t*right.T;MPReal refRecon=mpRel(mpReconstruct(refProduct),target);
    MatType gp,gt;bool fp=true,ft=true;try{product.onePlusInv(gp);}catch(...){fp=false;}try{gt=onePlusInv(left,right);}catch(...){ft=false;}
    MPReal rp=fp?mpBackwardResidual(target,gp):MPReal("nan"),oracleResidual=MPReal("nan"),oracleCondition=MPReal("nan"),gdiff=MPReal("nan");if(n==52&&range==2000&&seed==903520){MPSolve oracle=mpInverse(mpAddM(MPM::eye(n),target));oracleResidual=oracle.residual;oracleCondition=oracle.condition;MPM gmp=oracle.inverse;for(MPC&z:gmp.a)z=mpScale(z,MPReal(2));if(fp)gdiff=mpRel(mpFromDouble(gp),gmp);}
    MPM cyclic=mpMulM(mpReconstruct(right),mpAdj(mpReconstruct(left)));MPReal rt=ft?mpBackwardResidual(cyclic,gt):MPReal("nan");TContract lc=tContract(left.T),rc=tContract(right.T);
    refs<<generator<<','<<n<<','<<range<<','<<seed<<','<<ml.adj<<','<<ml.right<<','<<ml.trans<<','<<mr.adj<<','<<mr.right<<','<<mr.trans<<','<<lc.uniqueTerminal<<','<<rc.uniqueTerminal<<','<<lc.maxPivotMagnitudeError<<','<<rc.maxPivotMagnitudeError<<','<<mq.adj<<','<<mo.adj<<','<<mpString(ref.orth)<<','<<ref.minOverlap<<','<<mpString(prodRecon)<<','<<mpString(refRecon)<<','<<mpString(ref.recon)<<','<<ref.maxDdiff<<','<<ref.tRelative<<','<<mpString(rp)<<','<<mpString(oracleResidual)<<','<<mpString(oracleCondition)<<','<<mpString(gdiff)<<','<<mpString(rt)<<','<<(fp&&ft?1:0)<<'\n';
}
}

int main(int argc,char**argv)try{
    if(argc!=3&&argc!=4)throw std::runtime_error("usage: udt_orthogonality_driver stages.csv reference.csv [smoke]");
    std::ofstream stages(argv[1]),refs(argv[2]);stages<<"generator,n,exponent_range,seed,stage,k,pivot_position,selected_exponent,selected_norm,adjoint_left_residual,adjoint_right_residual,transpose_residual,max_offdiag_inner_product,min_column_norm,max_column_norm,min_singular_value,max_singular_value,max_projection_after_first,max_projection_after_second,min_remaining_norm,max_remaining_norm,finite,reconstruction_residual\n";
    refs<<"generator,n,exponent_range,seed,left_adjoint_residual,left_right_residual,left_transpose_residual,right_adjoint_residual,right_right_residual,right_transpose_residual,left_T_unique_terminal_pivots,right_T_unique_terminal_pivots,left_T_pivot_magnitude_error,right_T_pivot_magnitude_error,production_q_adjoint_residual,output_U_adjoint_residual,mp_q_adjoint_residual,min_production_mp_column_overlap,production_reconstruction_residual,mp_reference_reconstruction_residual,mp_local_reconstruction_residual,max_D_exponent_difference,T_relative_difference,production_onePlusInv_backward_residual,mp_oracle_solve_residual,mp_oracle_condition,production_mp_green_relative,two_udt_backward_residual,finite\n";
    const bool smoke=argc==4;std::vector<std::tuple<int,int,int>>cases;
    if(smoke)cases.push_back(std::make_tuple(12,500,900012));
    else {cases.push_back(std::make_tuple(52,1500,903520));cases.push_back(std::make_tuple(52,2000,903517));cases.push_back(std::make_tuple(52,2000,903518));cases.push_back(std::make_tuple(52,2000,903519));cases.push_back(std::make_tuple(52,2000,903520));cases.push_back(std::make_tuple(52,2000,903521));cases.push_back(std::make_tuple(52,2000,903522));cases.push_back(std::make_tuple(24,2000,903520));}
    int count=0;for(const auto&c:cases){int n,r,s;std::tie(n,r,s)=c;if(smoke||n==52&&r==2000&&s==903520){process("legacy_invalid_T",n,r,s,stages,refs);++count;}process("canonical_T",n,r,s,stages,refs);++count;}std::cout<<"{\"status\":\"complete\",\"cases\":"<<count<<"}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
