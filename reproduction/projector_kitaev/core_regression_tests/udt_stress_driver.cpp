#include "test_common.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>
#include <random>

namespace {
using Real=boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>>;
struct C { Real r=0,i=0; C(){} C(Real x):r(x){} C(Real x,Real y):r(x),i(y){} };
struct M { int n=0; std::vector<C>a; explicit M(int d=0):n(d),a(std::size_t(d)*d){} C&operator()(int i,int j){return a[std::size_t(i)*n+j];}const C&operator()(int i,int j)const{return a[std::size_t(i)*n+j];}static M eye(int n){M x(n);for(int i=0;i<n;++i)x(i,i)=C(1);return x;} };
C add(C a,const C&b){a.r+=b.r;a.i+=b.i;return a;} C sub(C a,const C&b){a.r-=b.r;a.i-=b.i;return a;} C mul(const C&a,const C&b){return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};} C conj(C a){a.i=-a.i;return a;}
Real abs2(const C&a){return a.r*a.r+a.i*a.i;} Real absv(const C&a){using boost::multiprecision::sqrt;return sqrt(abs2(a));}
M mmul(const M&a,const M&b){M c(a.n);for(int i=0;i<a.n;++i)for(int k=0;k<a.n;++k){C z=a(i,k);for(int j=0;j<a.n;++j)c(i,j)=add(c(i,j),mul(z,b(k,j)));}return c;}
M madd(const M&a,const M&b){M c(a.n);for(std::size_t k=0;k<c.a.size();++k)c.a[k]=add(a.a[k],b.a[k]);return c;}
M msub(const M&a,const M&b){M c(a.n);for(std::size_t k=0;k<c.a.size();++k)c.a[k]=sub(a.a[k],b.a[k]);return c;}
M adj(const M&a){M c(a.n);for(int i=0;i<a.n;++i)for(int j=0;j<a.n;++j)c(i,j)=conj(a(j,i));return c;}
Real norm(const M&a){using boost::multiprecision::sqrt;Real x=0;for(const C&z:a.a)x+=abs2(z);return sqrt(x);} Real rel(const M&a,const M&b){return norm(msub(a,b))/std::max(norm(b),Real("1e-90"));}
M fromDouble(const MatType&a){M x(a.rows());for(int i=0;i<a.rows();++i)for(int j=0;j<a.cols();++j)x(i,j)=C(Real(a(i,j).real()),Real(a(i,j).imag()));return x;}
Real pow2(int e){using boost::multiprecision::pow;return pow(Real(2),e);}
M reconstruct(const UDT&u){M U=fromDouble(u.U),T=fromDouble(u.T),d(u.nDim);for(int i=0;i<u.nDim;++i)d(i,i)=C(Real(u.D(i))*pow2(u.Dexp(i)));return mmul(mmul(U,d),T);}
std::string s(const Real&x){std::ostringstream o;o<<std::setprecision(24)<<std::scientific<<x;return o.str();}
MatType randomUnitary(int n,std::mt19937&rng){std::normal_distribution<double>d(0,1);MatType a(n,n);for(int j=0;j<n;++j)for(int i=0;i<n;++i)a(i,j)=DataType(d(rng),d(rng));Eigen::HouseholderQR<MatType>qr(a);return qr.householderQ()*MatType::Identity(n,n);}
MatType randomDense(int n,std::mt19937&rng){std::normal_distribution<double>d(0,.15/std::sqrt(double(n)));MatType a=MatType::Identity(n,n);for(int j=0;j<n;++j)for(int i=0;i<n;++i)a(i,j)+=DataType(d(rng),d(rng));return a;}
UDT makeUDT(int n,int range,std::mt19937&rng,bool cancelling,const UDT*other=nullptr){UDT u(n);u.U=randomUnitary(n,rng);u.T=randomUnitary(n,rng);std::uniform_real_distribution<double>m(.5,.999999);for(int i=0;i<n;++i){int e=-range+int((2LL*range*i)/std::max(1,n-1));if(cancelling&&other)e=-other->Dexp(i)+(i%3-1);u.D(i)=m(rng);u.Dexp(i)=e;}return u;}
double orth(const MatType&q){return (q.adjoint()*q-MatType::Identity(q.rows(),q.cols())).norm()/std::sqrt(double(q.rows()));}
double singleCoreRcond(const UDT&u){int n=u.nDim;MatType I=MatType::Identity(n,n),x=scaleSafeCheckedSolve(u.T,I,"stress/T");dVecType dp(n),dm(n);for(int i=0;i<n;++i){dp(i)=u.dLargeInverse(i);dm(i)=u.dSmallPart(i);}MatType core=x*dp.asDiagonal()+u.U*dm.asDiagonal();return core.fullPivLu().rcond();}
Real solveEquationResidual(const M&A,const MatType&g){M lhs=mmul(madd(M::eye(A.n),A),fromDouble(g));for(C&z:lhs.a){z.r*=Real("0.5");z.i*=Real("0.5");}return rel(lhs,M::eye(A.n));}
Real solveResidual(const M&A,const MatType&g){
 M system=madd(M::eye(A.n),A),gm=fromDouble(g),lhs=mmul(system,gm);
 for(C&z:lhs.a){z.r*=Real("0.5");z.i*=Real("0.5");}
 const M identity=M::eye(A.n);const Real denominator=norm(system)*norm(gm)*Real("0.5")+norm(identity);
 return norm(msub(lhs,identity))/std::max(denominator,Real("1e-90"));
}
}

int main(int argc,char**argv)try{
 if(argc!=5)throw std::runtime_error("usage: udt_stress_driver size exponent_range seed output.csv");int n=std::stoi(argv[1]),range=std::stoi(argv[2]),seed=std::stoi(argv[3]);std::mt19937 rng(seed);
 UDT right=makeUDT(n,range,rng,false),left=makeUDT(n,range,rng,true,&right);MatType B=randomDense(n,rng);
 M R=reconstruct(right),L=reconstruct(left),BM=fromDouble(B);
 UDT matProduct=B*right,udtProduct=left*right;
 Real matRecon=rel(reconstruct(matProduct),mmul(BM,R)),udtRecon=rel(reconstruct(udtProduct),mmul(L,R));
 MatType g1;right.onePlusInv(g1);Real oneResidual=solveResidual(R,g1),oneEquationResidual=solveEquationResidual(R,g1);
 MatType g2=onePlusInv(left,right);M cyclic=mmul(R,adj(L));Real twoResidual=solveResidual(cyclic,g2),twoEquationResidual=solveEquationResidual(cyclic,g2);
 MatType work=randomDense(n,rng);iVecType ex(n);for(int i=0;i<n;++i)ex(i)=-range+int((2LL*range*i)/std::max(1,n-1));UDT qr=UDT::factorPreparedScaled(work,ex);M scaled=fromDouble(work);for(int j=0;j<n;++j){Real z=pow2(ex(j));for(int i=0;i<n;++i){scaled(i,j).r*=z;scaled(i,j).i*=z;}}Real qrRecon=rel(reconstruct(qr),scaled);
 bool finite=scaleSafeFinite(matProduct.U)&&scaleSafeFinite(matProduct.T)&&scaleSafeFinite(udtProduct.U)&&scaleSafeFinite(udtProduct.T)&&matrixFinite(g1)&&matrixFinite(g2);
 int emin=qr.Dexp.minCoeff(),emax=qr.Dexp.maxCoeff();double rcond=singleCoreRcond(right);
 std::ofstream out(argv[4]);out<<"size,exponent_range,seed,mat_udt_reconstruction_error,udt_udt_reconstruction_error,onePlusInv_solve_residual,udtRdivudtL_solve_residual,onePlusInv_unscaled_equation_residual,udtRdivudtL_unscaled_equation_residual,exponent_qr_reconstruction_error,mat_udt_orthogonality,udt_udt_orthogonality,qr_orthogonality,finite,output_exponent_min,output_exponent_max,scaled_core_rcond,cancellation_case\n";
 out<<n<<','<<range<<','<<seed<<','<<s(matRecon)<<','<<s(udtRecon)<<','<<s(oneResidual)<<','<<s(twoResidual)<<','<<s(oneEquationResidual)<<','<<s(twoEquationResidual)<<','<<s(qrRecon)<<','<<std::setprecision(17)<<orth(matProduct.U)<<','<<orth(udtProduct.U)<<','<<orth(qr.U)<<','<<(finite?1:0)<<','<<emin<<','<<emax<<','<<rcond<<",1\n";
 std::cout<<"{\"status\":\"complete\",\"finite\":"<<(finite?"true":"false")<<",\"max_reconstruction\":\""<<s(std::max({matRecon,udtRecon,qrRecon}))<<"\",\"max_solve_residual\":\""<<s(std::max(oneResidual,twoResidual))<<"\"}\n";return finite?0:3;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
