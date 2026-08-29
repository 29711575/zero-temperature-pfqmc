#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "gaussian_trial_state.h"

namespace {
constexpr double kPi=3.141592653589793238462643383279502884;

MatType gammaMatrix(int modes,int which){int dimension=1<<modes,site=which%modes;bool second=which>=modes;
    MatType gamma=MatType::Zero(dimension,dimension);for(int state=0;state<dimension;++state){int flipped=state^(1<<site);
        int odd=__builtin_popcount(unsigned(state&((1<<site)-1)))&1;double sign=odd?-1:1;
        gamma(flipped,state)=second?DataType(0,(state&(1<<site))?sign:-sign):DataType(sign,0);}return gamma;}

std::vector<MatType> gammas(int L){std::vector<MatType> result;for(int i=0;i<2*L;++i)result.push_back(gammaMatrix(L,i));return result;}

MatType kinetic(int L,int boundary,double t,double delta,double mu){auto one=[&](double d,double m){
    SpinlessTvChainUtils c(L,1,0,2,boundary,d,m,0);MatType h=MatType::Zero(2*L,2*L);c.KineticGenerator(h);return h;};
    MatType base=one(0,0);return t*base+delta*(one(1,0)-base)+mu*(one(0,1)-base);}

MatType denseQuadratic(const MatType&h,const std::vector<MatType>&gamma){MatType result=MatType::Zero(gamma[0].rows(),gamma[0].cols());
    for(int i=0;i<h.rows();++i)for(int j=i+1;j<h.cols();++j)result+=.5*h(i,j)*gamma[i]*gamma[j];return result;}

DataType structure(const MatType&g,int L,double q){DataType sum=0;for(int i=0;i<L;++i)for(int j=0;j<L;++j){
    DataType corr;if(i==j)corr=.25;else{int ai=i,bi=L+i,aj=j,bj=L+j;
        corr=-.25*(g(ai,bi)*g(aj,bj)-g(ai,aj)*g(bi,bj)+g(ai,bj)*g(bi,aj));}
    sum+=std::exp(DataType(0,q*(i-j)))*corr;}return sum/double(L);}

std::map<std::string,std::string> arguments(int argc,char**argv){std::map<std::string,std::string> out;
    for(int i=1;i<argc;++i){std::string key=argv[i];if(key.rfind("--",0)!=0||i+1>=argc)
        throw std::invalid_argument("arguments must be pairs");out[key.substr(2)]=argv[++i];}return out;}

} // namespace

int main(int argc,char**argv){try{auto a=arguments(argc,argv);auto get=[&](const char*k){auto it=a.find(k);
    if(it==a.end())throw std::invalid_argument(std::string("missing --")+k);return it->second;};
    int L=std::stoi(get("L"));double V=std::stod(get("V")),t=std::stod(get("t")),delta=std::stod(get("delta")),mu=std::stod(get("mu"));
    int boundary=get("boundary")=="pbc"?0:1;double tt=std::stod(get("trial-t")),td=std::stod(get("trial-delta")),tm=std::stod(get("trial-mu"));
    int policy=std::stoi(get("trial-parity"));double split=std::stod(get("edge-splitting"));if(L>10)throw std::invalid_argument("ED limited to L<=10");
    std::vector<MatType> gamma=gammas(L);MatType hTrial=kinetic(L,boundary,tt,td,tm);
    if(split!=0){SpinlessTvChainUtils coordinates(L,1,V,2,boundary,td,tm,0);int left=coordinates.majoranaCoord2Idx(0,1),right=coordinates.majoranaCoord2Idx(L-1,1);
        DataType z(0,split);hTrial(left,right)+=z;hTrial(right,left)-=z;}
    GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(hTrial);if(trial.fermionParity()!=policy)
        throw std::invalid_argument("trial parity mismatch");
    MatType denseTrial=denseQuadratic(hTrial,gamma);Eigen::SelfAdjointEigenSolver<MatType> ts(denseTrial);
    cVecType trialVector=ts.eigenvectors().col(0);double denseTrialParity=0;for(int state=0;state<(1<<L);++state)
        denseTrialParity+=(__builtin_popcount(unsigned(state))&1?-1:1)*std::norm(trialVector(state));
    int fockSector=denseTrialParity>=0?1:-1;

    MatType h=kinetic(L,boundary,t,delta,mu);MatType physical=denseQuadratic(h,gamma);
    std::vector<MatType> density(L);for(int i=0;i<L;++i)density[i]=DataType(0,-.5)*gamma[i]*gamma[L+i];
    int bonds=boundary==0?L:L-1;for(int i=0;i<bonds;++i)physical+=V*density[i]*density[(i+1)%L];
    std::vector<int> basis;for(int state=0;state<(1<<L);++state)if(((__builtin_popcount(unsigned(state))&1)?-1:1)==fockSector)basis.push_back(state);
    MatType sector(basis.size(),basis.size());for(int i=0;i<int(basis.size());++i)for(int j=0;j<int(basis.size());++j)sector(i,j)=physical(basis[i],basis[j]);
    Eigen::SelfAdjointEigenSolver<MatType> solver(sector);if(solver.info()!=Eigen::Success)throw std::runtime_error("ED eigensolver failed");
    cVecType state=cVecType::Zero(1<<L);for(int i=0;i<int(basis.size());++i)state(basis[i])=solver.eigenvectors()(i,0);
    MatType green=MatType::Zero(2*L,2*L);for(int i=0;i<2*L;++i)for(int j=0;j<2*L;++j)if(i!=j)
        green(i,j)=-(state.adjoint()*gamma[i]*gamma[j]*state)(0);
    DataType spi=structure(green,L,kPi),sdq=structure(green,L,kPi-2*kPi/L);
    std::cout<<std::setprecision(17)<<"{\"status\":\"complete\",\"L\":"<<L<<",\"V\":"<<V
        <<",\"boundary\":\""<<(boundary?"obc":"pbc")<<"\",\"trial_parity\":"<<policy
        <<",\"dense_fock_parity_sector\":"<<fockSector<<",\"energy\":"<<solver.eigenvalues()(0)
        <<",\"fermion_parity\":"<<policy<<",\"S_pi\":"<<spi.real()<<",\"S_pi_dq\":"<<sdq.real()
        <<",\"R_CDW\":"<<(1-sdq/spi).real()<<"}\n";return 0;
}catch(const std::exception&e){std::cerr<<"pure_projector_ed: "<<e.what()<<'\n';return 1;}}
