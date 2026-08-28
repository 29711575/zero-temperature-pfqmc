#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"

namespace {
class Walker : public Spinless_tV {
public:
    int center=-1, ntrial=0, nphysical=0;
    Walker(const SpinlessTvChainUtils *config, rdGenerator *random,
           double theta, double beta) {
        build_projector_static_contour(*this,config,random,theta,beta,
                                       center,ntrial,nphysical);
    }
};

double oldConditionProxy(const PfQMC &q, int boundaryIndex) {
    int boundary=boundaryIndex%q.op_length;
    if(boundary<0) boundary+=q.op_length;
    UDT product(q.nDim);
    MatType block=MatType::Identity(q.nDim,q.nDim),tmp;
    for(int offset=0;offset<q.op_length;++offset){
        q.op_array[(boundary+offset)%q.op_length]->left_multiply(block,tmp);
        block.swap(tmp);product=block*product;block.setIdentity();
    }
    MatType lhs=product.T.inverse();
    dVecType dplus(q.nDim),dminus(q.nDim);
    for(int i=0;i<q.nDim;++i){
        dplus(i)=product.dLargeInverse(i);dminus(i)=product.dSmallPart(i);
    }
    MatType core=lhs*dplus.asDiagonal()+product.U*dminus.asDiagonal();
    Eigen::JacobiSVD<MatType> svd(core);
    const auto s=svd.singularValues();
    return s(0)/std::max(s(s.size()-1),std::numeric_limits<double>::min());
}
}

int main() try {
    MatType regular=MatType::Zero(2,2);
    regular(0,1)=DataType(2.0,0.0);regular(1,0)=-regular(0,1);
    MatType regularCopy=regular;
    const PfaffianResult regularResult=signOfPfafWithStatus(regularCopy);
    if(!regularResult.ok()||std::abs(regularResult.value-DataType(1,0))>1e-14)
        throw std::runtime_error("regular Pfaffian status failed");

    MatType zero=MatType::Zero(2,2);
    const PfaffianResult zeroResult=signOfPfafWithStatus(zero);
    if(zeroResult.ok()) throw std::runtime_error("zero Pfaffian did not fail closed");
    bool wrapperThrew=false;
    try { MatType z=MatType::Zero(2,2); (void)signOfPfaf(z); }
    catch(const std::exception&) { wrapperThrew=true; }
    if(!wrapperThrew) throw std::runtime_error("compatibility wrapper did not fail closed");

    SpinlessTvChainUtils config(6,.1,2,40,0,1,0,0);
    rdGenerator random(910001);Walker walker(&config,&random,2,2);PfQMC q(&walker,10);
    const double oldCondition=oldConditionProxy(q,walker.center);
    const double newCondition=q.fullContourCoreConditionAtBoundary(walker.center);
    const double conditionRelative=std::abs(oldCondition-newCondition)/
        std::max(std::abs(oldCondition),std::numeric_limits<double>::min());
    if(!std::isfinite(newCondition)||conditionRelative>1e-10)
        throw std::runtime_error("condition proxy solve does not match inverse reference");

    resetQRDiagonalPhaseDiagnostics();
    for(int sweep=0;sweep<5;++sweep){q.rightSweep();q.leftSweep();}
    const QRDiagonalPhaseDiagnostics qmcAudit=qrDiagonalPhaseDiagnostics();

    resetQRDiagonalPhaseDiagnostics();
    std::mt19937_64 generator(910002);
    std::normal_distribution<double> normal(0.0,1.0);
    for(int sample=0;sample<100;++sample){
        MatType matrix(12,12);
        for(int j=0;j<12;++j)for(int i=0;i<12;++i)
            matrix(i,j)=DataType(normal(generator),normal(generator));
        UDT factor(matrix);(void)factor;
    }
    const QRDiagonalPhaseDiagnostics randomAudit=qrDiagonalPhaseDiagnostics();
    std::cout<<std::setprecision(17)
      <<"{\"status\":\"complete\",\"pfaffian_regular_status\":\""
      <<pfaffianStatusName(regularResult.status)
      <<"\",\"pfaffian_zero_status\":\""<<pfaffianStatusName(zeroResult.status)
      <<"\",\"pfaffian_zero_lapack_info\":"<<zeroResult.lapack_info
      <<",\"pfaffian_zero_min_pivot\":";
    if(std::isfinite(zeroResult.min_pivot))std::cout<<zeroResult.min_pivot;else std::cout<<"null";
    std::cout<<",\"compatibility_wrapper_threw\":"<<(wrapperThrew?"true":"false")
      <<",\"condition_proxy_old\":"<<oldCondition
      <<",\"condition_proxy_new\":"<<newCondition
      <<",\"condition_proxy_relative_difference\":"<<conditionRelative
      <<",\"qmc_qr_diagonal_samples\":"<<qmcAudit.samples
      <<",\"qmc_qr_diagonal_max_imag_ratio\":"<<qmcAudit.maximum
      <<",\"qmc_qr_diagonal_above_1e12\":"<<qmcAudit.above_1e12
      <<",\"random_qr_diagonal_samples\":"<<randomAudit.samples
      <<",\"random_qr_diagonal_max_imag_ratio\":"<<randomAudit.maximum
      <<",\"random_qr_diagonal_above_1e12\":"<<randomAudit.above_1e12
      <<"}\n";
    return 0;
} catch(const std::exception &error) {
    std::cerr<<error.what()<<'\n';return 2;
}
