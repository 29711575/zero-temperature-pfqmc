#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "qr_udt.h"
#include "driven_comparator_reconstruction.h"

static double relativeError(const MatType &actual,const MatType &expected)
{
    return (actual-expected).norm()/std::max(expected.norm(),
        std::numeric_limits<double>::min());
}

int main(int argc,char **argv)
{
    if(argc!=2) {
        std::cerr<<"usage: dexp_comparator_regression output.csv\n";
        return 2;
    }

    UDT factor(2);
    factor.U=MatType::Identity(2,2);
    factor.T=MatType::Identity(2,2);
    factor.D(0)=0.5;
    factor.D(1)=0.5;
    factor.Dexp(0)=21;
    factor.Dexp(1)=-19;

    MatType expected=MatType::Zero(2,2);
    expected(0,0)=std::ldexp(0.5,21);
    expected(1,1)=std::ldexp(0.5,-19);
    const MatType fixedDense=drivenComparatorDenseProduct(factor);
    const MatType legacyDense=factor.U*factor.D.asDiagonal()*factor.T;

    MatType productionGreen;
    factor.onePlusInv(productionGreen);
    const MatType comparatorGreen=drivenComparatorGreenReconstruction(factor);
    const MatType directGreen=2.0*(MatType::Identity(2,2)+expected).inverse();

    const std::complex<double> comparatorLogdet=drivenComparatorLogdet(factor);
    const std::complex<double> directLogdet=logDet(MatType::Identity(2,2)+expected);

    const double legacyReconstructionError=relativeError(legacyDense,expected);
    const double fixedReconstructionError=relativeError(fixedDense,expected);
    const double comparatorVsOnePlusInv=relativeError(comparatorGreen,productionGreen);
    const double comparatorVsDirect=relativeError(comparatorGreen,directGreen);
    const double logdetError=std::abs(comparatorLogdet-directLogdet);
    const bool pass=fixedReconstructionError<=1e-15 &&
        comparatorVsOnePlusInv<=1e-13 && comparatorVsDirect<=1e-13 &&
        logdetError<=1e-13;

    std::ofstream out(argv[1]);
    out<<"case,D0,Dexp0,D1,Dexp1,legacy_reconstruction_relative_error,"
          "fixed_reconstruction_relative_error,comparator_vs_onePlusInv_relative_error,"
          "comparator_vs_direct_relative_error,logdet_abs_error,status\n";
    out<<std::setprecision(17)<<"synthetic_nonzero_dexp,"<<factor.D(0)<<','
       <<factor.Dexp(0)<<','<<factor.D(1)<<','<<factor.Dexp(1)<<','
       <<legacyReconstructionError<<','<<fixedReconstructionError<<','
       <<comparatorVsOnePlusInv<<','<<comparatorVsDirect<<','<<logdetError<<','
       <<(pass?"PASS":"FAIL")<<'\n';
    std::cout<<"synthetic_nonzero_dexp="<<(pass?"PASS":"FAIL")<<'\n';
    return pass?0:1;
}
