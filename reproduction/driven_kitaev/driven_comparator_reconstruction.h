#ifndef DRIVEN_COMPARATOR_RECONSTRUCTION_H
#define DRIVEN_COMPARATOR_RECONSTRUCTION_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

// Validation-only reconstruction helpers.  UDT::D is a mantissa; every
// diagnostic reconstruction must also account for UDT::Dexp.
struct DrivenComparatorSystem {
    MatType left;
    MatType core;
};

inline DrivenComparatorSystem drivenComparatorSystem(const UDT &factor)
{
    const MatType identity=MatType::Identity(factor.nDim,factor.nDim);
    const MatType inverseT=scaleSafeCheckedSolve(
        factor.T,identity,"driven comparator/T");
    dVecType largeInverse(factor.nDim),smallPart(factor.nDim);
    for(int i=0;i<factor.nDim;++i) {
        largeInverse(i)=factor.dLargeInverse(i);
        smallPart(i)=factor.dSmallPart(i);
    }
    DrivenComparatorSystem system;
    system.left=inverseT*largeInverse.asDiagonal();
    system.core=system.left+factor.U*smallPart.asDiagonal();
    return system;
}

inline MatType drivenComparatorGreenReconstruction(const UDT &factor)
{
    const DrivenComparatorSystem system=drivenComparatorSystem(factor);
    return 2.0*scaleSafeCheckedRightSolve(
        system.left,system.core,"driven comparator/core");
}

inline MatType drivenComparatorDenseProduct(const UDT &factor)
{
    dVecType scale(factor.nDim);
    for(int i=0;i<factor.nDim;++i) scale(i)=factor.actualD(i);
    return factor.U*scale.asDiagonal()*factor.T;
}

inline std::complex<double> drivenComparatorLogdet(const UDT &factor)
{
    const DrivenComparatorSystem system=drivenComparatorSystem(factor);
    Eigen::PartialPivLU<MatType> lu(system.core);
    const auto diagonal=lu.matrixLU().diagonal();
    std::complex<double> value=lu.permutationP().determinant()<0
        ? std::complex<double>(0.0,M_PI) : std::complex<double>(0.0,0.0);
    for(int i=0;i<diagonal.size();++i) {
        value+=std::log(diagonal(i));
        if(factor.dGreaterThanOne(i)) value+=std::log(factor.actualD(i));
    }
    return value;
}

inline double drivenComparatorScaleSpread(const UDT &factor)
{
    double minimum=std::numeric_limits<double>::infinity();
    double maximum=0.0;
    for(int i=0;i<factor.nDim;++i) {
        const double scale=factor.actualD(i);
        minimum=std::min(minimum,scale);
        maximum=std::max(maximum,scale);
    }
    return maximum/std::max(minimum,std::numeric_limits<double>::min());
}

#endif
