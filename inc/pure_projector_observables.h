#ifndef PURE_PROJECTOR_OBSERVABLES_H
#define PURE_PROJECTOR_OBSERVABLES_H

#include "skewMatUtils.h"
#include "types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

inline std::pair<int, int> pureProjectorCheckerboardBondCounts(
    int length, int boundaryType) {
    if (boundaryType == 0) return {length / 2, length / 2};
    return {length / 2, (length - 1) / 2};
}

inline DataType pureProjectorStructureFactorUnnormalized(
    const MatType &green, int length, double momentum) {
    DataType sum = 0.0;
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            DataType correlation;
            if (i == j) {
                correlation = 0.25;
            } else {
                const int ai = i, bi = length + i;
                const int aj = j, bj = length + j;
                correlation = -0.25 *
                    (green(ai, bi) * green(aj, bj) -
                     green(ai, aj) * green(bi, bj) +
                     green(ai, bj) * green(bi, aj));
            }
            sum += std::exp(DataType(0.0, momentum * (i - j))) * correlation;
        }
    }
    return sum;
}

// Paper convention: S(q)=(1/L^2) sum_ij exp(iq(i-j)) <rho_i rho_j>.
inline DataType pureProjectorStructureFactor(
    const MatType &green, int length, double momentum) {
    if (length <= 0 || green.rows() != 2 * length || green.cols() != 2 * length)
        return DataType(std::numeric_limits<double>::quiet_NaN(), 0.0);
    return pureProjectorStructureFactorUnnormalized(green, length, momentum) /
           double(length * length);
}

struct PurePhysicalParityResult {
    PfaffianStatus pfaffian_status = PfaffianStatus::invalid_dimension;
    DataType internal_pfaffian = DataType(0.0, 0.0);
    int internal_pfaffian_sign = 0;
    int block_reordering_sign = 0;
    int physical_parity = 0;
    double skew_residual = std::numeric_limits<double>::infinity();
    double diagonal_residual = std::numeric_limits<double>::infinity();
    double reality_error = std::numeric_limits<double>::infinity();
    bool ok() const {
        return pfaffian_status == PfaffianStatus::success &&
               internal_pfaffian_sign != 0 && block_reordering_sign != 0 &&
               physical_parity != 0;
    }
};

// Coordinates are block ordered (gamma_0...gamma_{L-1},
// gamma'_0...gamma'_{L-1}).  sign Pf(-i G) is the internal Pfaffian
// convention.  The physical Fock parity additionally contains the permutation
// from block to site-interleaved Majoranas:
//   (-i)^L (-1)^[L(L-1)/2] Pf(G).
inline PurePhysicalParityResult pureProjectorPhysicalParity(
    const MatType &green, double tolerance = 1e-8) {
    PurePhysicalParityResult result;
    if (green.rows() <= 0 || green.rows() != green.cols() ||
        (green.rows() % 2) != 0) return result;
    for (int col = 0; col < green.cols(); ++col)
        for (int row = 0; row < green.rows(); ++row)
            if (!std::isfinite(green(row, col).real()) ||
                !std::isfinite(green(row, col).imag())) return result;

    const int length = green.rows() / 2;
    const double scale = std::max(1.0, std::sqrt(double(green.rows())));
    result.skew_residual = (green + green.transpose()).norm() / scale;
    result.diagonal_residual = green.diagonal().norm() / scale;
    if (result.skew_residual > tolerance || result.diagonal_residual > tolerance)
        return result;

    MatType antisymmetric = 0.5 * (green - green.transpose());
    MatType internalMatrix = DataType(0.0, -1.0) * antisymmetric;
    const PfaffianResult pfaffian = signOfPfafWithStatus(internalMatrix);
    result.pfaffian_status = pfaffian.status;
    result.internal_pfaffian = pfaffian.value;
    result.reality_error = std::abs(pfaffian.value.imag()) /
        std::max(1.0, std::abs(pfaffian.value.real()));
    if (!pfaffian.ok() || result.reality_error > tolerance ||
        std::abs(pfaffian.value.real()) < 1.0 - tolerance) return result;
    result.internal_pfaffian_sign = pfaffian.value.real() >= 0.0 ? 1 : -1;
    result.block_reordering_sign =
        ((length * (length - 1) / 2) & 1) ? -1 : 1;
    result.physical_parity =
        result.block_reordering_sign * result.internal_pfaffian_sign;
    return result;
}

#endif
