#ifndef PURE_PROJECTOR_GREEN_H
#define PURE_PROJECTOR_GREEN_H

#include "operator.h"
#include "types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

enum class PureProjectorStatus {
    success,
    invalid_dimension,
    nonfinite_input,
    subspace_rank_deficient,
    overlap_rank_deficient,
    overlap_ill_conditioned,
    solve_residual_exceeded,
    green_residual_exceeded,
    green_structure_exceeded
};

inline const char *pureProjectorStatusName(PureProjectorStatus status) {
    switch (status) {
    case PureProjectorStatus::success: return "success";
    case PureProjectorStatus::invalid_dimension: return "invalid_dimension";
    case PureProjectorStatus::nonfinite_input: return "nonfinite_input";
    case PureProjectorStatus::subspace_rank_deficient: return "subspace_rank_deficient";
    case PureProjectorStatus::overlap_rank_deficient: return "overlap_rank_deficient";
    case PureProjectorStatus::overlap_ill_conditioned: return "overlap_ill_conditioned";
    case PureProjectorStatus::solve_residual_exceeded: return "solve_residual_exceeded";
    case PureProjectorStatus::green_residual_exceeded: return "green_residual_exceeded";
    case PureProjectorStatus::green_structure_exceeded: return "green_structure_exceeded";
    }
    return "unknown";
}

struct PureProjectorOptions {
    double rank_tolerance = 1e-12;
    double minimum_overlap_rcond = 1e-12;
    double solve_residual_tolerance = 1e-10;
    double green_residual_tolerance = 1e-10;
    // Always compute skew/diagonal diagnostics.  Enabling this flag promotes
    // them to an additional fail-closed gate; it is off by default so this
    // review does not silently change the frozen Phase 3C trust policy.
    bool fail_on_green_structure = false;
};

struct ThinQrResult {
    PureProjectorStatus status = PureProjectorStatus::invalid_dimension;
    MatType q;
    int rank = 0;
    double orthonormal_residual = std::numeric_limits<double>::infinity();
    bool ok() const { return status == PureProjectorStatus::success; }
};

struct PureProjectorGreenResult {
    PureProjectorStatus status = PureProjectorStatus::invalid_dimension;
    MatType green;
    int overlap_rank = 0;
    double overlap_rcond = 0.0;
    double solve_residual = std::numeric_limits<double>::infinity();
    double green_residual = std::numeric_limits<double>::infinity();
    double green_skew_residual = std::numeric_limits<double>::infinity();
    double green_diagonal_residual = std::numeric_limits<double>::infinity();
    bool ok() const { return status == PureProjectorStatus::success; }
};

inline bool pureProjectorMatrixFinite(const MatType &matrix) {
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = 0; row < matrix.rows(); ++row) {
            if (!std::isfinite(matrix(row, col).real()) ||
                !std::isfinite(matrix(row, col).imag())) return false;
        }
    }
    return true;
}

inline ThinQrResult thinQrSubspace(
    const MatType &subspace,
    const PureProjectorOptions &options = PureProjectorOptions()) {
    ThinQrResult result;
    if (subspace.rows() <= 0 || subspace.cols() <= 0 ||
        subspace.rows() < subspace.cols()) return result;
    if (!pureProjectorMatrixFinite(subspace)) {
        result.status = PureProjectorStatus::nonfinite_input;
        return result;
    }
    Eigen::ColPivHouseholderQR<MatType> qr(subspace);
    qr.setThreshold(options.rank_tolerance);
    result.rank = qr.rank();
    if (result.rank != subspace.cols()) {
        result.status = PureProjectorStatus::subspace_rank_deficient;
        return result;
    }
    result.q = qr.householderQ() *
        MatType::Identity(subspace.rows(), subspace.cols());
    result.orthonormal_residual =
        (result.q.adjoint() * result.q -
         MatType::Identity(subspace.cols(), subspace.cols())).norm() /
        std::max(1.0, std::sqrt(double(subspace.cols())));
    if (!pureProjectorMatrixFinite(result.q) ||
        result.orthonormal_residual > options.solve_residual_tolerance) {
        result.status = PureProjectorStatus::solve_residual_exceeded;
        return result;
    }
    result.status = PureProjectorStatus::success;
    return result;
}

inline PureProjectorGreenResult pureProjectorGreen(
    const MatType &phiRight, const MatType &phiLeft,
    const PureProjectorOptions &options = PureProjectorOptions()) {
    PureProjectorGreenResult result;
    if (phiRight.rows() <= 0 || phiRight.cols() <= 0 ||
        phiRight.rows() != phiLeft.rows() || phiRight.cols() != phiLeft.cols() ||
        phiRight.rows() < phiRight.cols()) return result;
    if (!pureProjectorMatrixFinite(phiRight) || !pureProjectorMatrixFinite(phiLeft)) {
        result.status = PureProjectorStatus::nonfinite_input;
        return result;
    }

    const MatType overlap = phiLeft.adjoint() * phiRight;
    Eigen::JacobiSVD<MatType> svd(
        overlap, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(options.rank_tolerance);
    result.overlap_rank = svd.rank();
    if (result.overlap_rank != overlap.rows()) {
        result.status = PureProjectorStatus::overlap_rank_deficient;
        return result;
    }
    const double largest = svd.singularValues()(0);
    const double smallest = svd.singularValues()(svd.singularValues().size() - 1);
    result.overlap_rcond = largest > 0.0 ? smallest / largest : 0.0;
    if (!std::isfinite(result.overlap_rcond) ||
        result.overlap_rcond < options.minimum_overlap_rcond) {
        result.status = PureProjectorStatus::overlap_ill_conditioned;
        return result;
    }

    const MatType rightHandSide = phiLeft.adjoint();
    const MatType solved = svd.solve(rightHandSide);
    result.solve_residual = (overlap * solved - rightHandSide).norm() /
        std::max(1.0, rightHandSide.norm());
    if (!pureProjectorMatrixFinite(solved) ||
        result.solve_residual > options.solve_residual_tolerance) {
        result.status = PureProjectorStatus::solve_residual_exceeded;
        return result;
    }
    result.green = MatType::Identity(phiRight.rows(), phiRight.rows()) -
                   2.0 * phiRight * solved;
    result.green_residual =
        (result.green * result.green -
         MatType::Identity(result.green.rows(), result.green.cols())).norm() /
        std::max(1.0, std::sqrt(double(result.green.rows())));
    if (!pureProjectorMatrixFinite(result.green) ||
        result.green_residual > options.green_residual_tolerance) {
        result.status = PureProjectorStatus::green_residual_exceeded;
        return result;
    }
    const double greenScale =
        std::max(1.0, std::sqrt(double(result.green.rows())));
    result.green_skew_residual =
        (result.green + result.green.transpose()).norm() / greenScale;
    result.green_diagonal_residual =
        result.green.diagonal().norm() / greenScale;
    if (!std::isfinite(result.green_skew_residual) ||
        !std::isfinite(result.green_diagonal_residual) ||
        (options.fail_on_green_structure &&
         (result.green_skew_residual > options.green_residual_tolerance ||
          result.green_diagonal_residual > options.green_residual_tolerance))) {
        result.status = PureProjectorStatus::green_structure_exceeded;
        return result;
    }
    result.status = PureProjectorStatus::success;
    return result;
}

inline PureProjectorGreenResult pureProjectorGreenThinQr(
    const MatType &phiRight, const MatType &phiLeft,
    const PureProjectorOptions &options = PureProjectorOptions()) {
    const ThinQrResult rightQr = thinQrSubspace(phiRight, options);
    if (!rightQr.ok()) {
        PureProjectorGreenResult result;
        result.status = rightQr.status;
        return result;
    }
    const ThinQrResult leftQr = thinQrSubspace(phiLeft, options);
    if (!leftQr.ok()) {
        PureProjectorGreenResult result;
        result.status = leftQr.status;
        return result;
    }
    return pureProjectorGreen(rightQr.q, leftQr.q, options);
}

inline bool propagateKet(const MatType &factor, MatType &subspace) {
    if (factor.rows() <= 0 || factor.rows() != factor.cols() ||
        factor.cols() != subspace.rows() || !pureProjectorMatrixFinite(factor) ||
        !pureProjectorMatrixFinite(subspace)) return false;
    subspace = factor * subspace;
    return pureProjectorMatrixFinite(subspace);
}

inline bool propagateBra(const MatType &factor, MatType &subspace) {
    if (factor.rows() <= 0 || factor.rows() != factor.cols() ||
        factor.cols() != subspace.rows() || !pureProjectorMatrixFinite(factor) ||
        !pureProjectorMatrixFinite(subspace)) return false;
    subspace = factor.adjoint() * subspace;
    return pureProjectorMatrixFinite(subspace);
}

inline bool propagateKet(Operator &factor, MatType &subspace) {
    if (!pureProjectorMatrixFinite(subspace)) return false;
    MatType output;
    factor.left_multiply(subspace, output);
    if (output.rows() != subspace.rows() || output.cols() != subspace.cols() ||
        !pureProjectorMatrixFinite(output)) return false;
    subspace.swap(output);
    return true;
}

inline bool propagateBra(Operator &factor, MatType &subspace) {
    if (!pureProjectorMatrixFinite(subspace)) return false;
    MatType output;
    factor.adjoint_left_multiply(subspace, output);
    if (output.rows() != subspace.rows() || output.cols() != subspace.cols() ||
        !pureProjectorMatrixFinite(output)) return false;
    subspace.swap(output);
    return true;
}

#endif
