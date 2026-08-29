#ifndef PURE_PROJECTOR_WEIGHT_H
#define PURE_PROJECTOR_WEIGHT_H

#include "gaussian_trial_state.h"
#include "pure_projector_green.h"
#include "skewMatUtils.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <utility>
#include <vector>

enum class PureProjectorWeightMode { RealZ2, GenericComplex };

enum class PureProjectorWeightStatus {
    success,
    invalid_slice,
    propagation_failure,
    overlap_untrusted,
    pfaffian_untrusted,
    complex_ratio,
    zero_weight,
    determinant_identity_failure
};

inline const char *pureProjectorWeightStatusName(PureProjectorWeightStatus status) {
    switch (status) {
    case PureProjectorWeightStatus::success: return "success";
    case PureProjectorWeightStatus::invalid_slice: return "invalid_slice";
    case PureProjectorWeightStatus::propagation_failure: return "propagation_failure";
    case PureProjectorWeightStatus::overlap_untrusted: return "overlap_untrusted";
    case PureProjectorWeightStatus::pfaffian_untrusted: return "pfaffian_untrusted";
    case PureProjectorWeightStatus::complex_ratio: return "complex_ratio";
    case PureProjectorWeightStatus::zero_weight: return "zero_weight";
    case PureProjectorWeightStatus::determinant_identity_failure:
        return "determinant_identity_failure";
    }
    return "unknown";
}

struct PureProjectorSlice {
    MatType matrix;
    DataType eta = DataType(1.0, 0.0);
    std::string label;

    PureProjectorSlice() = default;
    PureProjectorSlice(const MatType &factor, DataType scalar, std::string name)
        : matrix(factor), eta(scalar), label(std::move(name)) {}
};

struct PureProjectorWeightOptions {
    PureProjectorWeightMode mode = PureProjectorWeightMode::GenericComplex;
    PureProjectorOptions green_options;
    int thin_qr_interval = 8;
    double reality_tolerance = 1e-10;
    double zero_tolerance = 1e-14;
    double determinant_tolerance = 1e-9;
};

struct PureProjectorWeightResult {
    PureProjectorWeightStatus status = PureProjectorWeightStatus::invalid_slice;
    double log_abs_weight = -std::numeric_limits<double>::infinity();
    DataType complex_phase = DataType(0.0, 0.0);
    int z2_sign = 0;
    MatType green;
    int overlap_rank = 0;
    double overlap_rcond = 0.0;
    double overlap_residual = std::numeric_limits<double>::infinity();
    double green_residual = std::numeric_limits<double>::infinity();
    PfaffianStatus pfaffian_status = PfaffianStatus::success;
    int first_failing_slice = -1;
    double determinant_identity_error = 0.0;
    DataType weight = DataType(0.0, 0.0);
    DataType determinant_rhs = DataType(1.0, 0.0);
    bool ok() const { return status == PureProjectorWeightStatus::success; }
};

inline DataType pureProjectorUnitPhase(DataType value) {
    const double magnitude = std::abs(value);
    return magnitude > 0.0 ? value / magnitude : DataType(0.0, 0.0);
}

class PureProjectorWeightEvaluator {
public:
    explicit PureProjectorWeightEvaluator(
        PureProjectorWeightOptions options = PureProjectorWeightOptions())
        : options_(options) {}

    PureProjectorWeightResult evaluate(
        const GaussianTrialState &trial,
        const std::vector<PureProjectorSlice> &slices) const {
        PureProjectorWeightResult result;
        result.status = PureProjectorWeightStatus::success;
        result.log_abs_weight = 0.0;
        result.complex_phase = DataType(1.0, 0.0);
        result.z2_sign = 1;
        result.weight = DataType(1.0, 0.0);
        result.green = trial.G_T;
        result.overlap_rank = trial.Phi.cols();
        result.overlap_rcond = 1.0;
        result.overlap_residual = 0.0;
        result.green_residual = trial.diagnostics.green_involution_error;

        MatType rawPhi = trial.Phi;
        MatType stablePhi = trial.Phi;
        DataType accumulatedEta(1.0, 0.0);
        DataType previousDet(1.0, 0.0);

        for (int index = 0; index < int(slices.size()); ++index) {
            const PureProjectorSlice &slice = slices[index];
            if (slice.matrix.rows() != trial.Phi.rows() ||
                slice.matrix.cols() != trial.Phi.rows() ||
                !pureProjectorMatrixFinite(slice.matrix) ||
                !std::isfinite(slice.eta.real()) || !std::isfinite(slice.eta.imag())) {
                return fail(result, PureProjectorWeightStatus::invalid_slice, index);
            }

            MatType nextRaw = slice.matrix * rawPhi;
            MatType nextStable = slice.matrix * stablePhi;
            if (!pureProjectorMatrixFinite(nextRaw) || !pureProjectorMatrixFinite(nextStable))
                return fail(result, PureProjectorWeightStatus::propagation_failure, index);

            const MatType overlap = trial.Phi.adjoint() * nextRaw;
            Eigen::JacobiSVD<MatType> overlapSvd(
                overlap, Eigen::ComputeThinU | Eigen::ComputeThinV);
            overlapSvd.setThreshold(options_.green_options.rank_tolerance);
            result.overlap_rank = overlapSvd.rank();
            if (result.overlap_rank != overlap.rows()) {
                result.pfaffian_status = PfaffianStatus::zero_pivot;
                return zero(result, index);
            }
            const double largest = overlapSvd.singularValues()(0);
            const double smallest = overlapSvd.singularValues()(overlap.rows() - 1);
            result.overlap_rcond = largest > 0.0 ? smallest / largest : 0.0;
            if (!std::isfinite(result.overlap_rcond) ||
                result.overlap_rcond < options_.green_options.minimum_overlap_rcond)
                return fail(result, PureProjectorWeightStatus::overlap_untrusted, index);

            const DataType nextDet = overlap.determinant();
            if (!std::isfinite(nextDet.real()) || !std::isfinite(nextDet.imag()) ||
                std::abs(previousDet) <= options_.zero_tolerance)
                return zero(result, index);
            DataType ratioSquared = slice.eta * slice.eta * nextDet / previousDet;
            if (std::abs(ratioSquared) <= options_.zero_tolerance) return zero(result, index);
            DataType ratio = std::sqrt(ratioSquared);

            // The Pfaffian multiplication formula fixes the otherwise ambiguous
            // square-root branch.  Identity factors have a singular Cayley Green
            // and are anchored continuously at eta instead.
            const double identityResidual =
                (slice.matrix - MatType::Identity(slice.matrix.rows(), slice.matrix.cols())).norm();
            if (identityResidual <= 10.0 * options_.zero_tolerance) {
                ratio = slice.eta;
            } else {
                MatType onePlus = MatType::Identity(slice.matrix.rows(), slice.matrix.cols()) +
                                  slice.matrix;
                Eigen::FullPivLU<MatType> lu(onePlus);
                if (!lu.isInvertible()) {
                    result.pfaffian_status = PfaffianStatus::zero_pivot;
                    return fail(result, PureProjectorWeightStatus::pfaffian_untrusted, index);
                }
                MatType sliceGreen = lu.solve(2.0 * MatType::Identity(
                    slice.matrix.rows(), slice.matrix.cols())) -
                    MatType::Identity(slice.matrix.rows(), slice.matrix.cols());
                const PureProjectorGreenResult partial =
                    pureProjectorGreenThinQr(stablePhi, trial.Phi, options_.green_options);
                if (!partial.ok())
                    return fail(result, PureProjectorWeightStatus::overlap_untrusted, index);
                const PfaffianResult pf =
                    pfaffianForSignOfProductWithStatus(sliceGreen, partial.green);
                result.pfaffian_status = pf.status;
                if (!pf.ok())
                    return fail(result, PureProjectorWeightStatus::pfaffian_untrusted, index);
                const DataType target = pureProjectorUnitPhase(slice.eta * pf.value);
                if ((std::conj(target) * ratio).real() < 0.0) ratio = -ratio;
            }

            if (!std::isfinite(ratio.real()) || !std::isfinite(ratio.imag()) ||
                std::abs(ratio) <= options_.zero_tolerance)
                return zero(result, index);
            if (options_.mode == PureProjectorWeightMode::RealZ2) {
                if (std::abs(ratio.imag()) >
                    options_.reality_tolerance * std::max(1.0, std::abs(ratio.real())))
                    return fail(result, PureProjectorWeightStatus::complex_ratio, index);
                if (std::abs(ratio.real()) <= options_.zero_tolerance) return zero(result, index);
                result.z2_sign *= ratio.real() >= 0.0 ? 1 : -1;
                result.complex_phase = DataType(double(result.z2_sign), 0.0);
            } else {
                result.complex_phase *= pureProjectorUnitPhase(ratio);
                result.complex_phase = pureProjectorUnitPhase(result.complex_phase);
                if (std::abs(result.complex_phase.imag()) <= options_.reality_tolerance)
                    result.z2_sign = result.complex_phase.real() >= 0.0 ? 1 : -1;
                else result.z2_sign = 0;
            }
            result.log_abs_weight += std::log(std::abs(ratio));
            rawPhi.swap(nextRaw);
            stablePhi.swap(nextStable);
            accumulatedEta *= slice.eta;
            previousDet = nextDet;

            if (options_.thin_qr_interval > 0 &&
                ((index + 1) % options_.thin_qr_interval) == 0) {
                const ThinQrResult qr = thinQrSubspace(stablePhi, options_.green_options);
                if (!qr.ok())
                    return fail(result, PureProjectorWeightStatus::overlap_untrusted, index);
                stablePhi = qr.q;
            }
        }

        const PureProjectorGreenResult finalGreen =
            pureProjectorGreenThinQr(stablePhi, trial.Phi, options_.green_options);
        if (!finalGreen.ok())
            return fail(result, PureProjectorWeightStatus::overlap_untrusted,
                        slices.empty() ? -1 : int(slices.size()) - 1);
        result.green = finalGreen.green;
        result.overlap_rank = finalGreen.overlap_rank;
        result.overlap_rcond = finalGreen.overlap_rcond;
        result.overlap_residual = finalGreen.solve_residual;
        result.green_residual = finalGreen.green_residual;
        result.weight = std::exp(result.log_abs_weight) * result.complex_phase;
        result.determinant_rhs = accumulatedEta * accumulatedEta * previousDet;
        result.determinant_identity_error =
            std::abs(result.weight * result.weight - result.determinant_rhs) /
            std::max(options_.zero_tolerance, std::abs(result.determinant_rhs));
        if (result.determinant_identity_error > options_.determinant_tolerance)
            return fail(result, PureProjectorWeightStatus::determinant_identity_failure,
                        slices.empty() ? -1 : int(slices.size()) - 1);
        return result;
    }

private:
    PureProjectorWeightOptions options_;

    static PureProjectorWeightResult fail(PureProjectorWeightResult result,
                                          PureProjectorWeightStatus status,
                                          int slice) {
        result.status = status;
        result.first_failing_slice = slice;
        return result;
    }

    static PureProjectorWeightResult zero(PureProjectorWeightResult result, int slice) {
        result.status = PureProjectorWeightStatus::zero_weight;
        result.first_failing_slice = slice;
        result.log_abs_weight = -std::numeric_limits<double>::infinity();
        result.complex_phase = DataType(0.0, 0.0);
        result.z2_sign = 0;
        result.weight = DataType(0.0, 0.0);
        return result;
    }
};

#endif
