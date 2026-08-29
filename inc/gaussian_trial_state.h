#ifndef GAUSSIAN_TRIAL_STATE_H
#define GAUSSIAN_TRIAL_STATE_H

#include "kitaevChain.h"
#include "skewMatUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

enum class GaussianTrialStateStatus {
    success,
    invalid_dimension,
    nonfinite_input,
    nonhermitian_hamiltonian,
    nonskew_hamiltonian,
    eigensolver_failure,
    zero_mode_ambiguous,
    invalid_occupied_subspace,
    invariant_failure,
    parity_failure
};

inline const char *gaussianTrialStateStatusName(GaussianTrialStateStatus status) {
    switch (status) {
    case GaussianTrialStateStatus::success: return "success";
    case GaussianTrialStateStatus::invalid_dimension: return "invalid_dimension";
    case GaussianTrialStateStatus::nonfinite_input: return "nonfinite_input";
    case GaussianTrialStateStatus::nonhermitian_hamiltonian: return "nonhermitian_hamiltonian";
    case GaussianTrialStateStatus::nonskew_hamiltonian: return "nonskew_hamiltonian";
    case GaussianTrialStateStatus::eigensolver_failure: return "eigensolver_failure";
    case GaussianTrialStateStatus::zero_mode_ambiguous: return "zero_mode_ambiguous";
    case GaussianTrialStateStatus::invalid_occupied_subspace: return "invalid_occupied_subspace";
    case GaussianTrialStateStatus::invariant_failure: return "invariant_failure";
    case GaussianTrialStateStatus::parity_failure: return "parity_failure";
    }
    return "unknown";
}

class GaussianTrialStateError : public std::runtime_error {
public:
    GaussianTrialStateError(GaussianTrialStateStatus status, const std::string &message)
        : std::runtime_error(message), status_(status) {}
    GaussianTrialStateStatus code() const { return status_; }
private:
    GaussianTrialStateStatus status_;
};

struct GaussianTrialStateOptions {
    double validation_tolerance = 1e-10;
    double zero_mode_tolerance = 1e-10;
};

struct GaussianTrialStateDiagnostics {
    GaussianTrialStateStatus status = GaussianTrialStateStatus::invalid_dimension;
    int rows = 0;
    int occupied_columns = 0;
    int zero_mode_count = 0;
    double orthonormal_error = std::numeric_limits<double>::infinity();
    double isotropy_error = std::numeric_limits<double>::infinity();
    double green_skew_error = std::numeric_limits<double>::infinity();
    double green_involution_error = std::numeric_limits<double>::infinity();
    double hamiltonian_hermitian_error = 0.0;
    double hamiltonian_skew_error = 0.0;
    bool finite = false;
    bool valid = false;
};

class GaussianTrialState {
public:
    MatType Phi;
    MatType G_T;
    GaussianTrialStateDiagnostics diagnostics;

    static GaussianTrialState fromPhi(
        const MatType &phi,
        const GaussianTrialStateOptions &options = GaussianTrialStateOptions()) {
        GaussianTrialState state;
        state.Phi = phi;
        state.diagnostics = validatePhi(phi, options.validation_tolerance, &state.G_T);
        if (!state.diagnostics.valid) {
            throw GaussianTrialStateError(
                state.diagnostics.status,
                std::string("invalid pure Gaussian occupied subspace: ") +
                    gaussianTrialStateStatusName(state.diagnostics.status));
        }
        return state;
    }

    static GaussianTrialState fromMajoranaHamiltonian(
        const MatType &hamiltonian,
        const GaussianTrialStateOptions &options = GaussianTrialStateOptions()) {
        const int dimension = hamiltonian.rows();
        if (dimension <= 0 || hamiltonian.cols() != dimension || (dimension % 2) != 0) {
            throw GaussianTrialStateError(GaussianTrialStateStatus::invalid_dimension,
                                          "Majorana Hamiltonian must be even and square");
        }
        if (!matrixFinite(hamiltonian)) {
            throw GaussianTrialStateError(GaussianTrialStateStatus::nonfinite_input,
                                          "Majorana Hamiltonian contains nonfinite values");
        }
        const double scale = std::max(1.0, hamiltonian.norm());
        const double hermitianError = (hamiltonian - hamiltonian.adjoint()).norm() / scale;
        const double skewError = (hamiltonian + hamiltonian.transpose()).norm() / scale;
        if (hermitianError > options.validation_tolerance) {
            throw GaussianTrialStateError(GaussianTrialStateStatus::nonhermitian_hamiltonian,
                                          "Majorana Hamiltonian is not Hermitian");
        }
        if (skewError > options.validation_tolerance) {
            throw GaussianTrialStateError(GaussianTrialStateStatus::nonskew_hamiltonian,
                                          "Majorana Hamiltonian is not transpose-skew");
        }

        Eigen::SelfAdjointEigenSolver<MatType> solver(hamiltonian);
        if (solver.info() != Eigen::Success) {
            throw GaussianTrialStateError(GaussianTrialStateStatus::eigensolver_failure,
                                          "Majorana Hamiltonian eigensolver failed");
        }
        int zeroModes = 0;
        int negativeModes = 0;
        for (int index = 0; index < dimension; ++index) {
            const double eigenvalue = solver.eigenvalues()(index);
            if (std::abs(eigenvalue) <= options.zero_mode_tolerance) ++zeroModes;
            else if (eigenvalue < 0.0) ++negativeModes;
        }
        if (zeroModes != 0) {
            throw GaussianTrialStateError(
                GaussianTrialStateStatus::zero_mode_ambiguous,
                "Majorana Hamiltonian has zero-mode ambiguity; an explicit parity policy or "
                "edge splitting is required");
        }
        if (negativeModes != dimension / 2) {
            throw GaussianTrialStateError(GaussianTrialStateStatus::invalid_occupied_subspace,
                                          "negative-energy subspace is not half-dimensional");
        }

        GaussianTrialState state = fromPhi(
            solver.eigenvectors().leftCols(dimension / 2), options);
        state.diagnostics.zero_mode_count = zeroModes;
        state.diagnostics.hamiltonian_hermitian_error = hermitianError;
        state.diagnostics.hamiltonian_skew_error = skewError;
        return state;
    }

    static GaussianTrialState fromKitaevChain(
        const SpinlessTvChainUtils &config, double edgeSplitting = 0.0,
        const GaussianTrialStateOptions &options = GaussianTrialStateOptions()) {
        MatType hamiltonian = MatType::Zero(config.nDim, config.nDim);
        config.KineticGenerator(hamiltonian);
        if (edgeSplitting != 0.0) {
            const int leftEdge = config.majoranaCoord2Idx(0, 1);
            const int rightEdge = config.majoranaCoord2Idx(config.Lx - 1, 1);
            const DataType splitting(0.0, edgeSplitting);
            hamiltonian(leftEdge, rightEdge) += splitting;
            hamiltonian(rightEdge, leftEdge) -= splitting;
        }
        return fromMajoranaHamiltonian(hamiltonian, options);
    }

    int fermionParity(double realityTolerance = 1e-8) const {
        MatType realSkew = DataType(0.0, -1.0) * G_T;
        const PfaffianResult result = signOfPfafWithStatus(realSkew);
        if (!result.ok() || std::abs(result.value.imag()) > realityTolerance ||
            std::abs(result.value.real()) < 1.0 - realityTolerance) {
            throw GaussianTrialStateError(GaussianTrialStateStatus::parity_failure,
                                          "pure-state parity Pfaffian is untrusted");
        }
        return result.value.real() >= 0.0 ? 1 : -1;
    }

private:
    static bool matrixFinite(const MatType &matrix) {
        for (int col = 0; col < matrix.cols(); ++col) {
            for (int row = 0; row < matrix.rows(); ++row) {
                if (!std::isfinite(matrix(row, col).real()) ||
                    !std::isfinite(matrix(row, col).imag())) return false;
            }
        }
        return true;
    }

    static GaussianTrialStateDiagnostics validatePhi(
        const MatType &phi, double tolerance, MatType *green) {
        GaussianTrialStateDiagnostics result;
        result.rows = phi.rows();
        result.occupied_columns = phi.cols();
        if (phi.rows() <= 0 || phi.cols() <= 0 || phi.rows() != 2 * phi.cols()) {
            result.status = GaussianTrialStateStatus::invalid_dimension;
            return result;
        }
        result.finite = matrixFinite(phi);
        if (!result.finite) {
            result.status = GaussianTrialStateStatus::nonfinite_input;
            return result;
        }
        const int dimension = phi.rows();
        const int occupied = phi.cols();
        const double phiScale = std::max(1.0, std::sqrt(double(occupied)));
        result.orthonormal_error =
            (phi.adjoint() * phi - MatType::Identity(occupied, occupied)).norm() / phiScale;
        result.isotropy_error = (phi.transpose() * phi).norm() / phiScale;
        *green = MatType::Identity(dimension, dimension) - 2.0 * phi * phi.adjoint();
        const double greenScale = std::max(1.0, std::sqrt(double(dimension)));
        result.green_skew_error = (green->transpose() + *green).norm() / greenScale;
        result.green_involution_error =
            ((*green) * (*green) - MatType::Identity(dimension, dimension)).norm() /
            greenScale;
        result.finite = result.finite && matrixFinite(*green);
        result.valid = result.finite &&
            std::max({result.orthonormal_error, result.isotropy_error,
                      result.green_skew_error, result.green_involution_error}) <= tolerance;
        result.status = result.valid ? GaussianTrialStateStatus::success
                                     : GaussianTrialStateStatus::invariant_failure;
        return result;
    }
};

#endif
