#ifndef PURE_STATE_DRIVEN_STAGE1_DENSE_H
#define PURE_STATE_DRIVEN_STAGE1_DENSE_H

#include "stage1_common.h"

namespace driven_stage1 {

inline MatType gammaMatrix(int modes, int which) {
    const int dimension = 1 << modes;
    const int site = which % modes;
    const bool second = which >= modes;
    MatType gamma = MatType::Zero(dimension, dimension);
    for (int state = 0; state < dimension; ++state) {
        const int flipped = state ^ (1 << site);
        const bool odd = (__builtin_popcount(
            unsigned(state & ((1 << site) - 1))) & 1) != 0;
        const double jordanWigner = odd ? -1.0 : 1.0;
        gamma(flipped, state) = second ?
            DataType(0.0, (state & (1 << site)) ?
                jordanWigner : -jordanWigner) :
            DataType(jordanWigner, 0.0);
    }
    return gamma;
}

inline std::vector<MatType> gammaMatrices(int length) {
    std::vector<MatType> result;
    for (int index = 0; index < 2 * length; ++index)
        result.push_back(gammaMatrix(length, index));
    return result;
}

inline MatType denseQuadratic(const MatType &generator,
                              const std::vector<MatType> &gamma) {
    MatType result = MatType::Zero(gamma.front().rows(), gamma.front().cols());
    for (int row = 0; row < generator.rows(); ++row)
        for (int column = row + 1; column < generator.cols(); ++column)
            result += 0.5 * generator(row, column) *
                      gamma[row] * gamma[column];
    return result;
}

inline cVecType normalized(cVecType state) {
    const double norm = state.norm();
    require(std::isfinite(norm) && norm > 0.0,
            "dense state has zero or nonfinite norm");
    return state / norm;
}

struct DenseInitialState {
    cVecType vector;
    int parity = 0;
    double energy = 0.0;
};

inline DenseInitialState denseInitialState(
        const ModelParameters &p, const std::vector<MatType> &gamma,
        const GaussianTrialState &gaussianInitial, double trialMu,
        int requiredParity, double trialT = 1.0,
        double trialDelta = 1.0, bool requireGround = true) {
    const MatType h = -denseQuadratic(
        kineticGenerator(p, trialT, trialDelta, trialMu), gamma);
    Eigen::SelfAdjointEigenSolver<MatType> solver(h);
    require(solver.info() == Eigen::Success,
            "dense Phi_0 eigensolver failed");
    DenseInitialState result;
    double bestGreenError = std::numeric_limits<double>::infinity();
    int selectedColumn = -1;
    for (int column = 0; column < solver.eigenvectors().cols(); ++column) {
        const cVecType candidate = solver.eigenvectors().col(column);
        MatType candidateGreen = MatType::Zero(2 * p.L, 2 * p.L);
        for (int row = 0; row < candidateGreen.rows(); ++row)
            for (int col = 0; col < candidateGreen.cols(); ++col)
                if (row != col)
                    candidateGreen(row, col) =
                        -(candidate.adjoint() * gamma[row] * gamma[col] *
                          candidate)(0);
        const double error = relativeError(
            candidateGreen, gaussianInitial.G_T);
        if (error < bestGreenError) {
            bestGreenError = error;
            selectedColumn = column;
        }
    }
    if (selectedColumn >= 0 && bestGreenError < 1e-10) {
        result.vector = solver.eigenvectors().col(selectedColumn);
    } else {
        // A Gaussian state can be a particular superposition inside a
        // degenerate eigenspace, while SelfAdjointEigenSolver returns another
        // arbitrary basis.  Its flattened covariance is an unambiguous gapped
        // parent Hamiltonian and reconstructs the matching Fock vector.
        const MatType parent = -denseQuadratic(gaussianInitial.G_T, gamma);
        Eigen::SelfAdjointEigenSolver<MatType> parentSolver(parent);
        require(parentSolver.info() == Eigen::Success,
                "dense covariance-parent eigensolver failed");
        result.vector = parentSolver.eigenvectors().col(0);
        MatType parentGreen = MatType::Zero(2 * p.L, 2 * p.L);
        for (int row = 0; row < parentGreen.rows(); ++row)
            for (int col = 0; col < parentGreen.cols(); ++col)
                if (row != col)
                    parentGreen(row, col) =
                        -(result.vector.adjoint() * gamma[row] * gamma[col] *
                          result.vector)(0);
        require(relativeError(parentGreen, gaussianInitial.G_T) < 1e-10,
                "no dense Fock vector matches Phi_0 covariance");
    }
    result.energy = (result.vector.adjoint() * h * result.vector)(0).real();
    require((h * result.vector - result.energy * result.vector).norm() < 1e-9,
            "Phi_0 covariance is not an eigenstate of the requested dense H_T");
    if (requireGround)
        require(std::abs(result.energy - solver.eigenvalues()(0)) < 1e-10,
                "Phi_0 covariance does not select the physical dense ground state");
    double parity = 0.0;
    for (int state = 0; state < result.vector.size(); ++state)
        parity += ((__builtin_popcount(unsigned(state)) & 1) ? -1.0 : 1.0) *
                  std::norm(result.vector(state));
    require(std::abs(std::abs(parity) - 1.0) < 1e-9,
            "dense Phi_0 has ambiguous fermion parity");
    result.parity = parity >= 0.0 ? 1 : -1;
    require(requiredParity == 0 || result.parity == requiredParity,
            "dense Phi_0 parity differs from Gaussian Phi_0");
    return result;
}

inline MatType denseFactorAtFlattenedLocation(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        const std::vector<MatType> &gamma, const PureSliceLocation &location,
        int field) {
    MatType generator;
    double scale = 1.0;
    if (location.aux < 0) {
        generator = kineticGenerator(p, p.t, p.delta, p.mu);
        scale = -0.5 * p.dt;
    } else {
        generator = localHsGenerator(
            p, protocol.midpointValue(location.slice), location.bond,
            location.aux, field);
    }
    // Keep the Fock mapping identical to the validated Phase 2 convention:
    // exp(-1/2 A_ij gamma_i gamma_j) implements B=exp(A).
    MatType dense = exponential(denseQuadratic(generator, gamma), -scale);
    if (location.branch == PureBranch::Bra) dense.adjointInPlace();
    return dense;
}

struct DenseContourResult {
    DataType weight = 0.0;
    double ket_norm = 0.0;
    double bra_norm = 0.0;
    MatType green;
};

inline DenseContourResult denseContour(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        const GaussianTrialState &gaussianInitial,
        const cVecType &initial, const PureFastConfiguration &configuration) {
    const std::vector<MatType> gamma = gammaMatrices(p.L);
    const int center = int(configuration.slices.size() / 2);
    MatType ketProduct = MatType::Identity(initial.size(), initial.size());
    MatType braProduct = MatType::Identity(initial.size(), initial.size());
    for (int index = 0; index < int(configuration.slices.size()); ++index) {
        const MatType stored = denseFactorAtFlattenedLocation(
            p, protocol, gamma, configuration.locations[index],
            configuration.hs_fields[index]);
        if (index < center) ketProduct = ketProduct * stored;
        else braProduct = braProduct * stored;
    }
    DenseContourResult result;
    result.weight =
        (initial.adjoint() * ketProduct * braProduct * initial)(0);
    result.ket_norm = (ketProduct * initial).norm();
    result.bra_norm = (braProduct * initial).norm();
    require(std::abs(result.weight) > 1e-14,
            "dense contour has zero overlap weight");
    result.green = MatType::Zero(2 * p.L, 2 * p.L);
    for (int row = 0; row < result.green.rows(); ++row)
        for (int column = 0; column < result.green.cols(); ++column)
            if (row != column)
                result.green(row, column) =
                    -(initial.adjoint() * ketProduct * gamma[row] *
                      gamma[column] * braProduct * initial)(0) /
                    result.weight;
    (void)gaussianInitial;
    return result;
}

inline MatType denseHamiltonian(const ModelParameters &p, double interaction,
                                const std::vector<MatType> &gamma) {
    MatType result = -denseQuadratic(
        kineticGenerator(p, p.t, p.delta, p.mu), gamma);
    std::vector<MatType> density(p.L);
    for (int site = 0; site < p.L; ++site)
        density[site] = DataType(0.0, -0.5) *
                        gamma[site] * gamma[p.L + site];
    const int bonds = p.boundary == 0 ? p.L : p.L - 1;
    for (int bond = 0; bond < bonds; ++bond)
        result += interaction * density[bond] * density[(bond + 1) % p.L];
    return result;
}

inline MatType denseInteractionLayer(
        const ModelParameters &p, const std::vector<MatType> &gamma,
        int layer) {
    std::vector<MatType> density(p.L);
    for (int site = 0; site < p.L; ++site)
        density[site] = DataType(0.0, -0.5) *
                        gamma[site] * gamma[p.L + site];
    MatType result = MatType::Zero(gamma.front().rows(), gamma.front().cols());
    const int bonds = p.boundary == 0 ? p.L : p.L - 1;
    for (int bond = 0; bond < bonds; ++bond)
        if ((bond & 1) == layer)
            result += density[bond] * density[(bond + 1) % p.L];
    return result;
}

inline cVecType denseSameContourState(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        const std::vector<MatType> &gamma, const cVecType &initial,
        double *denominator = nullptr) {
    const MatType kinetic = -denseQuadratic(
        kineticGenerator(p, p.t, p.delta, p.mu), gamma);
    const MatType even = denseInteractionLayer(p, gamma, 0);
    const MatType odd = denseInteractionLayer(p, gamma, 1);
    const MatType halfKinetic = exponential(
        kinetic, -0.5 * protocol.deltaTau());
    MatType contour = MatType::Identity(initial.size(), initial.size());
    for (double interaction : protocol.midpointValues()) {
        const MatType slice = halfKinetic * exponential(
            even, -protocol.deltaTau() * interaction) * exponential(
            odd, -protocol.deltaTau() * interaction) * halfKinetic;
        contour = slice * contour;
    }
    cVecType state = contour * initial;
    const double normSquared = state.squaredNorm();
    require(std::isfinite(normSquared) && normSquared > 0.0,
            "same-contour dense denominator failed");
    if (denominator) *denominator = normSquared;
    return state / std::sqrt(normSquared);
}

inline MatType denseHsSummedLayer(
        const ModelParameters &p, double interaction,
        const std::vector<MatType> &gamma, int layer) {
    const auto counts = pureProjectorCheckerboardBondCounts(p.L, p.boundary);
    const int count = layer == 0 ? counts.first : counts.second;
    MatType sum = MatType::Zero(gamma.front().rows(), gamma.front().cols());
    for (int mask = 0; mask < (1 << count); ++mask) {
        MatType product = MatType::Identity(sum.rows(), sum.cols());
        for (int auxiliary = 0; auxiliary < count; ++auxiliary) {
            const int sigma = (mask & (1 << auxiliary)) ? 1 : -1;
            const MatType generator = localHsGenerator(
                p, interaction, layer, auxiliary, sigma);
            product = product * exponential(
                denseQuadratic(generator, gamma), -1.0);
        }
        sum += product;
    }
    // Per bond:
    // exp[-Delta_tau V rho_i rho_j]
    //   = exp[-Delta_tau V/4] / 2 * sum_sigma O_sigma.
    return std::pow(0.5 * std::exp(-0.25 * p.dt * interaction), count) *
           sum;
}

inline cVecType denseHsSameContourState(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        const std::vector<MatType> &gamma, const cVecType &initial,
        double *denominator = nullptr) {
    const MatType kinetic = -denseQuadratic(
        kineticGenerator(p, p.t, p.delta, p.mu), gamma);
    const MatType halfKinetic = exponential(kinetic, -0.5 * p.dt);
    MatType contour = MatType::Identity(initial.size(), initial.size());
    for (double interaction : protocol.midpointValues()) {
        const MatType slice = halfKinetic * denseHsSummedLayer(
            p, interaction, gamma, 0) * denseHsSummedLayer(
            p, interaction, gamma, 1) * halfKinetic;
        contour = slice * contour;
    }
    cVecType state = contour * initial;
    const double normSquared = state.squaredNorm();
    require(std::isfinite(normSquared) && normSquared > 0.0,
            "HS-reconstructed dense denominator failed");
    if (denominator) *denominator = normSquared;
    return state / std::sqrt(normSquared);
}

inline MatType denseStateGreen(
        const cVecType &state, const std::vector<MatType> &gamma) {
    MatType green = MatType::Zero(gamma.size(), gamma.size());
    for (int row = 0; row < green.rows(); ++row)
        for (int column = 0; column < green.cols(); ++column)
            if (row != column)
                green(row, column) =
                    -(state.adjoint() * gamma[row] * gamma[column] * state)(0);
    return green;
}

}  // namespace driven_stage1

#endif
