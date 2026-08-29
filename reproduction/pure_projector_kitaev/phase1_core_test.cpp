#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "gaussian_trial_state.h"
#include "pure_projector_green.h"

namespace {

constexpr double kTolerance = 1e-10;

double relativeError(const MatType &actual, const MatType &expected) {
    return (actual - expected).norm() / std::max(1.0, expected.norm());
}

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

MatType canonicalPhi(int modes, bool filled = false) {
    MatType phi = MatType::Zero(2 * modes, modes);
    const DataType phase(0.0, filled ? -1.0 : 1.0);
    for (int mode = 0; mode < modes; ++mode) {
        phi(2 * mode, mode) = 1.0 / std::sqrt(2.0);
        phi(2 * mode + 1, mode) = phase / std::sqrt(2.0);
    }
    return phi;
}

MatType randomLegalPhi(int modes, unsigned seed) {
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    Eigen::MatrixXd raw(2 * modes, 2 * modes);
    for (int row = 0; row < raw.rows(); ++row)
        for (int col = 0; col < raw.cols(); ++col) raw(row, col) = normal(generator);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(raw);
    const Eigen::MatrixXd orthogonal =
        qr.householderQ() * Eigen::MatrixXd::Identity(2 * modes, 2 * modes);
    return orthogonal.cast<DataType>() * canonicalPhi(modes);
}

MatType moderateMajoranaFactor(int dimension, int first, int second, double theta) {
    MatType factor = MatType::Identity(dimension, dimension);
    const double c = std::cosh(theta);
    const double s = std::sinh(theta);
    factor(first, first) = c;
    factor(second, second) = c;
    factor(first, second) = DataType(0.0, -s);
    factor(second, first) = DataType(0.0, s);
    return factor;
}

class MatrixTestOperator : public Operator {
public:
    explicit MatrixTestOperator(const MatType &matrix) : matrix_(matrix) {}
    void left_multiply(const MatType &input, MatType &output) override {
        output = matrix_ * input;
    }
    void adjoint_left_multiply(const MatType &input, MatType &output) override {
        output = matrix_.adjoint() * input;
    }
    DataType update(MatType &) override { return 1.0; }
    void getGreensMat(MatType &) override {}
private:
    MatType matrix_;
};

void printResult(const std::string &test, const std::map<std::string, double> &metrics,
                 const std::map<std::string, std::string> &text = {}) {
    std::cout << std::setprecision(17) << "{\"test\":\"" << test
              << "\",\"status\":\"PASS\"";
    for (const auto &entry : metrics) {
        std::cout << ",\"" << entry.first << "\":" << entry.second;
    }
    for (const auto &entry : text) {
        std::cout << ",\"" << entry.first << "\":\"" << entry.second << "\"";
    }
    std::cout << "}\n";
}

void testSingleSite() {
    const GaussianTrialState empty = GaussianTrialState::fromPhi(canonicalPhi(1));
    const GaussianTrialState filled = GaussianTrialState::fromPhi(canonicalPhi(1, true));
    MatType expectedEmpty(2, 2), expectedFilled(2, 2);
    expectedEmpty << 0.0, DataType(0.0, 1.0), DataType(0.0, -1.0), 0.0;
    expectedFilled = -expectedEmpty;
    const double phiError = std::max(relativeError(empty.Phi, canonicalPhi(1)),
                                     relativeError(filled.Phi, canonicalPhi(1, true)));
    const double greenError = std::max(relativeError(empty.G_T, expectedEmpty),
                                       relativeError(filled.G_T, expectedFilled));
    require(empty.diagnostics.valid && filled.diagnostics.valid, "single-site diagnostics failed");
    require(phiError < kTolerance && greenError < kTolerance, "single-site analytic mismatch");
    require(empty.fermionParity() == -filled.fermionParity(), "single-site parity mismatch");
    printResult("single_site", {{"phi_error", phiError}, {"green_error", greenError}});
}

void testIdentityPropagation() {
    const GaussianTrialState trial = GaussianTrialState::fromPhi(randomLegalPhi(4, 41001));
    const PureProjectorGreenResult direct = pureProjectorGreen(trial.Phi, trial.Phi);
    const PureProjectorGreenResult stable = pureProjectorGreenThinQr(trial.Phi, trial.Phi);
    require(direct.ok() && stable.ok(), "identity Green construction failed");
    const double directError = relativeError(direct.green, trial.G_T);
    const double stableError = relativeError(stable.green, trial.G_T);
    require(directError < kTolerance && stableError < kTolerance,
            "identity propagation did not reproduce trial Green");
    printResult("identity_propagation",
                {{"direct_error", directError}, {"thin_qr_error", stableError},
                 {"overlap_rcond", stable.overlap_rcond},
                 {"solve_residual", stable.solve_residual}});
}

void testRandomGaussian() {
    double maxOrthonormal = 0.0, maxIsotropic = 0.0, maxSkew = 0.0, maxInvolution = 0.0;
    for (int modes : {2, 4, 7}) {
        const GaussianTrialState trial =
            GaussianTrialState::fromPhi(randomLegalPhi(modes, 42000 + modes));
        require(trial.diagnostics.valid, "random Gaussian diagnostics failed");
        maxOrthonormal = std::max(maxOrthonormal, trial.diagnostics.orthonormal_error);
        maxIsotropic = std::max(maxIsotropic, trial.diagnostics.isotropy_error);
        maxSkew = std::max(maxSkew, trial.diagnostics.green_skew_error);
        maxInvolution = std::max(maxInvolution, trial.diagnostics.green_involution_error);
    }
    require(std::max({maxOrthonormal, maxIsotropic, maxSkew, maxInvolution}) < kTolerance,
            "random Gaussian invariant tolerance exceeded");
    printResult("random_gaussian",
                {{"max_orthonormal_error", maxOrthonormal},
                 {"max_isotropy_error", maxIsotropic},
                 {"max_green_skew_error", maxSkew},
                 {"max_green_involution_error", maxInvolution}});
}

void testFiniteLambda() {
    const GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(3));
    const MatType projector = trial.Phi * trial.Phi.adjoint();
    const MatType identity = MatType::Identity(projector.rows(), projector.cols());
    double previous = std::numeric_limits<double>::infinity();
    double finalError = previous;
    for (double lambda : {4.0, 8.0, 12.0, 24.0}) {
        const MatType boundary = std::exp(lambda) * projector +
                                 std::exp(-lambda) * (identity - projector);
        Eigen::FullPivLU<MatType> factor(identity + boundary);
        require(factor.isInvertible(), "finite-Lambda boundary solve is singular");
        const MatType finiteGreen = factor.solve(2.0 * identity) - identity;
        finalError = relativeError(finiteGreen, trial.G_T);
        require(finalError < previous, "finite-Lambda convergence is not monotone");
        previous = finalError;
    }
    require(finalError < 1e-9, "finite-Lambda boundary did not converge sufficiently");
    printResult("finite_lambda", {{"lambda_max", 24.0}, {"max_error", finalError}});
}

void testKitaevZeroModeParity() {
    double maxInvariant = 0.0;
    int ambiguityCount = 0;
    int oppositeParityCount = 0;
    for (int length : {4, 6}) {
        const SpinlessTvChainUtils config(length, 0.1, 0.0, 2, 1, 1.0, 0.0, 0);
        try {
            (void)GaussianTrialState::fromKitaevChain(config, 0.0);
        } catch (const GaussianTrialStateError &error) {
            if (error.code() == GaussianTrialStateStatus::zero_mode_ambiguous) ++ambiguityCount;
        }
        const GaussianTrialState plus = GaussianTrialState::fromKitaevChain(config, 1e-8);
        const GaussianTrialState minus = GaussianTrialState::fromKitaevChain(config, -1e-8);
        require(plus.diagnostics.valid && minus.diagnostics.valid,
                "split Kitaev trial state failed diagnostics");
        if (plus.fermionParity() == -minus.fermionParity()) ++oppositeParityCount;
        maxInvariant = std::max(
            maxInvariant,
            std::max({plus.diagnostics.orthonormal_error, plus.diagnostics.isotropy_error,
                      plus.diagnostics.green_skew_error,
                      plus.diagnostics.green_involution_error,
                      minus.diagnostics.orthonormal_error, minus.diagnostics.isotropy_error,
                      minus.diagnostics.green_skew_error,
                      minus.diagnostics.green_involution_error}));
    }
    require(ambiguityCount == 2, "unsplit Kitaev zero modes were not rejected");
    require(oppositeParityCount == 2, "edge splitting did not select opposite parity");
    require(maxInvariant < kTolerance, "Kitaev trial-state invariant tolerance exceeded");
    printResult("kitaev_zero_mode_parity",
                {{"ambiguity_count", double(ambiguityCount)},
                 {"opposite_parity_count", double(oppositeParityCount)},
                 {"max_invariant_error", maxInvariant}},
                {{"zero_mode_policy", "reject_without_explicit_splitting"}});
}

void testPropagationThinQr() {
    const int modes = 6;
    const int dimension = 2 * modes;
    const GaussianTrialState trial =
        GaussianTrialState::fromPhi(randomLegalPhi(modes, 44006));
    std::vector<MatType> ketFactors, braFactors;
    for (int step = 0; step < 24; ++step) {
        ketFactors.push_back(moderateMajoranaFactor(
            dimension, step % dimension, (step * 5 + 3) % dimension,
            0.025 + 0.002 * (step % 5)));
        braFactors.push_back(moderateMajoranaFactor(
            dimension, (step * 7 + 1) % dimension, (step * 11 + 4) % dimension,
            -0.02 - 0.001 * (step % 7)));
    }

    MatType directRight = trial.Phi, directLeft = trial.Phi;
    MatType qrRight = trial.Phi, qrLeft = trial.Phi;
    for (int step = 0; step < 24; ++step) {
        require(propagateKet(ketFactors[step], directRight), "direct ket propagation failed");
        require(propagateBra(braFactors[step], directLeft), "direct bra propagation failed");
        require(propagateKet(ketFactors[step], qrRight), "QR ket propagation failed");
        require(propagateBra(braFactors[step], qrLeft), "QR bra propagation failed");
        if ((step + 1) % 4 == 0) {
            const ThinQrResult rightQr = thinQrSubspace(qrRight);
            const ThinQrResult leftQr = thinQrSubspace(qrLeft);
            require(rightQr.ok() && leftQr.ok(), "periodic thin QR failed");
            qrRight = rightQr.q;
            qrLeft = leftQr.q;
        }
    }

    MatrixTestOperator ketOperator(ketFactors.front());
    MatrixTestOperator braOperator(braFactors.front());
    MatType matrixKet = trial.Phi, operatorKet = trial.Phi;
    MatType matrixBra = trial.Phi, operatorBra = trial.Phi;
    require(propagateKet(ketFactors.front(), matrixKet) &&
            propagateKet(ketOperator, operatorKet), "Operator ket helper failed");
    require(propagateBra(braFactors.front(), matrixBra) &&
            propagateBra(braOperator, operatorBra), "Operator bra helper failed");
    require(relativeError(matrixKet, operatorKet) < kTolerance &&
            relativeError(matrixBra, operatorBra) < kTolerance,
            "matrix and Operator propagation helpers disagree");

    const PureProjectorGreenResult direct = pureProjectorGreen(directRight, directLeft);
    const PureProjectorGreenResult rebuilt = pureProjectorGreenThinQr(qrRight, qrLeft);
    require(direct.ok() && rebuilt.ok(), "propagated Green construction failed");
    const double greenError = relativeError(direct.green, rebuilt.green);
    require(greenError < kTolerance, "direct and periodic thin-QR Green disagree");
    printResult("propagation_thin_qr",
                {{"green_error", greenError},
                 {"direct_solve_residual", direct.solve_residual},
                 {"thin_qr_solve_residual", rebuilt.solve_residual},
                 {"thin_qr_green_residual", rebuilt.green_residual}});
}

}  // namespace

int main(int argc, char **argv) try {
    if (argc != 2) {
        throw std::runtime_error(
            "usage: phase1_core_test single_site|identity_propagation|random_gaussian|"
            "finite_lambda|kitaev_zero_mode_parity|propagation_thin_qr");
    }
    const std::string test = argv[1];
    if (test == "single_site") testSingleSite();
    else if (test == "identity_propagation") testIdentityPropagation();
    else if (test == "random_gaussian") testRandomGaussian();
    else if (test == "finite_lambda") testFiniteLambda();
    else if (test == "kitaev_zero_mode_parity") testKitaevZeroModeParity();
    else if (test == "propagation_thin_qr") testPropagationThinQr();
    else throw std::runtime_error("unknown test: " + test);
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
