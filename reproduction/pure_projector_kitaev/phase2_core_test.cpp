#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "gaussian_trial_state.h"
#include "pure_projector_static.h"
#include "pure_projector_weight.h"

namespace {

void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}

double relativeError(const MatType &a, const MatType &b) {
    return (a - b).norm() / std::max(1.0, b.norm());
}

MatType canonicalPhi(int modes) {
    MatType phi = MatType::Zero(2 * modes, modes);
    for (int mode = 0; mode < modes; ++mode) {
        phi(2 * mode, mode) = 1.0 / std::sqrt(2.0);
        phi(2 * mode + 1, mode) = DataType(0.0, 1.0) / std::sqrt(2.0);
    }
    return phi;
}

MatType randomUnitary(int n, unsigned seed) {
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    MatType raw(n, n);
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < n; ++col)
            raw(row, col) = DataType(normal(generator), normal(generator));
    Eigen::HouseholderQR<MatType> qr(raw);
    return qr.householderQ() * MatType::Identity(n, n);
}

PureProjectorSlice rotationSlice(int dimension, int a, int b, double x,
                                 const std::string &label) {
    MatType h = MatType::Zero(dimension, dimension);
    h(a, b) = DataType(0.0, x);
    h(b, a) = -h(a, b);
    MatType exponent = h;
    return PureProjectorSlice(expm(exponent, 1.0), DataType(1.0), label);
}

void identityAndGaugeTest() {
    GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(3));
    PureProjectorWeightEvaluator evaluator;
    PureProjectorWeightResult identity = evaluator.evaluate(trial, {});
    require(identity.ok(), "identity evaluation failed");
    require(std::abs(identity.log_abs_weight) < 1e-14, "identity weight is not one");
    require(identity.z2_sign == 1, "identity Z2 is not +1");
    require(relativeError(identity.green, trial.G_T) < 1e-12,
            "identity Green differs from trial Green");

    std::vector<PureProjectorSlice> slices{
        rotationSlice(6, 0, 3, 0.07, "a"),
        rotationSlice(6, 1, 4, -0.04, "b")};
    PureProjectorWeightResult original = evaluator.evaluate(trial, slices);
    GaussianTrialState gauged = GaussianTrialState::fromPhi(
        trial.Phi * randomUnitary(trial.Phi.cols(), 9917));
    PureProjectorWeightResult transformed = evaluator.evaluate(gauged, slices);
    require(original.ok() && transformed.ok(), "gauge evaluation failed");
    require(std::abs(original.log_abs_weight - transformed.log_abs_weight) < 1e-11,
            "weight is not trial-gauge invariant");
    require(std::abs(original.complex_phase - transformed.complex_phase) < 1e-11,
            "phase is not trial-gauge invariant");
    require(original.z2_sign == transformed.z2_sign, "Z2 is not gauge invariant");
    require(relativeError(original.green, transformed.green) < 1e-11,
            "Green is not gauge invariant");
    std::cout << "PASS identity_contour\nPASS trial_state_gauge_invariance\n";
}

void contourOrderTest() {
    GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(2));
    PureStaticProjectorContour contour(trial, 0.2, 0.1, 1, 1e-8);
    PureProjectorSlice a = rotationSlice(4, 0, 2, 0.23, "ket_a");
    PureProjectorSlice b = rotationSlice(4, 1, 2, -0.17, "ket_b");
    contour.setKetActionOrder({a, b});
    contour.setBraProtocolOrder({a, b});
    const PureStaticMeasurement measurement = contour.measurementGreen();
    require(measurement.ok(), "ordered contour Green failed");

    MatType ket = b.matrix * a.matrix * trial.Phi;
    MatType bra = a.matrix.adjoint() * b.matrix.adjoint() * trial.Phi;
    PureProjectorGreenResult expected = pureProjectorGreenThinQr(ket, bra);
    require(expected.ok() && relativeError(measurement.green, expected.green) < 1e-11,
            "contour did not use true ket/strict-adjoint bra order");

    MatType wrongBra = b.matrix.adjoint() * a.matrix.adjoint() * trial.Phi;
    PureProjectorGreenResult wrong = pureProjectorGreenThinQr(ket, wrongBra);
    require(wrong.ok() && relativeError(measurement.green, wrong.green) > 1e-5,
            "noncommuting order test lacks sensitivity");
    require(contour.flattenedActionOrder().size() == 4,
            "flattened action order was not recorded");
    std::cout << "PASS noncommuting_contour_order\n";
}

void failurePolicyTest() {
    GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(2));
    PureProjectorWeightOptions options;
    options.mode = PureProjectorWeightMode::RealZ2;
    PureProjectorWeightEvaluator evaluator(options);
    PureProjectorSlice complexPhase(MatType::Identity(4, 4), DataType(0.0, 1.0),
                                    "complex_eta");
    PureProjectorWeightResult result = evaluator.evaluate(trial, {complexPhase});
    require(!result.ok() && result.status == PureProjectorWeightStatus::complex_ratio,
            "RealZ2 accepted a clearly complex ratio");
    require(result.first_failing_slice == 0, "first failing slice was not retained");

    PureProjectorSlice zero(MatType::Zero(4, 4), DataType(1.0), "zero");
    result = evaluator.evaluate(trial, {zero});
    require(result.status == PureProjectorWeightStatus::zero_weight && result.z2_sign == 0,
            "zero weight did not fail closed");
    require(!std::isnan(result.log_abs_weight), "zero weight produced NaN");
    std::cout << "PASS real_z2_and_zero_weight_fail_closed\n";
}

void slowWalkerTransactionTest() {
    GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(2));
    PureProjectorConfiguration initial;
    initial.slices.push_back(rotationSlice(4, 0, 2, 0.08, "hs0"));
    initial.hs_fields = {1};
    PureProjectorSlowWalker walker(trial, initial, 12345);
    const PureProjectorWalkerSnapshot before = walker.snapshot();
    PureProjectorProposal proposal{0, -1, 0.3141592653589793};
    const PureProjectorProposalResult outcome = walker.propose(proposal);
    require(outcome.uniform == proposal.uniform, "proposal uniform changed");
    require(outcome.proposal_index == proposal.index, "proposal identity changed");
    require(outcome.predecision_hash == before.state_hash,
            "walker mutated before Metropolis decision");
    require(walker.currentWeight().ok(), "walker installed an invalid state");
    std::cout << "PASS slow_walker_transaction\n";
}

} // namespace

int main() {
    try {
        identityAndGaugeTest();
        contourOrderTest();
        failurePolicyTest();
        slowWalkerTransactionTest();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
