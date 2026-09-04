#include "stage1_dense.h"

#include <fstream>
#include <iostream>

using namespace driven_stage1;

namespace {

struct DenseOperators {
    std::vector<MatType> gamma;
    std::vector<MatType> density;
    MatType kinetic;
    MatType evenInteraction;
    MatType oddInteraction;
    MatType parity;
};

DenseOperators operators(const ModelParameters &p) {
    DenseOperators result;
    result.gamma = gammaMatrices(p.L);
    result.kinetic = -denseQuadratic(
        kineticGenerator(p, p.t, p.delta, p.mu), result.gamma);
    result.density.resize(p.L);
    for (int site = 0; site < p.L; ++site)
        result.density[site] = DataType(0.0, -0.5) *
            result.gamma[site] * result.gamma[p.L + site];
    result.evenInteraction = MatType::Zero(1 << p.L, 1 << p.L);
    result.oddInteraction = MatType::Zero(1 << p.L, 1 << p.L);
    const int bonds = p.boundary == 0 ? p.L : p.L - 1;
    for (int bond = 0; bond < bonds; ++bond) {
        const MatType term = result.density[bond] *
                             result.density[(bond + 1) % p.L];
        if ((bond & 1) == 0) result.evenInteraction += term;
        else result.oddInteraction += term;
    }
    result.parity = MatType::Zero(1 << p.L, 1 << p.L);
    for (int state = 0; state < (1 << p.L); ++state)
        result.parity(state, state) =
            (__builtin_popcount(unsigned(state)) & 1) ? -1.0 : 1.0;
    return result;
}

struct PropagatedState {
    cVecType state;
    double log_norm_squared = 0.0;
};

void normalizeTracking(PropagatedState &state) {
    const double norm = state.state.norm();
    require(std::isfinite(norm) && norm > 0.0,
            "ED propagation norm failed");
    state.log_norm_squared += 2.0 * std::log(norm);
    state.state /= norm;
}

PropagatedState sameContour(
        const DenseOperators &dense, const PureImaginaryTimeProtocol &protocol,
        const cVecType &initial) {
    MatType contour = MatType::Identity(initial.size(), initial.size());
    const MatType halfKinetic = exponential(
        dense.kinetic, -0.5 * protocol.deltaTau());
    for (double interaction : protocol.midpointValues()) {
        const MatType slice = halfKinetic * exponential(
            dense.evenInteraction, -protocol.deltaTau() * interaction) *
            exponential(
                dense.oddInteraction, -protocol.deltaTau() * interaction) *
            halfKinetic;
        contour = slice * contour;
    }
    PropagatedState result{contour * initial, 0.0};
    normalizeTracking(result);
    return result;
}

double protocolValue(const PureImaginaryTimeProtocol &protocol, double tau) {
    if (protocol.kind() == PureImaginaryTimeProtocolKind::LinearRamp)
        return protocol.initialValue() + protocol.rate() * tau;
    if (protocol.kind() == PureImaginaryTimeProtocolKind::SuddenQuench)
        return protocol.finalValue();
    return protocol.initialValue();
}

PropagatedState continuousReference(
        const DenseOperators &dense, const PureImaginaryTimeProtocol &protocol,
        const cVecType &initial, double referenceDt) {
    const int slices = int(std::llround(protocol.tauFinal() / referenceDt));
    require(slices > 0 &&
                std::abs(slices * referenceDt - protocol.tauFinal()) < 1e-12,
            "continuous-reference step does not divide tau_f");
    PropagatedState result{initial, 0.0};
    for (int slice = 0; slice < slices; ++slice) {
        const double tauMidpoint = (slice + 0.5) * referenceDt;
        const double interaction = protocolValue(protocol, tauMidpoint);
        const MatType hamiltonian = dense.kinetic + interaction *
            (dense.evenInteraction + dense.oddInteraction);
        result.state = exponential(hamiltonian, -referenceDt) * result.state;
        normalizeTracking(result);
    }
    return result;
}

struct DenseMeasurement {
    double energy = 0.0;
    double S_pi = 0.0;
    double S_pi_dq = 0.0;
    double R_CDW = 0.0;
    double parity = 0.0;
    MatType green;
};

DataType expectation(const cVecType &state, const MatType &op) {
    return (state.adjoint() * op * state)(0);
}

DenseMeasurement measure(const DenseOperators &dense, const cVecType &state,
                         double interaction, int length) {
    DenseMeasurement result;
    result.energy = expectation(
        state, dense.kinetic + interaction *
            (dense.evenInteraction + dense.oddInteraction)).real();
    auto structure = [&](double momentum) {
        MatType rho = MatType::Zero(state.size(), state.size());
        for (int site = 0; site < length; ++site)
            rho += std::exp(DataType(0.0, momentum * site)) * dense.density[site];
        return expectation(state, rho.adjoint() * rho).real() /
               double(length * length);
    };
    result.S_pi = structure(kPi);
    result.S_pi_dq = structure(kPi - 2.0 * kPi / length);
    result.R_CDW = 1.0 - result.S_pi_dq / result.S_pi;
    result.parity = expectation(state, dense.parity).real();
    result.green = MatType::Zero(2 * length, 2 * length);
    for (int row = 0; row < 2 * length; ++row)
        for (int column = 0; column < 2 * length; ++column)
            if (row != column)
                result.green(row, column) = -expectation(
                    state, dense.gamma[row] * dense.gamma[column]);
    return result;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        if (argc != 3) throw std::invalid_argument(
            "usage: driven_ed same_contour_reference.csv timestep.csv");
        ModelParameters p;
        p.L = 4;
        p.boundary = 0;
        p.hs = 0;
        const GaussianTrialState gaussianInitial = makeInitialState(p, 0.0);
        const int initialParity = gaussianInitial.fermionParity();
        const DenseOperators dense = operators(p);
        const DenseInitialState initial = denseInitialState(
            p, dense.gamma, gaussianInitial, 0.0, initialParity);
        const double evenInteractionActionNorm =
            (dense.evenInteraction * initial.vector).norm();
        const double oddInteractionActionNorm =
            (dense.oddInteraction * initial.vector).norm();
        const double fullInteractionActionNorm =
            ((dense.evenInteraction + dense.oddInteraction) *
             initial.vector).norm();

        std::ofstream csv(argv[1]);
        require(bool(csv), "cannot open driven ED CSV");
        csv << "protocol,tau_f,Delta_tau,V_i,V_f,R,Phi_0_parity,"
               "contour_log_norm_squared,contour_denominator,"
               "hs_reconstructed_denominator,hs_denominator_relative_error,"
               "hs_state_max_abs_error,contour_energy,contour_S_pi,"
               "contour_R_CDW,contour_parity,contour_G_0_1_real,"
               "contour_G_0_1_imag,contour_G_0_L_real,contour_G_0_L_imag,"
               "reference_dt,reference_energy,reference_S_pi,reference_R_CDW,"
               "reference_parity,reference_G_0_1_real,reference_G_0_1_imag,"
               "reference_G_0_L_real,reference_G_0_L_imag,energy_abs_error,"
               "S_pi_abs_error,R_CDW_abs_error,green_max_abs_error\n";

        double maxTrotterEnergy = 0.0, maxTrotterGreen = 0.0;
        for (const std::string kind : {"quench", "ramp"}) {
            for (double tau : {0.1, 0.2, 0.4}) {
                for (double dt : (kind == "ramp" && tau == 0.4 ?
                                  std::vector<double>{0.1, 0.05} :
                                  std::vector<double>{0.1})) {
                    p.dt = dt;
                    const PureImaginaryTimeProtocol protocol = kind == "quench" ?
                        PureImaginaryTimeProtocol::suddenQuench(
                            0.0, 2.0, tau, dt) :
                        PureImaginaryTimeProtocol::linearRamp(
                            0.0, 10.0, tau, dt);
                    const PropagatedState contour = sameContour(
                        dense, protocol, initial.vector);
                    double hsDenominator = 0.0;
                    const cVecType hsState = denseHsSameContourState(
                        p, protocol, dense.gamma, initial.vector,
                        &hsDenominator);
                    const double contourDenominator =
                        std::exp(contour.log_norm_squared);
                    const double hsDenominatorError = std::abs(
                        hsDenominator - contourDenominator) /
                        contourDenominator;
                    const double hsStateError =
                        (hsState - contour.state).cwiseAbs().maxCoeff();
                    require(hsDenominatorError < 1e-11 &&
                                hsStateError < 1e-11,
                            "HS reconstruction differs from same-contour ED");
                    const double referenceDt = 0.00625;
                    const PropagatedState reference = continuousReference(
                        dense, protocol, initial.vector, referenceDt);
                    const DenseMeasurement contourMeasurement = measure(
                        dense, contour.state, protocol.finalValue(), p.L);
                    const DenseMeasurement referenceMeasurement = measure(
                        dense, reference.state, protocol.finalValue(), p.L);
                    const double energyError = std::abs(
                        contourMeasurement.energy - referenceMeasurement.energy);
                    const double spiError = std::abs(
                        contourMeasurement.S_pi - referenceMeasurement.S_pi);
                    const double rError = std::abs(
                        contourMeasurement.R_CDW - referenceMeasurement.R_CDW);
                    const double greenError =
                        (contourMeasurement.green - referenceMeasurement.green)
                            .cwiseAbs().maxCoeff();
                    maxTrotterEnergy = std::max(maxTrotterEnergy, energyError);
                    maxTrotterGreen = std::max(maxTrotterGreen, greenError);
                    csv << protocol.name() << ',' << std::setprecision(17)
                        << tau << ',' << dt << ',' << protocol.initialValue() << ','
                        << protocol.finalValue() << ',' << protocol.rate() << ','
                        << initialParity << ',' << contour.log_norm_squared << ','
                        << contourDenominator << ',' << hsDenominator << ','
                        << hsDenominatorError << ',' << hsStateError << ','
                        << contourMeasurement.energy << ',' << contourMeasurement.S_pi
                        << ',' << contourMeasurement.R_CDW << ','
                        << contourMeasurement.parity << ','
                        << contourMeasurement.green(0, 1).real() << ','
                        << contourMeasurement.green(0, 1).imag() << ','
                        << contourMeasurement.green(0, p.L).real() << ','
                        << contourMeasurement.green(0, p.L).imag() << ','
                        << referenceDt << ',' << referenceMeasurement.energy << ','
                        << referenceMeasurement.S_pi << ','
                        << referenceMeasurement.R_CDW << ','
                        << referenceMeasurement.parity << ','
                        << referenceMeasurement.green(0, 1).real() << ','
                        << referenceMeasurement.green(0, 1).imag() << ','
                        << referenceMeasurement.green(0, p.L).real() << ','
                        << referenceMeasurement.green(0, p.L).imag() << ','
                        << energyError << ',' << spiError << ',' << rError << ','
                        << greenError << '\n';
                }
            }
        }
        csv.flush();
        require(bool(csv), "driven ED CSV write failed");

        // The required odd-parity L=4 sector is annihilated by the full
        // nearest-neighbour interaction.  Retain it above as the physical
        // benchmark.  For a separately labelled timestep-only probe, occupy
        // the lowest positive BdG mode instead of its negative partner.  This
        // explicitly selects even parity and makes V(tau) nontrivial.
        const double probeMu = 0.0;
        const double probeDelta = 1.0;
        const MatType probeHamiltonian = kineticGenerator(
            p, 1.0, probeDelta, probeMu);
        Eigen::SelfAdjointEigenSolver<MatType> probeSolver(probeHamiltonian);
        require(probeSolver.info() == Eigen::Success,
                "timestep probe eigensolver failed");
        MatType probePhi = probeSolver.eigenvectors().leftCols(p.L);
        probePhi.col(p.L - 1) = probePhi.col(p.L - 1).conjugate();
        const GaussianTrialState gaussianProbe =
            GaussianTrialState::fromPhi(probePhi);
        require(gaussianProbe.fermionParity() == -initialParity,
                "timestep probe did not flip the explicit parity sector");
        const DenseInitialState probe = denseInitialState(
            p, dense.gamma, gaussianProbe, probeMu,
            gaussianProbe.fermionParity(), 1.0, probeDelta, false);
        struct TimestepRow {
            double dt = 0.0;
            double energyError = 0.0;
            double spiError = 0.0;
            double rError = 0.0;
            double greenError = 0.0;
        };
        std::vector<TimestepRow> timestepRows;
        for (double dt : {0.1, 0.05}) {
            p.dt = dt;
            const PureImaginaryTimeProtocol protocol =
                PureImaginaryTimeProtocol::linearRamp(0.0, 10.0, 0.4, dt);
            const PropagatedState contour = sameContour(
                dense, protocol, probe.vector);
            const PropagatedState reference = continuousReference(
                dense, protocol, probe.vector, 0.003125);
            const DenseMeasurement cm = measure(
                dense, contour.state, protocol.finalValue(), p.L);
            const DenseMeasurement rm = measure(
                dense, reference.state, protocol.finalValue(), p.L);
            timestepRows.push_back({dt, std::abs(cm.energy - rm.energy),
                std::abs(cm.S_pi - rm.S_pi),
                std::abs(cm.R_CDW - rm.R_CDW),
                (cm.green - rm.green).cwiseAbs().maxCoeff()});
        }
        std::ofstream timestep(argv[2]);
        require(bool(timestep), "cannot open timestep convergence CSV");
        timestep << "L,boundary,protocol,V_i,R,tau_f,Phi_0_definition,"
                    "Phi_0_parity,Delta_tau,reference_dt,energy_abs_error,"
                    "S_pi_abs_error,R_CDW_abs_error,green_max_abs_error,"
                    "energy_error_reduction_factor,green_error_reduction_factor\n";
        for (std::size_t index = 0; index < timestepRows.size(); ++index) {
            const TimestepRow &row = timestepRows[index];
            const bool fine = index > 0;
            timestep << "4,pbc,linear_ramp,0,10,0.4,"
                        "one_BdG_quasiparticle_even_parity_timestep_probe,"
                     << gaussianProbe.fermionParity() << ','
                     << std::setprecision(17)
                     << row.dt << ",0.003125," << row.energyError << ','
                     << row.spiError << ',' << row.rError << ','
                     << row.greenError << ',';
            if (fine)
                timestep << timestepRows[0].energyError / row.energyError;
            timestep << ',';
            if (fine)
                timestep << timestepRows[0].greenError / row.greenError;
            timestep << '\n';
        }
        timestep.flush();
        require(bool(timestep), "timestep convergence CSV write failed");
        std::cout << std::setprecision(17)
                  << "{\"status\":\"complete\",\"L\":4,"
                  << "\"boundary\":\"pbc\",\"Phi_0_parity\":"
                  << initialParity
                  << ",\"same_contour_order\":\"Khalf_even_odd_Khalf\","
                  << "\"continuous_reference_dt\":0.00625,"
                  << "\"even_interaction_action_norm_on_Phi_0\":"
                  << evenInteractionActionNorm << ','
                  << "\"odd_interaction_action_norm_on_Phi_0\":"
                  << oddInteractionActionNorm << ','
                  << "\"full_interaction_action_norm_on_Phi_0\":"
                  << fullInteractionActionNorm << ','
                  << "\"probe_energy_error_reduction_factor\":"
                  << timestepRows[0].energyError /
                         timestepRows[1].energyError << ','
                  << "\"probe_green_error_reduction_factor\":"
                  << timestepRows[0].greenError /
                         timestepRows[1].greenError << ','
                  << "\"maximum_dt_0p1_trotter_energy_error\":"
                  << maxTrotterEnergy
                  << ",\"maximum_dt_0p1_trotter_green_error\":"
                  << maxTrotterGreen << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "driven_ed: " << error.what() << '\n';
        return 1;
    }
}
