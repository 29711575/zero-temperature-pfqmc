#ifndef PURE_STATE_DRIVEN_STAGE1_COMMON_H
#define PURE_STATE_DRIVEN_STAGE1_COMMON_H

#include "kitaevChain.h"
#include "pure_projector_driven.h"
#include "pure_projector_observables.h"
#include "pure_projector_protocol.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace driven_stage1 {

constexpr double kPi = 3.141592653589793238462643383279502884;

inline void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

struct ModelParameters {
    int L = 4;
    int boundary = 0;
    int hs = 0;
    double dt = 0.1;
    double t = 1.0;
    double delta = 1.0;
    double mu = 0.0;
};

inline MatType kineticGenerator(const ModelParameters &p, double t,
                                double delta, double mu) {
    auto component = [&](double d, double m) {
        SpinlessTvChainUtils config(
            p.L, 1.0, 0.0, 2, p.boundary, d, m, p.hs);
        MatType generator = MatType::Zero(2 * p.L, 2 * p.L);
        config.KineticGenerator(generator);
        return generator;
    };
    const MatType base = component(0.0, 0.0);
    return t * base + delta * (component(1.0, 0.0) - base) +
           mu * (component(0.0, 1.0) - base);
}

inline MatType exponential(MatType generator, double scale = 1.0) {
    return expm(generator, scale);
}

inline MatType localHsGenerator(const ModelParameters &p, double interaction,
                                int layer, int auxiliary, int sigma) {
    SpinlessTvChainUtils config(
        p.L, p.dt, interaction, 2, p.boundary, p.delta, p.mu, p.hs);
    MatType generator = MatType::Zero(2 * p.L, 2 * p.L);
    const double lambda = std::acosh(std::exp(0.5 * interaction * p.dt));
    int a, b, c, d;
    config.aux2MajoranaIdx(auxiliary, 0, layer, a, b);
    config.aux2MajoranaIdx(auxiliary, 1, layer, c, d);
    const DataType value(0.0, lambda * sigma);
    if (p.hs == 0) {
        generator(a, b) = value;
        generator(b, a) = -value;
        generator(c, d) = value;
        generator(d, c) = -value;
    } else {
        generator(a, c) = value;
        generator(c, a) = -value;
        generator(b, d) = -value;
        generator(d, b) = value;
    }
    return generator;
}

inline MatType localHsFactor(const ModelParameters &p, double interaction,
                             int layer, int auxiliary, int sigma) {
    return exponential(localHsGenerator(
        p, interaction, layer, auxiliary, sigma));
}

inline GaussianTrialState makeInitialState(
        const ModelParameters &p, double trialMu = 0.0,
        double edgeSplitting = 0.0, int expectedParity = 0) {
    MatType h = kineticGenerator(p, 1.0, 1.0, trialMu);
    if (edgeSplitting != 0.0) {
        SpinlessTvChainUtils coordinates(
            p.L, p.dt, 0.0, 2, p.boundary, 1.0, trialMu, p.hs);
        const int left = coordinates.majoranaCoord2Idx(0, 1);
        const int right = coordinates.majoranaCoord2Idx(p.L - 1, 1);
        const DataType value(0.0, edgeSplitting);
        h(left, right) += value;
        h(right, left) -= value;
    }
    GaussianTrialState state = GaussianTrialState::fromMajoranaHamiltonian(h);
    if (expectedParity != 0 && state.fermionParity() != expectedParity)
        throw std::invalid_argument(
            "Phi_0 parity does not match the explicit parity policy");
    return state;
}

inline PureFastConfiguration makeProtocolBranch(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        PureBranch branch, std::mt19937_64 &rng) {
    require(protocol.slices() > 0 &&
                std::abs(protocol.deltaTau() - p.dt) < 1e-12,
            "protocol Delta_tau does not match model dt");
    const MatType kinetic = kineticGenerator(p, p.t, p.delta, p.mu);
    const MatType halfKinetic = exponential(kinetic, -0.5 * p.dt);
    const auto counts = pureProjectorCheckerboardBondCounts(p.L, p.boundary);
    std::uniform_int_distribution<int> bit(0, 1);
    PureFastConfiguration result;
    int factor = 0;
    auto push = [&](const MatType &matrix, int field, int slice, int layer,
                    int auxiliary, const std::string &label) {
        result.slices.emplace_back(matrix, 1.0, label);
        result.hs_fields.push_back(field);
        result.locations.push_back(
            {branch, slice, factor++, layer, auxiliary});
    };
    // The validated Majorana-to-Fock representation is order reversing:
    // flattened one-particle factors B_0,...,B_n represent the dense Fock
    // product U_0...U_n.  Store late imaginary-time slices first so that the
    // represented ket is the time-ordered U(tau_f,0)=U_{n}...U_0.
    for (int slice = protocol.slices() - 1; slice >= 0; --slice) {
        const double interaction = protocol.midpointValue(slice);
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:first");
        for (int layer = 0; layer < 2; ++layer) {
            const int count = layer == 0 ? counts.first : counts.second;
            for (int auxiliary = 0; auxiliary < count; ++auxiliary) {
                const int sigma = bit(rng) ? 1 : -1;
                push(localHsFactor(p, interaction, layer, auxiliary, sigma),
                     sigma, slice, layer, auxiliary,
                     "tau:" + std::to_string(slice) + ":V" +
                         std::to_string(layer) + ":" +
                         std::to_string(auxiliary));
            }
        }
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:last");
    }
    return result;
}

inline PureFastConfiguration makeStaticBranch(
        const ModelParameters &p, double interaction, int slices,
        PureBranch branch, std::mt19937_64 &rng) {
    require(slices > 0, "static branch requires positive slice count");
    const MatType kinetic = kineticGenerator(p, p.t, p.delta, p.mu);
    const MatType halfKinetic = exponential(kinetic, -0.5 * p.dt);
    const auto counts = pureProjectorCheckerboardBondCounts(p.L, p.boundary);
    std::uniform_int_distribution<int> bit(0, 1);
    PureFastConfiguration result;
    int factor = 0;
    auto push = [&](const MatType &matrix, int field, int slice, int layer,
                    int auxiliary, const std::string &label) {
        result.slices.emplace_back(matrix, 1.0, label);
        result.hs_fields.push_back(field);
        result.locations.push_back(
            {branch, slice, factor++, layer, auxiliary});
    };
    for (int slice = 0; slice < slices; ++slice) {
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:first");
        for (int layer = 0; layer < 2; ++layer) {
            const int count = layer == 0 ? counts.first : counts.second;
            for (int auxiliary = 0; auxiliary < count; ++auxiliary) {
                const int sigma = bit(rng) ? 1 : -1;
                push(localHsFactor(p, interaction, layer, auxiliary, sigma),
                     sigma, slice, layer, auxiliary,
                     "tau:" + std::to_string(slice) + ":V" +
                         std::to_string(layer) + ":" +
                         std::to_string(auxiliary));
            }
        }
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:last");
    }
    return result;
}

inline PureDrivenContourResult makeDrivenContour(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        std::mt19937_64 &rng) {
    PureFastConfiguration ket =
        makeProtocolBranch(p, protocol, PureBranch::Ket, rng);
    PureFastConfiguration bra =
        makeProtocolBranch(p, protocol, PureBranch::Bra, rng);
    return pureProjectorIndependentBraKetContour(
        std::move(ket), std::move(bra));
}

inline PureDrivenContourResult makeStaticIndependentContour(
        const ModelParameters &p, double interaction, int slices,
        std::mt19937_64 &rng) {
    PureFastConfiguration ket =
        makeStaticBranch(p, interaction, slices, PureBranch::Ket, rng);
    PureFastConfiguration bra =
        makeStaticBranch(p, interaction, slices, PureBranch::Bra, rng);
    // Exact compatibility copy of the retained production sequential-audit
    // initializer.  Production appends the reverse bra protocol without an
    // explicit adjoint because every current static K/2 and HS factor is
    // Hermitian.  The driven builder independently enforces the general
    // reverse-adjoint rule, so the static-limit test compares the two paths.
    PureDrivenContourResult result;
    result.ket_factor_count = int(ket.slices.size());
    result.bra_factor_count = int(bra.slices.size());
    result.configuration = std::move(ket);
    for (int index = int(bra.slices.size()) - 1; index >= 0; --index) {
        result.configuration.slices.push_back(bra.slices[index]);
        result.configuration.hs_fields.push_back(bra.hs_fields[index]);
        result.configuration.locations.push_back(bra.locations[index]);
    }
    return result;
}

inline std::vector<int> proposalIndices(const PureFastConfiguration &configuration) {
    std::vector<int> result;
    for (int index = 0; index < int(configuration.locations.size()); ++index)
        if (configuration.locations[index].aux >= 0) result.push_back(index);
    return result;
}

inline PureFastProposal flipProposal(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        const PureFastConfiguration &configuration, int index, double uniform) {
    const PureSliceLocation &location = configuration.locations[index];
    const int newHs = -configuration.hs_fields[index];
    MatType factor = localHsFactor(
        p, protocol.midpointValue(location.slice), location.bond,
        location.aux, newHs);
    if (location.branch == PureBranch::Bra) factor.adjointInPlace();
    PureFastProposal proposal;
    proposal.index = index;
    proposal.new_hs = newHs;
    proposal.new_factor = std::move(factor);
    proposal.new_eta = 1.0;
    proposal.uniform = uniform;
    return proposal;
}

inline DataType energy(const ModelParameters &p, double interaction,
                       const MatType &green) {
    auto value = [&](double V, double delta, double mu) {
        SpinlessTvChainUtils config(
            p.L, p.dt, V, 2, p.boundary, delta, mu, p.hs);
        return config.energyFromGreensFunc(green);
    };
    const DataType base = value(0.0, 0.0, 0.0);
    return p.t * base + p.delta * (value(0.0, 1.0, 0.0) - base) +
           p.mu * (value(0.0, 0.0, 1.0) - base) +
           (value(interaction, 0.0, 0.0) - base);
}

struct Observables {
    double energy = 0.0;
    double S_pi = 0.0;
    double S_pi_dq = 0.0;
    double R_CDW = 0.0;
    double parity = 0.0;
};

inline Observables observe(const ModelParameters &p, double interaction,
                           const MatType &green) {
    Observables result;
    result.energy = energy(p, interaction, green).real();
    result.S_pi = pureProjectorStructureFactor(green, p.L, kPi).real();
    result.S_pi_dq = pureProjectorStructureFactor(
        green, p.L, kPi - 2.0 * kPi / p.L).real();
    result.R_CDW = 1.0 - result.S_pi_dq / result.S_pi;
    const PurePhysicalParityResult parity = pureProjectorPhysicalParity(green);
    require(parity.ok(), "fermion parity measurement is untrusted");
    result.parity = parity.physical_parity;
    return result;
}

inline std::uint64_t hashRng(const std::mt19937_64 &rng) {
    std::ostringstream state;
    state << rng;
    std::uint64_t hash = 1469598103934665603ULL;
    for (unsigned char byte : state.str()) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    return hash;
}

inline double relativeError(const MatType &left, const MatType &right) {
    return (left - right).norm() / std::max(1.0, right.norm());
}

}  // namespace driven_stage1

#endif
