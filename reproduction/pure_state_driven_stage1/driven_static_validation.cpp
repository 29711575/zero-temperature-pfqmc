#include "stage1_common.h"

#include <fstream>
#include <iostream>

using namespace driven_stage1;

namespace {

PureFastProposal staticFlip(
        const ModelParameters &p, double interaction,
        const PureFastConfiguration &configuration, int index, double uniform) {
    const PureSliceLocation &location = configuration.locations[index];
    const int newHs = -configuration.hs_fields[index];
    MatType factor = localHsFactor(
        p, interaction, location.bond, location.aux, newHs);
    if (location.branch == PureBranch::Bra) factor.adjointInPlace();
    PureFastProposal proposal;
    proposal.index = index;
    proposal.new_hs = newHs;
    proposal.new_factor = std::move(factor);
    proposal.uniform = uniform;
    return proposal;
}

double sliceError(const PureFastConfiguration &left,
                  const PureFastConfiguration &right) {
    require(left.slices.size() == right.slices.size(),
            "static/driven contour size mismatch");
    double error = 0.0;
    for (std::size_t index = 0; index < left.slices.size(); ++index)
        error = std::max(error,
            (left.slices[index].matrix - right.slices[index].matrix).norm());
    return error;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        if (argc != 2) throw std::invalid_argument(
            "usage: driven_static_validation output.csv");
        std::ofstream csv(argv[1]);
        require(bool(csv), "cannot open driven_static_limit.csv");
        csv << "L,event,proposal,weight_log_error,weight_phase_error,z2_mismatch,"
               "green_relative_error,ratio_error,accept_mismatch,configuration_hash_mismatch,"
               "rng_hash_mismatch,energy_error,S_pi_error,R_CDW_error,parity_error,status\n";

        double maxWeight = 0.0, maxPhase = 0.0, maxGreen = 0.0;
        double maxRatio = 0.0, maxObservable = 0.0;
        long long z2Mismatch = 0, acceptMismatch = 0;
        long long hashMismatch = 0, rngMismatch = 0;
        for (int length : {2, 4}) {
            ModelParameters p;
            p.L = length;
            p.boundary = 0;
            p.hs = 0;
            p.dt = 0.1;
            const double interaction = 2.0;
            const double tauFinal = 0.3;
            const int slices = 3;
            const std::uint64_t seed = 710000 + length;
            GaussianTrialState initial = makeInitialState(p, 0.3);
            const int initialParity = initial.fermionParity();
            require(initialParity == 1 || initialParity == -1,
                    "static-limit initial parity is not explicit");

            // Continue the very same engine after initialization, as the
            // production static driver does.  Separate copies only permit the
            // path-by-path comparison.
            std::mt19937_64 staticRng(seed), drivenRng(seed);
            PureDrivenContourResult staticContour = makeStaticIndependentContour(
                p, interaction, slices, staticRng);
            const PureImaginaryTimeProtocol protocol =
                PureImaginaryTimeProtocol::constant(interaction, tauFinal, p.dt);
            PureDrivenContourResult drivenContour =
                makeDrivenContour(p, protocol, drivenRng);
            require(staticContour.configuration.hs_fields ==
                        drivenContour.configuration.hs_fields,
                    "constant protocol changed initial HS stream");
            require(sliceError(staticContour.configuration,
                               drivenContour.configuration) < 1e-14,
                    "constant protocol changed initial factors");

            PureFastOptions options;
            options.weight_mode = PureProjectorWeightMode::RealZ2;
            PureProjectorFastWalker staticWalker(
                initial, std::move(staticContour.configuration), 4,
                PureFastRunMode::FastStrict, options,
                PureFastInitializationPolicy::SequentialAudit);
            PureProjectorFastWalker drivenWalker(
                initial, std::move(drivenContour.configuration), 4,
                PureFastRunMode::FastStrict, options,
                PureFastInitializationPolicy::SequentialAudit);
            std::vector<int> indices = proposalIndices(staticWalker.configuration());
            require(indices == proposalIndices(drivenWalker.configuration()),
                    "constant protocol changed proposal identity list");
            std::uniform_real_distribution<double> uniform(0.0, 1.0);

            auto compare = [&](const std::string &event, int proposal,
                               double ratioError, bool decisionMismatch) {
                const auto &sw = staticWalker.currentWeight();
                const auto &dw = drivenWalker.currentWeight();
                const double weightError = std::abs(
                    sw.log_abs_weight - dw.log_abs_weight);
                const double phaseError = std::abs(
                    sw.complex_phase - dw.complex_phase);
                const int signMismatch = staticWalker.z2Sign() != drivenWalker.z2Sign();
                const PureProjectorGreenResult sg = staticWalker.measurementGreen();
                const PureProjectorGreenResult dg = drivenWalker.measurementGreen();
                require(sg.ok() && dg.ok(), "static-limit Green reconstruction failed");
                const double greenError = relativeError(sg.green, dg.green);
                const Observables so = observe(p, interaction, sg.green);
                const Observables d = observe(p, interaction, dg.green);
                const double energyError = std::abs(so.energy - d.energy);
                const double spiError = std::abs(so.S_pi - d.S_pi);
                const double rError = std::abs(so.R_CDW - d.R_CDW);
                const double parityError = std::abs(so.parity - d.parity);
                const bool configurationMismatch =
                    staticWalker.configurationHash() != drivenWalker.configurationHash();
                const bool randomMismatch =
                    hashRng(staticRng) != hashRng(drivenRng);
                maxWeight = std::max(maxWeight, weightError);
                maxPhase = std::max(maxPhase, phaseError);
                maxGreen = std::max(maxGreen, greenError);
                maxRatio = std::max(maxRatio, ratioError);
                maxObservable = std::max({maxObservable, energyError, spiError,
                                          rError, parityError});
                z2Mismatch += signMismatch;
                acceptMismatch += decisionMismatch;
                hashMismatch += configurationMismatch;
                rngMismatch += randomMismatch;
                const bool pass = weightError < 1e-12 && phaseError < 1e-12 &&
                    !signMismatch && greenError < 1e-11 && ratioError < 1e-11 &&
                    !decisionMismatch && !configurationMismatch && !randomMismatch &&
                    std::max({energyError, spiError, rError, parityError}) < 1e-11;
                csv << length << ',' << event << ',' << proposal << ','
                    << std::setprecision(17) << weightError << ',' << phaseError << ','
                    << signMismatch << ',' << greenError << ',' << ratioError << ','
                    << decisionMismatch << ',' << configurationMismatch << ','
                    << randomMismatch << ',' << energyError << ',' << spiError << ','
                    << rError << ',' << parityError << ',' << (pass ? "PASS" : "FAIL")
                    << '\n';
                require(pass, "driven constant protocol differs from static projector");
            };
            compare("initial", -1, 0.0, false);
            for (int proposal = 0; proposal < 64; ++proposal) {
                const int index = indices[proposal % indices.size()];
                const double us = uniform(staticRng);
                const double ud = uniform(drivenRng);
                require(us == ud, "constant protocol changed proposal uniform");
                const PureFastProposalResult sr = staticWalker.propose(
                    staticFlip(p, interaction, staticWalker.configuration(), index, us));
                const PureFastProposalResult dr = drivenWalker.propose(
                    flipProposal(p, protocol, drivenWalker.configuration(), index, ud));
                require(!sr.terminated && !dr.terminated &&
                            sr.ratio.ok() && dr.ratio.ok(),
                        "static-limit proposal failed closed");
                const double ratioError = std::abs(sr.ratio.ratio - dr.ratio.ratio) /
                    std::max(1e-14, std::abs(sr.ratio.ratio));
                compare("proposal", proposal, ratioError,
                        sr.accepted != dr.accepted);
            }
        }
        csv.flush();
        require(bool(csv), "driven_static_limit.csv write failed");
        std::cout << std::setprecision(17)
                  << "{\"status\":\"PASS\",\"max_weight_log_error\":"
                  << maxWeight << ",\"max_weight_phase_error\":" << maxPhase
                  << ",\"max_green_relative_error\":" << maxGreen
                  << ",\"max_ratio_error\":" << maxRatio
                  << ",\"max_observable_error\":" << maxObservable
                  << ",\"z2_mismatch\":" << z2Mismatch
                  << ",\"accept_mismatch\":" << acceptMismatch
                  << ",\"configuration_hash_mismatch\":" << hashMismatch
                  << ",\"rng_hash_mismatch\":" << rngMismatch << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "driven_static_validation: " << error.what() << '\n';
        return 1;
    }
}
