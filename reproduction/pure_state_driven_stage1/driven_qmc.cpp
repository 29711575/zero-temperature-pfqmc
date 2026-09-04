#include "stage1_common.h"

#include <fstream>
#include <iostream>

using namespace driven_stage1;

namespace {

std::map<std::string, std::string> arguments(int argc, char **argv) {
    std::map<std::string, std::string> result;
    for (int index = 1; index < argc; ++index) {
        const std::string key = argv[index];
        if (key.rfind("--", 0) != 0 || index + 1 >= argc)
            throw std::invalid_argument("arguments must be --name value pairs");
        result[key.substr(2)] = argv[++index];
    }
    return result;
}

std::string get(const std::map<std::string, std::string> &args,
                const std::string &name) {
    const auto found = args.find(name);
    if (found == args.end()) throw std::invalid_argument("missing --" + name);
    return found->second;
}

PureImaginaryTimeProtocol protocolFromArguments(
        const std::map<std::string, std::string> &args, double dt) {
    const std::string kind = get(args, "protocol");
    const double initial = std::stod(get(args, "Vi"));
    const double tau = std::stod(get(args, "tau"));
    if (kind == "quench")
        return PureImaginaryTimeProtocol::suddenQuench(
            initial, std::stod(get(args, "Vf")), tau, dt);
    if (kind == "ramp")
        return PureImaginaryTimeProtocol::linearRamp(
            initial, std::stod(get(args, "R")), tau, dt);
    if (kind == "constant")
        return PureImaginaryTimeProtocol::constant(initial, tau, dt);
    throw std::invalid_argument("unknown protocol");
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const auto args = arguments(argc, argv);
        ModelParameters p;
        p.L = std::stoi(get(args, "L"));
        p.boundary = get(args, "boundary") == "pbc" ? 0 : 1;
        p.hs = std::stoi(get(args, "hs"));
        p.dt = std::stod(get(args, "dt"));
        const int explicitParity = std::stoi(get(args, "parity"));
        const int burn = std::stoi(get(args, "burn"));
        const int measurements = std::stoi(get(args, "measurements"));
        const int block = std::stoi(get(args, "block"));
        const std::uint64_t seed = std::stoull(get(args, "seed"));
        const std::string outputPath = get(args, "output");
        require(p.L >= 2 && p.L <= 6 && p.dt > 0.0 && burn >= 0 &&
                    measurements > 0 && block > 0 &&
                    (explicitParity == 1 || explicitParity == -1),
                "invalid driven QMC arguments");
        const PureImaginaryTimeProtocol protocol =
            protocolFromArguments(args, p.dt);
        const GaussianTrialState initial = makeInitialState(
            p, 0.0, 0.0, explicitParity);

        std::mt19937_64 rng(seed);
        PureDrivenContourResult contour = makeDrivenContour(p, protocol, rng);
        std::vector<int> indices = proposalIndices(contour.configuration);
        require(!indices.empty(), "driven contour has no HS fields");
        int ketCount = 0, braCount = 0;
        for (int index : indices) {
            if (contour.configuration.locations[index].branch == PureBranch::Ket)
                ++ketCount;
            else ++braCount;
        }
        require(ketCount == braCount && ketCount > 0,
                "bra/ket HS field counts are not independent and balanced");

        PureFastOptions options;
        options.weight_mode = PureProjectorWeightMode::RealZ2;
        PureProjectorFastWalker walker(
            initial, std::move(contour.configuration), block,
            PureFastRunMode::FastStrict, options,
            PureFastInitializationPolicy::SequentialAudit);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        long long accepted = 0, attempted = 0, sameProposalMismatch = 0;
        std::uint64_t trajectoryHash = 1469598103934665603ULL;
        auto mix = [&](std::uint64_t value) {
            trajectoryHash ^= value;
            trajectoryHash *= 1099511628211ULL;
        };
        int direction = 1;
        auto sweep = [&]() {
            indices = proposalIndices(walker.configuration());
            auto step = [&](int index) {
                const double u = uniform(rng);
                const PureFastProposalResult result = walker.propose(
                    flipProposal(p, protocol, walker.configuration(), index, u));
                ++attempted;
                accepted += result.accepted;
                if (result.ratio.used_reference)
                    sameProposalMismatch += result.z2_reference_mismatch;
                require(!result.terminated && result.ratio.ok(),
                        "driven proposal reference failed closed");
                mix(index);
                mix(result.accepted);
                mix(std::uint64_t(result.ratio.status));
                mix(walker.configurationHash());
                mix(std::uint64_t(std::int64_t(walker.z2Sign())));
            };
            if (direction > 0)
                for (int index : indices) step(index);
            else
                for (auto iterator = indices.rbegin(); iterator != indices.rend();
                     ++iterator) step(*iterator);
            direction = -direction;
        };
        for (int index = 0; index < burn; ++index) sweep();

        std::ofstream csv(outputPath);
        require(bool(csv), "cannot open driven QMC measurement CSV");
        csv << "seed,measurement,physical_sign,energy,S_pi,S_pi_dq,R_CDW,"
               "fermion_parity,G_0_1_real,G_0_1_imag,G_0_L_real,G_0_L_imag,"
               "configuration_hash,acceptance\n";
        double signDenominator = 0.0;
        int positive = 0, negative = 0;
        double endpointResidual = 0.0, phaseImaginary = 0.0;
        for (int measurement = 0; measurement < measurements; ++measurement) {
            sweep();
            const PureProjectorGreenResult green = walker.measurementGreen();
            require(green.ok(), "driven measurement Green failed");
            const int sign = walker.z2Sign();
            const Observables observable = observe(
                p, protocol.finalValue(), green.green);
            signDenominator += sign;
            positive += sign > 0;
            negative += sign < 0;
            endpointResidual = std::max(endpointResidual, green.green_residual);
            phaseImaginary = std::max(
                phaseImaginary,
                std::abs(walker.currentWeight().complex_phase.imag()));
            csv << seed << ',' << measurement << ',' << sign << ','
                << std::setprecision(17) << observable.energy << ','
                << observable.S_pi << ',' << observable.S_pi_dq << ','
                << observable.R_CDW << ',' << observable.parity << ','
                << green.green(0, 1).real() << ',' << green.green(0, 1).imag()
                << ',' << green.green(0, p.L).real() << ','
                << green.green(0, p.L).imag() << ','
                << walker.configurationHash() << ','
                << double(accepted) / attempted << '\n';
        }
        csv.flush();
        require(bool(csv), "driven QMC measurement CSV write failed");
        const PureFastDiagnostics &diagnostics = walker.diagnostics();
        std::cout << std::setprecision(17)
                  << "{\"status\":\"complete\",\"L\":" << p.L
                  << ",\"boundary\":\"" << (p.boundary ? "obc" : "pbc")
                  << "\",\"protocol\":\"" << protocol.name() << "\""
                  << ",\"V_i\":" << protocol.initialValue()
                  << ",\"V_f\":" << protocol.finalValue()
                  << ",\"R\":" << protocol.rate()
                  << ",\"tau_f\":" << protocol.tauFinal()
                  << ",\"Delta_tau\":" << p.dt
                  << ",\"Phi_0_parity\":" << explicitParity
                  << ",\"seed\":" << seed
                  << ",\"burn_sweeps\":" << burn
                  << ",\"measurements\":" << measurements
                  << ",\"bra_ket_hs_independent\":true"
                  << ",\"ket_hs_fields\":" << ketCount
                  << ",\"bra_hs_fields\":" << braCount
                  << ",\"positive_weights\":" << positive
                  << ",\"negative_weights\":" << negative
                  << ",\"sign_denominator\":" << signDenominator
                  << ",\"average_sign\":" << signDenominator / measurements
                  << ",\"acceptance\":" << double(accepted) / attempted
                  << ",\"trajectory_hash\":" << trajectoryHash
                  << ",\"final_configuration_hash\":"
                  << walker.configurationHash()
                  << ",\"rng_hash\":" << hashRng(rng)
                  << ",\"same_proposal_z2_mismatch\":"
                  << sameProposalMismatch
                  << ",\"mp_fallback_count\":"
                  << diagnostics.mp_fallback_count
                  << ",\"mp_fallback_failure_count\":"
                  << diagnostics.mp_fallback_failure_count
                  << ",\"endpoint_green_residual_max\":"
                  << endpointResidual
                  << ",\"complex_phase_imaginary_max\":"
                  << phaseImaginary << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "driven_qmc: " << error.what() << '\n';
        return 1;
    }
}
