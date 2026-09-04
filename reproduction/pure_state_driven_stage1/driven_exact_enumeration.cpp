#include "stage1_dense.h"

#include <fstream>
#include <iostream>

using namespace driven_stage1;

int main(int argc, char **argv) {
    try {
        if (argc != 2) throw std::invalid_argument(
            "usage: driven_exact_enumeration output.csv");
        ModelParameters p;
        p.L = 2;
        p.boundary = 0;
        p.hs = 0;
        p.dt = 0.1;
        const PureImaginaryTimeProtocol protocol =
            PureImaginaryTimeProtocol::linearRamp(0.0, 20.0, 0.2, p.dt);
        const double trialMu = 0.3;
        const GaussianTrialState initial = makeInitialState(p, trialMu);
        const int initialParity = initial.fermionParity();
        const std::vector<MatType> gamma = gammaMatrices(p.L);
        const DenseInitialState denseInitial = denseInitialState(
            p, gamma, initial, trialMu, initialParity);

        std::mt19937_64 rng(720002);
        PureDrivenContourResult built = makeDrivenContour(p, protocol, rng);
        const std::vector<int> indices = proposalIndices(built.configuration);
        require(indices.size() <= 20, "enumeration contour is too large");
        const std::uint64_t configurations = 1ULL << indices.size();
        PureProjectorWeightOptions weightOptions;
        weightOptions.mode = PureProjectorWeightMode::RealZ2;
        PureProjectorWeightEvaluator evaluator(weightOptions);

        std::ofstream csv(argv[1]);
        require(bool(csv), "cannot open driven_exact_enumeration.csv");
        csv << "configuration,qmc_weight_real,qmc_weight_imag,dense_weight_real,"
               "dense_weight_imag,weight_abs_error,qmc_z2,dense_z2,z2_mismatch,"
               "qmc_G_0_1_real,qmc_G_0_1_imag,dense_G_0_1_real,"
               "dense_G_0_1_imag,green_max_abs_error,local_ratio_max_abs_error,status\n";
        double maxWeightError = 0.0, maxGreenError = 0.0, maxRatioError = 0.0;
        long long signMismatch = 0;
        DataType qmcSignedSum = 0.0, denseSignedSum = 0.0;
        double qmcAbsoluteSum = 0.0, denseAbsoluteSum = 0.0;
        MatType weightedGreen = MatType::Zero(2 * p.L, 2 * p.L);

        for (std::uint64_t mask = 0; mask < configurations; ++mask) {
            PureFastConfiguration configuration = built.configuration;
            for (std::size_t bit = 0; bit < indices.size(); ++bit) {
                const int index = indices[bit];
                const int sigma = ((mask >> bit) & 1ULL) ? 1 : -1;
                const PureSliceLocation &location = configuration.locations[index];
                MatType factor = localHsFactor(
                    p, protocol.midpointValue(location.slice), location.bond,
                    location.aux, sigma);
                if (location.branch == PureBranch::Bra) factor.adjointInPlace();
                configuration.hs_fields[index] = sigma;
                configuration.slices[index].matrix = std::move(factor);
            }
            const PureProjectorWeightResult qmc =
                evaluator.evaluate(initial, configuration.slices);
            require(qmc.ok(), "enumerated PfQMC weight failed");
            PureProjectorStackManager stack(
                initial, configuration.slices, 4, PureProjectorOptions(),
                int(configuration.slices.size() / 2));
            require(stack.ok() && stack.green().ok(),
                    "enumerated PfQMC center Green failed");
            const MatType qmcGreen = stack.green().green;
            const DenseContourResult dense = denseContour(
                p, protocol, initial, denseInitial.vector, configuration);
            const double weightError = std::abs(qmc.weight - dense.weight);
            const double greenError =
                (qmcGreen - dense.green).cwiseAbs().maxCoeff();
            const int denseSign = dense.weight.real() >= 0.0 ? 1 : -1;
            const int mismatch = qmc.z2_sign != denseSign;
            double configurationRatioError = 0.0;
            for (int index : indices) {
                PureFastConfiguration candidate = configuration;
                const PureSliceLocation &location = candidate.locations[index];
                const int sigma = -candidate.hs_fields[index];
                MatType factor = localHsFactor(
                    p, protocol.midpointValue(location.slice), location.bond,
                    location.aux, sigma);
                if (location.branch == PureBranch::Bra) factor.adjointInPlace();
                candidate.hs_fields[index] = sigma;
                candidate.slices[index].matrix = std::move(factor);
                const PureProjectorWeightResult qmcCandidate =
                    evaluator.evaluate(initial, candidate.slices);
                require(qmcCandidate.ok(), "enumerated candidate PfQMC weight failed");
                const DenseContourResult denseCandidate = denseContour(
                    p, protocol, initial, denseInitial.vector, candidate);
                const DataType qmcRatio = qmcCandidate.weight / qmc.weight;
                const DataType denseRatio = denseCandidate.weight / dense.weight;
                configurationRatioError = std::max(
                    configurationRatioError, std::abs(qmcRatio - denseRatio));
            }
            maxWeightError = std::max(maxWeightError, weightError);
            maxGreenError = std::max(maxGreenError, greenError);
            maxRatioError = std::max(maxRatioError, configurationRatioError);
            signMismatch += mismatch;
            qmcSignedSum += qmc.weight;
            denseSignedSum += dense.weight;
            qmcAbsoluteSum += std::abs(qmc.weight);
            denseAbsoluteSum += std::abs(dense.weight);
            weightedGreen += qmc.weight * qmcGreen;
            const bool pass = weightError < 1e-10 && greenError < 1e-10 &&
                              configurationRatioError < 1e-10 && !mismatch;
            csv << mask << ',' << std::setprecision(17)
                << qmc.weight.real() << ',' << qmc.weight.imag() << ','
                << dense.weight.real() << ',' << dense.weight.imag() << ','
                << weightError << ',' << qmc.z2_sign << ',' << denseSign << ','
                << mismatch << ',' << qmcGreen(0, 1).real() << ','
                << qmcGreen(0, 1).imag() << ',' << dense.green(0, 1).real()
                << ',' << dense.green(0, 1).imag() << ',' << greenError << ','
                << configurationRatioError
                << ',' << (pass ? "PASS" : "FAIL") << '\n';
            require(pass, "enumeration differs from dense Fock reference");
        }
        csv.flush();
        require(bool(csv), "driven_exact_enumeration.csv write failed");
        const DataType qmcAverageSign = qmcSignedSum / qmcAbsoluteSum;
        const DataType denseAverageSign = denseSignedSum / denseAbsoluteSum;
        const double averageSignError =
            std::abs(qmcAverageSign - denseAverageSign);
        require(averageSignError < 1e-11,
                "enumerated average sign differs from dense Fock reference");
        const MatType enumeratedGreen = weightedGreen / qmcSignedSum;
        double sameContourDenominator = 0.0;
        const cVecType exactState = denseSameContourState(
            p, protocol, gamma, denseInitial.vector, &sameContourDenominator);
        const MatType exactGreen = denseStateGreen(exactState, gamma);
        const double ensembleGreenError =
            (enumeratedGreen - exactGreen).cwiseAbs().maxCoeff();
        require(ensembleGreenError < 1e-10,
                "enumerated auxiliary fields disagree with time-ordered ED");
        std::cout << std::setprecision(17)
                  << "{\"status\":\"PASS\",\"L\":2,"
                  << "\"boundary\":\"pbc\",\"protocol\":\"linear_ramp\","
                  << "\"V_i\":0,\"R\":20,\"tau_f\":0.2,\"dt\":0.1,"
                  << "\"Phi_0_parity\":" << initialParity
                  << ",\"configurations\":" << configurations
                  << ",\"bra_ket_hs_independent\":true,"
                  << "\"max_weight_abs_error\":" << maxWeightError
                  << ",\"max_green_abs_error\":" << maxGreenError
                  << ",\"max_local_ratio_abs_error\":" << maxRatioError
                  << ",\"z2_mismatch\":" << signMismatch
                  << ",\"qmc_average_sign\":" << qmcAverageSign.real()
                  << ",\"dense_average_sign\":" << denseAverageSign.real()
                  << ",\"average_sign_abs_error\":" << averageSignError
                  << ",\"ensemble_green_vs_time_ordered_ed_max_abs_error\":"
                  << ensembleGreenError
                  << ",\"same_contour_denominator\":"
                  << sameContourDenominator
                  << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "driven_exact_enumeration: " << error.what() << '\n';
        return 1;
    }
}
