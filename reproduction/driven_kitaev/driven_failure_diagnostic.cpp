#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pfqmc.h"
#include "../projector_kitaev/projector_contour.h"
#include "multiprecision_driven_rebuild.h"

namespace {

constexpr double kFailureTolerance = 1e-7;
constexpr double kJackknifeDenominatorTolerance = 1e-12;
constexpr int kSignRecomputeStride = 20;
constexpr double kSignCorrectionTolerance = 1e-2;
constexpr int kDefaultStabilizationInterval = 10;

struct DiagnosticStop {};

struct Args {
    int L;
    int boundary;
    int burn;
    int measurements;
    int seed;
    int ts_stride = 0;
    double V0;
    double Vf;
    double rate;
    double theta_init;
    double beta_trial;
    double dt;
    double delta;
    double mu;
    bool adaptive_guard = false;
    double guard_threshold = 0.8;
    std::string ts_path;
};

struct BinRecord {
    int sample_count = 0;
    double sign_sum = 0.0;
    double signed_S_pi_numerator = 0.0;
    double signed_S_pi_dq_numerator = 0.0;
};

Args parse(int count, char **values) {
    if (count != 14 && count != 16 && count != 18) {
        throw std::runtime_error(
            "usage: driven_driver L V0 Vf drive_rate theta_init beta_trial dt "
            "delta mu boundary burn measurements seed "
            "[adaptive_guard guard_threshold] [timeseries.csv stride]");
    }
    Args args{std::stoi(values[1]), std::stoi(values[10]), std::stoi(values[11]),
              std::stoi(values[12]), std::stoi(values[13]), 0,
              std::stod(values[2]), std::stod(values[3]), std::stod(values[4]),
              std::stod(values[5]), std::stod(values[6]), std::stod(values[7]),
              std::stod(values[8]), std::stod(values[9]), false, 0.8, ""};
    int next = 14;
    if (count >= 16 && (std::string(values[14]) == "0" || std::string(values[14]) == "1")) {
        args.adaptive_guard = std::stoi(values[14]) != 0;
        args.guard_threshold = std::stod(values[15]);
        next = 16;
    }
    if (count - next == 2) {
        args.ts_path = values[next];
        args.ts_stride = std::stoi(values[next + 1]);
    }
    if (args.L < 2 || args.V0 < 0 || args.Vf < 0 || args.theta_init <= 0 ||
        args.beta_trial <= 0 || args.dt <= 0 || args.guard_threshold <= 0 ||
        (args.boundary != 0 && args.boundary != 1) || args.burn < 0 ||
        args.measurements <= 0 || args.ts_stride < 0 ||
        (args.ts_stride > 0 && args.ts_path.empty())) {
        throw std::runtime_error("invalid parameter");
    }
    const double thetaSlices = args.theta_init / args.dt;
    if (std::abs(std::round(thetaSlices) - thetaSlices) > 1e-10) {
        throw std::runtime_error("theta_init/dt must be an integer");
    }
    const double deltaV = args.Vf - args.V0;
    if (std::abs(deltaV) > 1e-12 && (args.rate == 0 || deltaV / args.rate < 0)) {
        throw std::runtime_error("drive_rate must point from V0 to Vf");
    }
    if (std::abs(deltaV) > 1e-12) {
        const double rampSlices = deltaV / (args.rate * args.dt);
        if (std::abs(std::round(rampSlices) - rampSlices) > 1e-10) {
            throw std::runtime_error("(Vf-V0)/(drive_rate*dt) must be an integer");
        }
    }
    return args;
}

double standardError(const std::vector<double> &values) {
    if (values.size() < 2) return 0.0;
    double mean = 0.0;
    for (double value : values) mean += value;
    mean /= values.size();
    double varianceSum = 0.0;
    for (double value : values) varianceSum += (value - mean) * (value - mean);
    return std::sqrt(varianceSum / (values.size() * (values.size() - 1.0)));
}

double jackknifeError(const std::vector<double> &leaveOneOut) {
    if (leaveOneOut.size() < 2) return 0.0;
    double mean = 0.0;
    for (double value : leaveOneOut) mean += value;
    mean /= leaveOneOut.size();
    double sum = 0.0;
    for (double value : leaveOneOut) sum += (value - mean) * (value - mean);
    return std::sqrt((leaveOneOut.size() - 1.0) * sum / leaveOneOut.size());
}

bool finiteComplex(const DataType &value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

std::vector<int> fields(const Spinless_tV &walker) {
    std::vector<int> result;
    for (Operator *op : walker.op_array) {
        if (auto *aux = op->getAuxField()) {
            for (int index = 0; index < aux->size(); ++index) result.push_back((*aux)(index));
        }
    }
    return result;
}

long long changes(const std::vector<int> &before, const std::vector<int> &after) {
    long long result = 0;
    for (size_t index = 0; index < before.size(); ++index) result += before[index] != after[index];
    return result;
}

int configuredBinCount() {
    const char *configured = std::getenv("PFQMC_N_BINS");
    const int requested = configured ? std::stoi(configured) : 15;
    if (requested <= 0) throw std::runtime_error("PFQMC_N_BINS must be positive");
    return requested;
}

std::string binRecordsPath(int seed) {
    if (const char *configured = std::getenv("PFQMC_BIN_RECORDS_PATH")) return configured;
    return "bin_records_seed_" + std::to_string(seed) + ".csv";
}

std::string codeVersion() {
    if (const char *configured = std::getenv("PFQMC_CODE_VERSION")) return configured;
    return "source-tree-no-git";
}

void writeJsonString(std::ostream &out, const std::string &value) {
    out << '"';
    for (char ch : value) {
        if (ch == '"' || ch == '\\') out << '\\';
        out << ch;
    }
    out << '"';
}

class DrivenWalker : public Spinless_tV {
public:
    int center_boundary = -1;
    int n_trial = 0;
    int n_init = 0;
    int n_ramp = 0;
    std::vector<double> ket_schedule;
    std::vector<double> bra_schedule;
    std::vector<const iVecType *> ket_fields;
    std::vector<const iVecType *> bra_fields;
    bool independent = true;

    DrivenWalker(const SpinlessTvChainUtils *config, rdGenerator *random, const Args &args) {
        nDim = config->nDim;
        n_init = int(std::llround(args.theta_init / args.dt));
        n_ramp = std::abs(args.Vf - args.V0) < 1e-12
                     ? 0
                     : int(std::llround((args.Vf - args.V0) / (args.rate * args.dt)));
        if (n_ramp == 0) {
            int physicalSlices = 0;
            build_projector_static_contour(
                *this, config, random, args.theta_init, args.beta_trial,
                center_boundary, n_trial, physicalSlices);
            return;
        }
        int bonds[2];
        if (config->boundaryType == 0) {
            bonds[0] = (config->Lx + 1) / 2;
            bonds[1] = config->Lx / 2;
        } else {
            bonds[0] = config->Lx / 2;
            bonds[1] = (config->Lx - 1) / 2;
        }
        MatType kinetic(nDim, nDim);
        config->KineticGenerator(kinetic);
        for (double remaining = args.beta_trial; remaining > 1e-12;) {
            const double step = std::min(args.dt, remaining);
            MatType exponent = kinetic;
            MatType matrix = expm(exponent, -step);
            exponent = step * kinetic;
            op_array.push_back(new DenseOperator(matrix, signOfHamiltonian(exponent)));
            ++n_trial;
            remaining -= step;
        }
        MatType exponent = kinetic;
        MatType halfKinetic = expm(exponent, -args.dt / 2);
        exponent = (args.dt / 2) * kinetic;
        const DataType kineticSign = signOfHamiltonian(exponent);
        const auto addSlice = [&](double interaction, bool dagger,
                                  std::vector<const iVecType *> *trackedFields) {
            op_array.push_back(new DenseOperator(halfKinetic, kineticSign));
            for (int slot = 0; slot < 2; ++slot) {
                const int bond = dagger ? 1 - slot : slot;
                auto *aux = new iVecType(bonds[bond]);
                for (int index = 0; index < aux->size(); ++index) (*aux)(index) = random->rdZ2();
                if (trackedFields) trackedFields->push_back(aux);
                op_array.push_back(new SpinlessVOperator(config, aux, bond, random, interaction));
            }
            op_array.push_back(new DenseOperator(halfKinetic, kineticSign));
        };
        for (int slice = 0; slice < n_init; ++slice) addSlice(args.V0, false, nullptr);
        for (int slice = 0; slice < n_ramp; ++slice) {
            const double interaction = args.V0 + args.rate * (slice + 0.5) * args.dt;
            ket_schedule.push_back(interaction);
            addSlice(interaction, false, &ket_fields);
        }
        center_boundary = int(op_array.size());
        for (int slice = n_ramp - 1; slice >= 0; --slice) {
            const double interaction = args.V0 + args.rate * (slice + 0.5) * args.dt;
            bra_schedule.push_back(interaction);
            addSlice(interaction, true, &bra_fields);
        }
        for (int slice = 0; slice < n_init; ++slice) addSlice(args.V0, true, nullptr);
        for (const auto *ket : ket_fields) {
            for (const auto *bra : bra_fields) {
                if (ket == bra) independent = false;
            }
        }
    }
};

void writeArray(std::ostream &out, const std::vector<double> &values) {
    out << '[';
    for (size_t index = 0; index < values.size(); ++index) {
        if (index) out << ',';
        out << values[index];
    }
    out << ']';
}

double skewSymmetryError(const MatType &g) {
    MatType G = g - MatType::Identity(g.rows(), g.cols());
    return (G + G.transpose()).cwiseAbs().maxCoeff();
}

}  // namespace

int main(int argc, char **argv) try {
    const Args args = parse(argc, argv);
    mkl_set_num_threads(1);
    const auto started = std::chrono::steady_clock::now();
    const int nRamp = std::abs(args.Vf - args.V0) < 1e-12
                        ? 0
                        : int(std::llround((args.Vf - args.V0) / (args.rate * args.dt)));
    SpinlessTvChainUtils config(
        args.L, args.dt, args.V0,
        2 * (int(std::llround(args.theta_init / args.dt)) + nRamp),
        args.boundary, args.delta, args.mu, 0);
    rdGenerator random(args.seed);
    DrivenWalker walker(&config, &random, args);
    PfQMC qmc(&walker, kDefaultStabilizationInterval);
    qmc.configureAdaptiveGuard(args.adaptive_guard, args.guard_threshold, 100.0);

    const char *mpEnv = std::getenv("PFQMC_MULTIPRECISION_FALLBACK");
    const bool multiprecisionEnabled = mpEnv &&
        (std::string(mpEnv) == "ON" || std::string(mpEnv) == "1");
    const char *mpRecordEnv = std::getenv("PFQMC_MULTIPRECISION_RECORD_PROXY");
    const bool multiprecisionRecord = mpRecordEnv &&
        (std::string(mpRecordEnv) == "ON" || std::string(mpRecordEnv) == "1");
    if (multiprecisionEnabled || multiprecisionRecord) {
        if (nRamp == 0 || args.boundary != 1 || config.hsScheme != 0 ||
            std::abs(args.beta_trial / args.dt - std::round(args.beta_trial / args.dt)) > 1e-10) {
            throw std::runtime_error(
                "multiprecision fallback supports only driven OBC complete-slice contours");
        }
        double threshold = 1e4;
        if (const char *value = std::getenv("PFQMC_MULTIPRECISION_CORE_THRESHOLD")) {
            threshold = std::stod(value);
        }
        qmc.configureMultiprecisionFallback(
            multiprecisionEnabled, threshold, [&config, &walker](int boundary, MatType &out) {
                return driven_multiprecision::rebuild(
                    config, walker.op_array, boundary, walker.n_trial, out);
            });
    }

    for (int sweep = 0; sweep < args.burn; ++sweep) {
        qmc.rightSweep();
        qmc.leftSweep();
    }

    std::ofstream timeSeries;
    if (args.ts_stride > 0) {
        timeSeries.open(args.ts_path);
        if (!timeSeries) throw std::runtime_error("cannot open time-series output");
        timeSeries << "measurement,sign,sign_S_pi_numerator,sign_S_pi_dq_numerator,"
                   << "S_pi_estimator,S_pi_dq_estimator,R_cdw,acceptance\n"
                   << std::setprecision(17);
    }

    const int nBinsRequested = configuredBinCount();
    const int nBinsUsed = std::min(nBinsRequested, args.measurements);
    std::vector<BinRecord> bins(nBinsUsed);
    const std::string binsPath = binRecordsPath(args.seed);
    double signSum = 0.0;
    double signedSpiSum = 0.0;
    double signedSpidqSum = 0.0;
    double maxSignImag = 0.0;
    double maxObservableImag = 0.0;
    double maxGreenSkewSymmetryError = 0.0;
    long long accepted = 0;
    long long attempted = 0;
    double minUpdateDenominator = INFINITY;
    int negativeSigns = 0;
    int signRecomputes = 0;
    int signCorrections = 0;
    bool finiteSamples = true;
    int completedMeasurements = 0;
    int currentDiagnosticSample = -1;
    bool diagnosticStopped = false;
    const char *diagnosticPath = std::getenv("PFQMC_FAILURE_DIAGNOSTIC_PATH");
    std::ofstream diagnostic;
    if (diagnosticPath) {
        diagnostic.open(diagnosticPath);
        if (!diagnostic) throw std::runtime_error("cannot open failure diagnostic output");
        diagnostic << "measurement,boundary,trigger,threshold,tracked_sign_real,"
                      "tracked_sign_imag,tracked_eigen_sign,direct_sign_real,"
                      "direct_sign_imag,direct_eigen_sign,eigen_sign_mismatch,"
                      "fast_full_green_relative_error,fast_green_skew,full_green_skew,"
                      "fast_S_pi_real,fast_S_pi_imag,full_S_pi_real,full_S_pi_imag,"
                      "fast_S_pi_dq_real,fast_S_pi_dq_imag,full_S_pi_dq_real,"
                      "full_S_pi_dq_imag\n" << std::setprecision(17);
        qmc.capture_diagnostic_callback = [&](int boundary, const MatType &fastGreen,
                                              DataType trackedSign) {
            const DataType fastSpi = config.StructureFactorCDW(fastGreen);
            const DataType fastSpidq = config.StructureFactorCDWOffset(fastGreen);
            const double fastSkew = skewSymmetryError(fastGreen);
            const bool signTrigger = std::abs(trackedSign.imag()) > kFailureTolerance;
            const bool observableTrigger =
                std::max(std::abs(fastSpi.imag()), std::abs(fastSpidq.imag())) >
                kFailureTolerance;
            const bool greenTrigger = fastSkew > kFailureTolerance;
            if (!signTrigger && !observableTrigger && !greenTrigger) return;

            MatType fullGreen;
            qmc.rebuildGreenFromFullContourAtBoundary(boundary, fullGreen);
            const DataType directSign = qmc.getSignRaw();
            const DataType fullSpi = config.StructureFactorCDW(fullGreen);
            const DataType fullSpidq = config.StructureFactorCDWOffset(fullGreen);
            const double relativeGreenError = (fastGreen - fullGreen).norm() /
                std::max(fullGreen.norm(), std::numeric_limits<double>::min());
            const int trackedEigenSign = trackedSign.real() >= 0 ? 1 : -1;
            const int directEigenSign = directSign.real() >= 0 ? 1 : -1;
            std::string trigger;
            if (signTrigger) trigger += "SIGN_IMAG";
            if (observableTrigger) trigger += (trigger.empty() ? "" : "+") +
                                                std::string("OBSERVABLE_IMAG");
            if (greenTrigger) trigger += (trigger.empty() ? "" : "+") +
                                           std::string("GREEN_SKEW");
            diagnostic << currentDiagnosticSample << ',' << boundary << ',' << trigger
                       << ',' << kFailureTolerance << ',' << trackedSign.real() << ','
                       << trackedSign.imag() << ',' << trackedEigenSign << ','
                       << directSign.real() << ',' << directSign.imag() << ','
                       << directEigenSign << ',' << (trackedEigenSign != directEigenSign)
                       << ',' << relativeGreenError << ',' << fastSkew << ','
                       << skewSymmetryError(fullGreen) << ',' << fastSpi.real() << ','
                       << fastSpi.imag() << ',' << fullSpi.real() << ',' << fullSpi.imag()
                       << ',' << fastSpidq.real() << ',' << fastSpidq.imag() << ','
                       << fullSpidq.real() << ',' << fullSpidq.imag() << '\n';
            diagnostic.flush();
            throw DiagnosticStop{};
        };
    }

    for (int sample = 0; sample < args.measurements; ++sample) {
        currentDiagnosticSample = sample;
        if (sample % kSignRecomputeStride == 0) {
            const PfaffianResult raw = qmc.getSignRawWithStatus();
            ++signRecomputes;
            if (raw.ok() &&
                (qmc.sign.real() >= 0) != (raw.value.real() >= 0)) {
                ++signCorrections;
            }
        }
        const auto before = fields(walker);
        MatType g;
        DataType capturedSign;
        try {
            qmc.rightSweep(walker.center_boundary, &g, &capturedSign);
        } catch (const DiagnosticStop &) {
            diagnosticStopped = true;
            break;
        }
        const auto afterRight = fields(walker);
        long long sampleAccepted = changes(before, afterRight);
        long long sampleAttempted = before.size();
        accepted += sampleAccepted;
        attempted += sampleAttempted;

        const DataType rawSpi = config.StructureFactorCDW(g);
        const DataType rawSpidq = config.StructureFactorCDWOffset(g);
        const bool finite = finiteComplex(capturedSign) && finiteComplex(rawSpi) &&
                            finiteComplex(rawSpidq) &&
                            std::isfinite(skewSymmetryError(g));
        if (!finite) {
            finiteSamples = false;
            break;
        }
        for (Operator *op : walker.op_array) {
            if (auto *interaction = dynamic_cast<SpinlessVOperator *>(op)) {
                const double denominator = std::min(
                    std::abs(interaction->last_denom1), std::abs(interaction->last_denom2));
                minUpdateDenominator = std::min(minUpdateDenominator, denominator);
            }
        }

        const double sampledSign = capturedSign.real() >= 0 ? 1.0 : -1.0;
        const double signedSpi = -sampledSign * rawSpi.real();
        const double signedSpidq = -sampledSign * rawSpidq.real();
        signSum += sampledSign;
        signedSpiSum += signedSpi;
        signedSpidqSum += signedSpidq;
        if (sampledSign < 0) ++negativeSigns;
        const int bin = std::min(
            nBinsUsed - 1,
            int((static_cast<long long>(sample) * nBinsUsed) / args.measurements));
        ++bins[bin].sample_count;
        bins[bin].sign_sum += sampledSign;
        bins[bin].signed_S_pi_numerator += signedSpi;
        bins[bin].signed_S_pi_dq_numerator += signedSpidq;
        maxSignImag = std::max(maxSignImag, std::abs(capturedSign.imag()));
        maxObservableImag = std::max(
            maxObservableImag, std::max(std::abs(rawSpi.imag()), std::abs(rawSpidq.imag())));
        maxGreenSkewSymmetryError = std::max(maxGreenSkewSymmetryError, skewSymmetryError(g));
        ++completedMeasurements;

        const auto beforeLeft = afterRight;
        qmc.leftSweep();
        const auto afterLeft = fields(walker);
        const long long leftAccepted = changes(beforeLeft, afterLeft);
        const long long leftAttempted = beforeLeft.size();
        accepted += leftAccepted;
        attempted += leftAttempted;
        sampleAccepted += leftAccepted;
        sampleAttempted += leftAttempted;
        if (args.ts_stride > 0 && sample % args.ts_stride == 0) {
            const double spi = -rawSpi.real();
            const double spidq = -rawSpidq.real();
            timeSeries << sample << ',' << sampledSign << ',' << signedSpi << ','
                       << signedSpidq << ',' << spi << ',' << spidq << ','
                       << (1.0 - spidq / spi) << ','
                       << double(sampleAccepted) / sampleAttempted << '\n';
        }
    }

    if (diagnosticStopped) finiteSamples = false;

    std::ofstream binRecords(binsPath);
    if (!binRecords) throw std::runtime_error("cannot open bin-record output");
    binRecords << "bin_index,sample_count,sign_sum,signed_S_pi_numerator,"
               << "signed_S_pi_dq_numerator\n" << std::setprecision(17);
    for (int index = 0; index < nBinsUsed; ++index) {
        const BinRecord &bin = bins[index];
        binRecords << index << ',' << bin.sample_count << ',' << bin.sign_sum << ','
                   << bin.signed_S_pi_numerator << ','
                   << bin.signed_S_pi_dq_numerator << '\n';
    }

    bool jackknifeValid = finiteSamples && completedMeasurements == args.measurements &&
                          std::isfinite(signSum) &&
                          std::abs(signSum) > kJackknifeDenominatorTolerance &&
                          std::isfinite(signedSpiSum) &&
                          std::isfinite(signedSpidqSum) &&
                          std::abs(signedSpiSum) > kJackknifeDenominatorTolerance;
    std::string failureReason;
    if (!finiteSamples) failureReason = "nonfinite_sign_or_observable";
    if (completedMeasurements != args.measurements && failureReason.empty()) {
        failureReason = "incomplete_measurement_loop";
    }
    if (!jackknifeValid && failureReason.empty()) {
        failureReason = "pooled_jackknife_denominator_below_tolerance";
    }

    std::vector<double> jackknifeSpi;
    std::vector<double> jackknifeSpidq;
    std::vector<double> jackknifeR;
    if (jackknifeValid) {
        for (const BinRecord &bin : bins) {
            const double leaveSign = signSum - bin.sign_sum;
            const double leaveSpiNumerator = signedSpiSum - bin.signed_S_pi_numerator;
            const double leaveSpidqNumerator = signedSpidqSum - bin.signed_S_pi_dq_numerator;
            if (!std::isfinite(leaveSign) ||
                std::abs(leaveSign) <= kJackknifeDenominatorTolerance ||
                !std::isfinite(leaveSpiNumerator) ||
                std::abs(leaveSpiNumerator) <= kJackknifeDenominatorTolerance ||
                !std::isfinite(leaveSpidqNumerator)) {
                jackknifeValid = false;
                failureReason = "leave_one_out_denominator_below_tolerance";
                break;
            }
            const double leaveSpi = leaveSpiNumerator / leaveSign;
            const double leaveSpidq = leaveSpidqNumerator / leaveSign;
            jackknifeSpi.push_back(leaveSpi);
            jackknifeSpidq.push_back(leaveSpidq);
            jackknifeR.push_back(1.0 - leaveSpidqNumerator / leaveSpiNumerator);
        }
    }

    const double contact = 1.0 / (4.0 * args.L);
    const double spi = jackknifeValid ? signedSpiSum / signSum : 0.0;
    const double spidq = jackknifeValid ? signedSpidqSum / signSum : 0.0;
    const double rCdw = jackknifeValid ? 1.0 - signedSpidqSum / signedSpiSum : 0.0;
    const double spiError = jackknifeValid ? jackknifeError(jackknifeSpi) : 0.0;
    const double spidqError = jackknifeValid ? jackknifeError(jackknifeSpidq) : 0.0;
    const double rError = jackknifeValid ? jackknifeError(jackknifeR) : 0.0;

    const bool numericalTolerancePass = finiteSamples &&
        maxSignImag <= kFailureTolerance &&
        maxObservableImag <= kFailureTolerance &&
        maxGreenSkewSymmetryError <= kFailureTolerance;
    if (!numericalTolerancePass && failureReason.empty()) {
        failureReason = "numerical_diagnostic_exceeds_tolerance";
    }
    const int exitCode = jackknifeValid && numericalTolerancePass ? 0 : 3;
    const std::string status = exitCode == 0 ? "complete" : "failed";
    const double runtime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    const double reportedMinDenominator =
        args.adaptive_guard ? qmc.min_update_denominator : minUpdateDenominator;
    const double guardFrequency = qmc.proposal_attempt_count
        ? double(qmc.pre_decision_rebuild_count) / qmc.proposal_attempt_count
        : 0.0;
    const auto conditionQuantile = [&qmc](double quantile) {
        if (qmc.multiprecision_condition_samples.empty()) return 0.0;
        auto values = qmc.multiprecision_condition_samples;
        std::sort(values.begin(), values.end());
        const size_t index = std::min(
            values.size() - 1, size_t(std::ceil(quantile * values.size()) - 1));
        return values[index];
    };

    std::cout << std::setprecision(17)
              << "{\"output_schema_version\":\"2.0\","
              << "\"observable_convention_version\":\"driven_cdw_ratio_of_means_v1\","
              << "\"status\":";
    writeJsonString(std::cout, status);
    std::cout << ",\"failure_reason\":";
    writeJsonString(std::cout, failureReason);
    std::cout << ",\"code_version\":";
    writeJsonString(std::cout, codeVersion());
    std::cout << ",\"mode\":\"driven_projector\",\"L\":" << args.L
              << ",\"V0\":" << args.V0 << ",\"Vf\":" << args.Vf
              << ",\"drive_rate\":" << args.rate << ",\"ramp_time\":" << nRamp * args.dt
              << ",\"n_ramp\":" << nRamp << ",\"theta_init\":" << args.theta_init
              << ",\"beta_trial\":" << args.beta_trial << ",\"dt\":" << args.dt
              << ",\"delta\":" << args.delta << ",\"mu\":" << args.mu
              << ",\"boundary\":" << args.boundary << ",\"seed\":" << args.seed
              << ",\"burn\":" << args.burn << ",\"measurements\":" << args.measurements
              << ",\"measurements_completed\":" << completedMeasurements
              << ",\"trial_slices\":" << walker.n_trial
              << ",\"initial_slices_per_side\":" << walker.n_init
              << ",\"center_operator_boundary\":" << walker.center_boundary
              << ",\"operator_count\":" << walker.op_array.size()
              << ",\"stabilization_interval\":" << qmc.stb
              << ",\"sign_recompute_stride\":" << kSignRecomputeStride
              << ",\"sign_correction_tolerance\":" << kSignCorrectionTolerance
              << ",\"n_bins_requested\":" << nBinsRequested
              << ",\"n_bins_used\":" << nBinsUsed
              << ",\"jackknife_denominator_tolerance\":" << kJackknifeDenominatorTolerance
              << ",\"bin_records_path\":";
    writeJsonString(std::cout, binsPath);
    std::cout << ",\"ket_schedule\":";
    writeArray(std::cout, walker.ket_schedule);
    std::cout << ",\"bra_schedule\":";
    writeArray(std::cout, walker.bra_schedule);
    std::cout << ",\"ket_slice_bond_order\":\"even_then_odd\""
              << ",\"bra_slice_bond_order\":\"odd_then_even\""
              << ",\"bra_ket_hs_independent\":" << (walker.independent ? "true" : "false")
              << ",\"onsite_contact\":" << contact
              << ",\"onsite_contact_is_diagnostic_only\":true"
              << ",\"S_pi_offsite\":" << spi - contact
              << ",\"S_pi_offsite_mean\":" << spi - contact
              << ",\"S_pi_offsite_err\":" << spiError
              << ",\"S_pi_dq_offsite\":" << spidq - contact
              << ",\"S_pi_dq_offsite_mean\":" << spidq - contact
              << ",\"S_pi_dq_offsite_err\":" << spidqError
              << ",\"S_pi\":" << spi << ",\"S_pi_mean\":" << spi
              << ",\"S_pi_err\":" << spiError
              << ",\"S_pi_dq\":" << spidq << ",\"S_pi_dq_mean\":" << spidq
              << ",\"S_pi_dq_err\":" << spidqError
              << ",\"R_cdw\":" << rCdw << ",\"R_cdw_mean\":" << rCdw
              << ",\"R_cdw_err\":" << rError
              << ",\"average_sign\":" << signSum / std::max(1, completedMeasurements)
              << ",\"average_sign_mean\":" << signSum / std::max(1, completedMeasurements)
              << ",\"average_sign_err\":" << 0.0
              << ",\"acceptance\":" << (attempted ? double(accepted) / attempted : 0.0)
              << ",\"runtime_seconds\":" << runtime
              << ",\"negative_signs\":" << negativeSigns
              << ",\"sign_recomputes\":" << signRecomputes
              << ",\"sign_corrections\":" << signCorrections
              << ",\"max_sign_imag\":" << maxSignImag
              << ",\"max_observable_imag\":" << maxObservableImag
              << ",\"max_green_skew_symmetry_error\":" << maxGreenSkewSymmetryError
              << ",\"sign_imag_tolerance\":" << kFailureTolerance
              << ",\"observable_imag_tolerance\":" << kFailureTolerance
              << ",\"green_skew_symmetry_tolerance\":" << kFailureTolerance
              << ",\"adaptive_guard\":" << (args.adaptive_guard ? "true" : "false")
              << ",\"adaptive_rebuild_count\":" << qmc.adaptive_rebuild_count
              << ",\"pre_decision_rebuild_count\":" << qmc.pre_decision_rebuild_count
              << ",\"post_accept_rebuild_count\":" << qmc.post_accept_rebuild_count
              << ",\"multiprecision_fallback\":" << (multiprecisionEnabled ? "true" : "false")
              << ",\"multiprecision_fallback_count\":" << qmc.multiprecision_fallback_count
              << ",\"multiprecision_proxy_trigger_count\":" << qmc.multiprecision_proxy_trigger_count
              << ",\"multiprecision_condition_samples\":"
              << qmc.multiprecision_condition_samples.size()
              << ",\"multiprecision_condition_p99\":" << conditionQuantile(.99)
              << ",\"multiprecision_condition_p999\":" << conditionQuantile(.999)
              << ",\"multiprecision_condition_max\":" << conditionQuantile(1.0)
              << ",\"guard_trigger_frequency\":" << guardFrequency
              << ",\"min_update_denominator\":" << reportedMinDenominator
              << ",\"guard_threshold\":" << args.guard_threshold << "}\n";
    return exitCode;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
}
