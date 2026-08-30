#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include "kitaevChain.h"
#include "pure_projector_static.h"

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}

double relativeError(DataType a, DataType b) {
    return std::abs(a - b) / std::max(1e-14, std::abs(b));
}

double relativeError(const MatType &a, const MatType &b) {
    return (a - b).norm() / std::max(1.0, b.norm());
}

MatType gammaMatrix(int modes, int which) {
    const int dimension = 1 << modes;
    const int site = which % modes;
    const bool second = which >= modes;
    MatType gamma = MatType::Zero(dimension, dimension);
    for (int state = 0; state < dimension; ++state) {
        const int flipped = state ^ (1 << site);
        const int parity = __builtin_popcount(unsigned(state & ((1 << site) - 1))) & 1;
        const double sign = parity ? -1.0 : 1.0;
        gamma(flipped, state) = second
            ? DataType(0.0, (state & (1 << site)) ? sign : -sign)
            : DataType(sign, 0.0);
    }
    return gamma;
}

std::vector<MatType> gammaMatrices(int modes) {
    std::vector<MatType> result;
    for (int index = 0; index < 2 * modes; ++index)
        result.push_back(gammaMatrix(modes, index));
    return result;
}

MatType fockFactor(const MatType &majoranaExponent,
                   const std::vector<MatType> &gamma) {
    MatType generator = MatType::Zero(gamma[0].rows(), gamma[0].cols());
    for (int i = 0; i < majoranaExponent.rows(); ++i)
        for (int j = i + 1; j < majoranaExponent.cols(); ++j)
            // With the Jordan-Wigner Majoranas used here,
            // exp(-1/2 A_ij gamma_i gamma_j) implements the single-particle
            // left action B=exp(A).  The minus sign is fixed by
            // [gamma_i gamma_j,gamma_i]=-2 gamma_j.
            generator -= 0.5 * majoranaExponent(i, j) * gamma[i] * gamma[j];
    MatType copy = generator;
    return expm(copy, 1.0);
}

struct TestFactor {
    PureProjectorSlice slice;
    MatType fock;
};

TestFactor factorFromExponent(const MatType &exponent,
                              const std::vector<MatType> &gamma,
                              const std::string &label) {
    MatType copy = exponent;
    return {PureProjectorSlice(expm(copy, 1.0), DataType(1.0), label),
            fockFactor(exponent, gamma)};
}

MatType hsExponent(const SpinlessTvChainUtils &config, int bondType,
                   const std::vector<int> &fields) {
    MatType exponent = MatType::Zero(config.nDim, config.nDim);
    const double lambda = std::acosh(std::exp(0.5 * config.V * config.dt));
    for (int aux = 0; aux < int(fields.size()); ++aux) {
        int a, b, c, d;
        config.aux2MajoranaIdx(aux, 0, bondType, a, b);
        config.aux2MajoranaIdx(aux, 1, bondType, c, d);
        const DataType value(0.0, lambda * fields[aux]);
        if (config.hsScheme == 0) {
            exponent(a, b) = value; exponent(b, a) = -value;
            exponent(c, d) = value; exponent(d, c) = -value;
        } else {
            exponent(a, c) = value; exponent(c, a) = -value;
            exponent(b, d) = -value; exponent(d, b) = value;
        }
    }
    return exponent;
}

std::pair<int, int> bondCounts(const SpinlessTvChainUtils &config) {
    return pureProjectorCheckerboardBondCounts(config.Lx, config.boundaryType);
}

struct BuiltContour {
    std::vector<TestFactor> ket;
    std::vector<TestFactor> braAction;
    std::vector<PureProjectorSlice> fullSlices;
    MatType fockKet;
    MatType fockBra;
};

BuiltContour buildContour(const SpinlessTvChainUtils &config,
                          const std::vector<int> &bits,
                          const std::vector<MatType> &gamma) {
    const auto counts = bondCounts(config);
    const int fieldsPerBranch = counts.first + counts.second;
    require(int(bits.size()) == 2 * fieldsPerBranch, "field count mismatch");
    MatType kinetic(config.nDim, config.nDim);
    config.KineticGenerator(kinetic);
    const MatType halfExponent = -0.5 * config.dt * kinetic;
    const TestFactor half = factorFromExponent(halfExponent, gamma, "K/2");

    auto branch = [&](int offset, const std::string &prefix) {
        std::vector<int> first(bits.begin() + offset,
                               bits.begin() + offset + counts.first);
        std::vector<int> second(bits.begin() + offset + counts.first,
                                bits.begin() + offset + fieldsPerBranch);
        std::vector<TestFactor> factors;
        factors.push_back(half);
        factors.push_back(factorFromExponent(hsExponent(config, 0, first), gamma,
                                             prefix + ":V0"));
        factors.push_back(factorFromExponent(hsExponent(config, 1, second), gamma,
                                             prefix + ":V1"));
        factors.push_back(half);
        return factors;
    };

    BuiltContour result;
    result.ket = branch(0, "ket");
    // The bra protocol is generated independently and flattened in strict
    // adjoint/reverse bond-factor order.  Its action on the ket is U_bra.
    std::vector<TestFactor> braProtocol = branch(fieldsPerBranch, "bra");
    result.braAction.assign(braProtocol.rbegin(), braProtocol.rend());
    for (const TestFactor &factor : result.ket) result.fullSlices.push_back(factor.slice);
    for (const TestFactor &factor : result.braAction) result.fullSlices.push_back(factor.slice);
    result.fockKet = MatType::Identity(1 << config.Lx, 1 << config.Lx);
    for (const TestFactor &factor : result.ket) result.fockKet = result.fockKet * factor.fock;
    result.fockBra = MatType::Identity(1 << config.Lx, 1 << config.Lx);
    for (const TestFactor &factor : result.braAction)
        result.fockBra = result.fockBra * factor.fock;
    return result;
}

MatType denseGreen(const cVecType &trial, const MatType &fockProduct,
                   const std::vector<MatType> &gamma, DataType weight) {
    MatType green = MatType::Zero(gamma.size(), gamma.size());
    for (int i = 0; i < int(gamma.size()); ++i)
        for (int j = 0; j < int(gamma.size()); ++j) {
            if (i == j) continue;
            // pureProjectorGreen(phiRight,phiLeft) is the cut contraction
            // <left|gamma_i gamma_j|right>, hence gamma gamma precedes the
            // ket-side product at the end cut.
            green(i, j) = -(trial.adjoint() * gamma[i] * gamma[j] * fockProduct * trial)(0) /
                          weight;
        }
    return green;
}

cVecType denseTrialVector(const GaussianTrialState &trial, const MatType &hTrial,
                          const std::vector<MatType> &gamma) {
    MatType denseHamiltonian = MatType::Zero(gamma[0].rows(), gamma[0].cols());
    for (int i = 0; i < hTrial.rows(); ++i)
        for (int j = i + 1; j < hTrial.cols(); ++j)
            denseHamiltonian += 0.5 * hTrial(i, j) * gamma[i] * gamma[j];
    Eigen::SelfAdjointEigenSolver<MatType> solver(denseHamiltonian);
    require(solver.info() == Eigen::Success, "dense trial eigensolver failed");
    double best = std::numeric_limits<double>::infinity();
    cVecType selected;
    const MatType identity = MatType::Identity(denseHamiltonian.rows(), denseHamiltonian.cols());
    for (int column = 0; column < solver.eigenvectors().cols(); ++column) {
        cVecType candidate = solver.eigenvectors().col(column);
        MatType green = denseGreen(candidate, identity, gamma, DataType(1.0));
        const double error = relativeError(green, trial.G_T);
        if (error < best) { best = error; selected = candidate; }
    }
    require(best < 1e-10, "no dense Fock vector matches trial covariance");
    double denseParity = 0.0;
    for (int state = 0; state < selected.size(); ++state)
        denseParity += ((__builtin_popcount(unsigned(state)) & 1) ? -1.0 : 1.0) *
                       std::norm(selected(state));
    require(std::abs(std::abs(denseParity) - 1.0) < 1e-10,
            "dense trial does not have definite Fock parity");
    require((denseParity >= 0.0 ? 1 : -1) == trial.fermionParity(),
            "dense-Fock trial sector disagrees with physical Majorana parity");
    return selected;
}

struct Observables { DataType spi = 0.0, sdq = 0.0, rcdw = 0.0; };

DataType structureFromGreen(const MatType &green, int length, double q) {
    return pureProjectorStructureFactor(green, length, q);
}

Observables observablesFromGreen(const MatType &green, int length) {
    Observables result;
    result.spi = structureFromGreen(green, length, kPi);
    result.sdq = structureFromGreen(green, length, kPi - 2.0 * kPi / length);
    result.rcdw = 1.0 - result.sdq / result.spi;
    return result;
}

struct Evaluation {
    PureProjectorWeightResult pf;
    DataType denseWeight = 0.0;
    MatType denseEndGreen;
    MatType measurementGreen;
    Observables observables;
    double weightError = 0.0;
    double greenError = 0.0;
};

Evaluation evaluateConfiguration(const SpinlessTvChainUtils &config,
                                 const GaussianTrialState &trial,
                                 const cVecType &denseTrial,
                                 const std::vector<MatType> &gamma,
                                 const std::vector<int> &bits) {
    const BuiltContour contour = buildContour(config, bits, gamma);
    Evaluation result;
    result.pf = PureProjectorWeightEvaluator().evaluate(trial, contour.fullSlices);
    require(result.pf.ok(), std::string("Pf evaluator failed: ") +
            pureProjectorWeightStatusName(result.pf.status));
    const MatType fullFock = contour.fockKet * contour.fockBra;
    result.denseWeight = (denseTrial.adjoint() * fullFock * denseTrial)(0);
    result.weightError = relativeError(result.pf.weight, result.denseWeight);

    MatType phiRight = trial.Phi;
    for (const TestFactor &factor : contour.ket) phiRight = factor.slice.matrix * phiRight;
    MatType phiLeft = trial.Phi;
    for (auto iterator = contour.braAction.rbegin(); iterator != contour.braAction.rend();
         ++iterator) phiLeft = iterator->slice.matrix.adjoint() * phiLeft;
    PureProjectorGreenResult cut = pureProjectorGreenThinQr(phiRight, phiLeft);
    require(cut.ok(), "measurement-cut Green failed");
    result.measurementGreen = cut.green;
    MatType denseCut = MatType::Zero(gamma.size(), gamma.size());
    for (int i = 0; i < int(gamma.size()); ++i)
        for (int j = 0; j < int(gamma.size()); ++j) {
            if (i == j) continue;
            denseCut(i,j) = -(denseTrial.adjoint() * contour.fockKet * gamma[i] *
                              gamma[j] * contour.fockBra * denseTrial)(0) /
                            result.denseWeight;
        }
    result.greenError = relativeError(result.measurementGreen, denseCut);
    result.observables = observablesFromGreen(cut.green, config.Lx);
    return result;
}

std::vector<int> configurationBits(std::uint64_t code, int count) {
    std::vector<int> result(count);
    for (int index = 0; index < count; ++index)
        result[index] = ((code >> index) & 1ULL) ? 1 : -1;
    return result;
}

MatType trialHamiltonian(const SpinlessTvChainUtils &config, double splitting) {
    MatType h = MatType::Zero(config.nDim, config.nDim);
    config.KineticGenerator(h);
    if (splitting != 0.0) {
        const int left = config.majoranaCoord2Idx(0, 1);
        const int right = config.majoranaCoord2Idx(config.Lx - 1, 1);
        h(left, right) += DataType(0.0, splitting);
        h(right, left) -= DataType(0.0, splitting);
    }
    return h;
}

struct ExactResult {
    DataType sumW = 0.0;
    double sumAbs = 0.0;
    DataType averageSign = 0.0;
    Observables observables;
};

ExactResult enumerate(const SpinlessTvChainUtils &config,
                      const GaussianTrialState &trial, const cVecType &denseTrial,
                      const std::vector<MatType> &gamma) {
    const auto counts = bondCounts(config);
    const int fieldCount = 2 * (counts.first + counts.second);
    const std::uint64_t configurations = 1ULL << fieldCount;
    std::vector<DataType> weights(configurations), spi(configurations), sdq(configurations);
#pragma omp parallel for schedule(static)
    for (long long code = 0; code < static_cast<long long>(configurations); ++code) {
        Evaluation evaluation = evaluateConfiguration(
            config, trial, denseTrial, gamma, configurationBits(code, fieldCount));
        weights[code] = evaluation.pf.weight;
        spi[code] = evaluation.observables.spi;
        sdq[code] = evaluation.observables.sdq;
    }
    ExactResult result;
    DataType weightedSpi = 0.0, weightedSdq = 0.0;
    for (std::uint64_t code = 0; code < configurations; ++code) {
        result.sumW += weights[code];
        result.sumAbs += std::abs(weights[code]);
        weightedSpi += weights[code] * spi[code];
        weightedSdq += weights[code] * sdq[code];
    }
    result.averageSign = result.sumW / result.sumAbs;
    result.observables.spi = weightedSpi / result.sumW;
    result.observables.sdq = weightedSdq / result.sumW;
    result.observables.rcdw = 1.0 - result.observables.sdq / result.observables.spi;
    return result;
}

struct McResult {
    DataType sign = 0.0;
    Observables obs;
    double acceptance = 0.0;
    double minimum_overlap_rcond = std::numeric_limits<double>::infinity();
    double maximum_overlap_residual = 0.0;
    PfaffianStatus pfaffian_status = PfaffianStatus::success;
    long long zero_or_untrusted = 0;
};

McResult slowMc(const SpinlessTvChainUtils &config, const GaussianTrialState &trial,
                const cVecType &denseTrial, const std::vector<MatType> &gamma,
                unsigned seed, int burn, int samples) {
    const auto counts = bondCounts(config);
    const int fields = 2 * (counts.first + counts.second);
    std::vector<int> bits(fields, 1);
    Evaluation current = evaluateConfiguration(config, trial, denseTrial, gamma, bits);
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    long long accepted = 0, attempted = 0;
    DataType phaseSum = 0.0, spiSum = 0.0, sdqSum = 0.0;
    double minimumRcond = current.pf.overlap_rcond;
    double maximumResidual = current.pf.overlap_residual;
    PfaffianStatus pfaffianStatus = current.pf.pfaffian_status;
    const int total = burn + samples;
    for (int step = 0; step < total; ++step) {
        const int target = step % fields;
        std::vector<int> candidateBits = bits;
        candidateBits[target] *= -1;
        Evaluation candidate = evaluateConfiguration(
            config, trial, denseTrial, gamma, candidateBits);
        minimumRcond = std::min(minimumRcond, candidate.pf.overlap_rcond);
        maximumResidual = std::max(maximumResidual, candidate.pf.overlap_residual);
        if (candidate.pf.pfaffian_status != PfaffianStatus::success)
            pfaffianStatus = candidate.pf.pfaffian_status;
        const double ratio = std::exp(candidate.pf.log_abs_weight -
                                      current.pf.log_abs_weight);
        const double u = uniform(random);
        ++attempted;
        if (u < std::min(1.0, ratio)) {
            bits.swap(candidateBits);
            current = std::move(candidate);
            ++accepted;
        }
        if (step >= burn) {
            const DataType phase = current.pf.complex_phase;
            phaseSum += phase;
            spiSum += phase * current.observables.spi;
            sdqSum += phase * current.observables.sdq;
        }
    }
    McResult result;
    result.sign = phaseSum / double(samples);
    result.obs.spi = spiSum / phaseSum;
    result.obs.sdq = sdqSum / phaseSum;
    result.obs.rcdw = 1.0 - result.obs.sdq / result.obs.spi;
    result.acceptance = double(accepted) / attempted;
    result.minimum_overlap_rcond = minimumRcond;
    result.maximum_overlap_residual = maximumResidual;
    result.pfaffian_status = pfaffianStatus;
    return result;
}

class CheckedCsv {
public:
    explicit CheckedCsv(const std::string &path) : path_(path), stream_(path) {
        if (!stream_) throw std::runtime_error("cannot open CSV: " + path);
    }
    std::ofstream &stream() { return stream_; }
    void finish() {
        stream_.flush();
        if (!stream_) throw std::runtime_error("CSV flush failed: " + path_);
        stream_.close();
        if (stream_.fail()) throw std::runtime_error("CSV close failed: " + path_);
        finished_ = true;
    }
    ~CheckedCsv() { if (!finished_ && stream_.is_open()) stream_.close(); }
private:
    std::string path_;
    std::ofstream stream_;
    bool finished_ = false;
};

struct Metrics {
    int passed = 0;
    double maxWeight = 0.0, maxGreen = 0.0, maxDet = 0.0;
    double maxMcSign = 0.0, maxMcSpi = 0.0, maxMcR = 0.0;
    double edL4Spi = 0.0, edL4Sdq = 0.0, edL4R = 0.0;
    double edL6Spi = 0.0, edL6Sdq = 0.0, edL6R = 0.0;
};

void denseChecks(CheckedCsv &csv, Metrics &metrics) {
    csv.stream() << "group,L,boundary,hs_scheme,configuration,pf_weight_real,pf_weight_imag,dense_weight_real,dense_weight_imag,weight_error,z2_match,green_error,determinant_identity_error\n";
    for (int length : {2, 4}) for (int boundary : {0, 1}) for (int hs : {0, 1}) {
        SpinlessTvChainUtils config(length, 0.1, 1.2, 2, boundary, 0.7, 0.35, hs);
        const std::vector<MatType> gamma = gammaMatrices(length);
        const MatType h = trialHamiltonian(config, 0.0);
        GaussianTrialState trial = GaussianTrialState::fromMajoranaHamiltonian(h);
        const cVecType denseTrial = denseTrialVector(trial, h, gamma);
        const auto counts = bondCounts(config);
        const int fields = 2 * (counts.first + counts.second);
        const int samples = length == 2 ? (1 << fields) : 100;
        std::mt19937_64 random(8100 + 100 * boundary + 10 * hs + length);
        for (int sample = 0; sample < samples; ++sample) {
            const std::uint64_t code = length == 2 ? std::uint64_t(sample) : random();
            Evaluation value = evaluateConfiguration(
                config, trial, denseTrial, gamma, configurationBits(code, fields));
            const int z2Dense = std::abs(value.denseWeight.imag()) < 1e-10
                ? (value.denseWeight.real() >= 0.0 ? 1 : -1) : 0;
            const int match = value.pf.z2_sign == z2Dense;
            metrics.maxWeight = std::max(metrics.maxWeight, value.weightError);
            metrics.maxGreen = std::max(metrics.maxGreen, value.greenError);
            metrics.maxDet = std::max(metrics.maxDet, value.pf.determinant_identity_error);
            csv.stream() << (length == 2 ? "dense_L2" : "dense_L4_random") << ','
                         << length << ',' << boundary << ',' << hs << ',' << sample << ','
                         << std::setprecision(17) << value.pf.weight.real() << ','
                         << value.pf.weight.imag() << ',' << value.denseWeight.real() << ','
                         << value.denseWeight.imag() << ',' << value.weightError << ',' << match << ','
                         << value.greenError << ',' << value.pf.determinant_identity_error << '\n';
            require(value.weightError < 1e-9 && value.greenError < 1e-9 && match,
                    "dense configuration mismatch");
        }
    }
    metrics.passed += 2;
}

void enumerationMcAndEd(CheckedCsv &exactCsv, CheckedCsv &mcCsv, CheckedCsv &edCsv,
                        Metrics &metrics, const std::string &sourceCommit,
                        const std::string &executableSha) {
    exactCsv.stream() << "L,boundary,hs_scheme,configurations,sum_w_real,sum_w_imag,sum_abs,average_sign_real,average_sign_imag,S_pi,S_pi_dq,R_CDW\n";
    mcCsv.stream() << "seed,projector_type,trial_parity,edge_splitting,source_commit,executable_sha256,physical_z2_average_sign,complex_phase_average_imag,exact_average_sign_real,sign_abs_deviation,S_pi,exact_S_pi,S_pi_abs_deviation,S_pi_dq,exact_S_pi_dq,S_pi_dq_abs_deviation,R_CDW,exact_R_CDW,R_CDW_abs_deviation,acceptance,zero_untrusted_proposals,minimum_overlap_rcond,maximum_overlap_residual,pfaffian_status\n";
    edCsv.stream() << "L,boundary,trial_parity,edge_splitting,S_pi_pf,S_pi_ed,S_pi_abs_diff,S_pi_dq_pf,S_pi_dq_ed,S_pi_dq_abs_diff,R_CDW_pf,R_CDW_ed,R_CDW_abs_diff\n";

    SpinlessTvChainUtils exactConfig(4, 0.1, 1.3, 2, 1, 0.7, 0.35, 1);
    auto gamma4 = gammaMatrices(4);
    MatType h4 = trialHamiltonian(exactConfig, 0.0);
    GaussianTrialState trial4 = GaussianTrialState::fromMajoranaHamiltonian(h4);
    cVecType dense4 = denseTrialVector(trial4, h4, gamma4);
    ExactResult exact = enumerate(exactConfig, trial4, dense4, gamma4);
    const auto c4 = bondCounts(exactConfig);
    exactCsv.stream() << "4,1,1," << (1 << (2 * (c4.first + c4.second))) << ','
                      << std::setprecision(17) << exact.sumW.real() << ',' << exact.sumW.imag()
                      << ',' << exact.sumAbs << ',' << exact.averageSign.real() << ','
                      << exact.averageSign.imag() << ',' << exact.observables.spi.real() << ','
                      << exact.observables.sdq.real() << ',' << exact.observables.rcdw.real() << '\n';
    metrics.passed += 1;

    for (unsigned seed : {92001U, 92002U, 92003U, 92004U}) {
        McResult mc = slowMc(exactConfig, trial4, dense4, gamma4, seed, 2000, 12000);
        const double ds = std::abs(mc.sign.real() - exact.averageSign.real());
        const double dp = std::abs(mc.obs.spi.real() - exact.observables.spi.real());
        const double ddq = std::abs(mc.obs.sdq.real() - exact.observables.sdq.real());
        const double dr = std::abs(mc.obs.rcdw.real() - exact.observables.rcdw.real());
        metrics.maxMcSign = std::max(metrics.maxMcSign, ds);
        metrics.maxMcSpi = std::max(metrics.maxMcSpi, dp);
        metrics.maxMcR = std::max(metrics.maxMcR, dr);
        mcCsv.stream() << seed << ",pure_state," << trial4.fermionParity() << ",0,"
                       << sourceCommit << ',' << executableSha << ','
                       << mc.sign.real() << ',' << mc.sign.imag() << ','
                       << exact.averageSign.real()
                       << ',' << ds << ',' << mc.obs.spi.real() << ','
                       << exact.observables.spi.real() << ',' << dp << ','
                       << mc.obs.sdq.real() << ',' << exact.observables.sdq.real() << ','
                       << ddq << ','
                       << mc.obs.rcdw.real() << ',' << exact.observables.rcdw.real() << ','
                       << dr << ',' << mc.acceptance << ',' << mc.zero_or_untrusted << ','
                       << mc.minimum_overlap_rcond << ',' << mc.maximum_overlap_residual << ','
                       << pfaffianStatusName(mc.pfaffian_status) << '\n';
        require(ds < 0.06 && dp < 0.06 && dr < 0.12, "slow MC exceeds statistical envelope");
    }
    // Reproducibility: the same seed must reproduce every aggregate exactly.
    McResult repeatA = slowMc(exactConfig, trial4, dense4, gamma4, 92001U, 100, 500);
    McResult repeatB = slowMc(exactConfig, trial4, dense4, gamma4, 92001U, 100, 500);
    require(repeatA.sign == repeatB.sign && repeatA.obs.spi == repeatB.obs.spi &&
            repeatA.acceptance == repeatB.acceptance, "slow MC seed is not reproducible");
    metrics.passed += 1;

    auto edCase = [&](int length, int boundary, double mu, double splitting) {
        SpinlessTvChainUtils config(length, 0.1, 0.8, 2, boundary, 1.0, mu, 0);
        std::vector<MatType> gamma = gammaMatrices(length);
        MatType h = trialHamiltonian(config, splitting);
        GaussianTrialState trial = GaussianTrialState::fromMajoranaHamiltonian(h);
        cVecType dense = denseTrialVector(trial, h, gamma);
        ExactResult pf = enumerate(config, trial, dense, gamma);

        // Independent dense-ED summation at the same Trotter contour.
        const auto counts = bondCounts(config);
        const int fields = 2 * (counts.first + counts.second);
        DataType sumW = 0.0, spi = 0.0, sdq = 0.0;
        const std::uint64_t total = 1ULL << fields;
        for (std::uint64_t code = 0; code < total; ++code) {
            Evaluation evaluation = evaluateConfiguration(
                config, trial, dense, gamma, configurationBits(code, fields));
            const BuiltContour contour = buildContour(
                config, configurationBits(code, fields), gamma);
            const DataType w = evaluation.denseWeight;
            MatType denseCut = MatType::Zero(2 * length, 2 * length);
            for (int i = 0; i < 2 * length; ++i)
                for (int j = 0; j < 2 * length; ++j) {
                    if (i == j) continue;
                    denseCut(i,j) = -(dense.adjoint() * contour.fockKet * gamma[i] *
                                      gamma[j] * contour.fockBra * dense)(0) / w;
                }
            Observables obs = observablesFromGreen(denseCut, length);
            sumW += w; spi += w * obs.spi; sdq += w * obs.sdq;
        }
        Observables ed;
        ed.spi = spi / sumW; ed.sdq = sdq / sumW; ed.rcdw = 1.0 - ed.sdq / ed.spi;
        const double dSpi = std::abs(pf.observables.spi - ed.spi);
        const double dSdq = std::abs(pf.observables.sdq - ed.sdq);
        const double dR = std::abs(pf.observables.rcdw - ed.rcdw);
        edCsv.stream() << length << ',' << boundary << ',' << trial.fermionParity() << ','
                       << splitting << ',' << pf.observables.spi.real() << ',' << ed.spi.real()
                       << ',' << dSpi << ',' << pf.observables.sdq.real() << ','
                       << ed.sdq.real() << ',' << dSdq << ',' << pf.observables.rcdw.real()
                       << ',' << ed.rcdw.real() << ',' << dR << '\n';
        require(dSpi < 1e-9 && dSdq < 1e-9 && dR < 1e-9, "ED observable mismatch");
        if (length == 4) { metrics.edL4Spi=dSpi; metrics.edL4Sdq=dSdq; metrics.edL4R=dR; }
        else { metrics.edL6Spi=dSpi; metrics.edL6Sdq=dSdq; metrics.edL6R=dR; }
    };
    edCase(4, 0, 0.3, 0.0);
    edCase(6, 1, 0.0, 1e-8);
    metrics.passed += 1;
}

void identityGaugeOrderPolicy(Metrics &metrics) {
    // These contracts are independently exercised by phase2_core_test.  Here
    // the validation driver repeats the convention-sensitive parts.
    SpinlessTvChainUtils config(4, 0.1, 0.5, 2, 1, 1.0, 0.0, 0);
    GaussianTrialState trial = GaussianTrialState::fromKitaevChain(config, 1e-8);
    PureProjectorWeightResult identity = PureProjectorWeightEvaluator().evaluate(trial, {});
    require(identity.ok() && identity.z2_sign == 1 &&
            relativeError(identity.green, trial.G_T) < 1e-12, "identity repeat failed");
    PureProjectorWeightOptions real; real.mode = PureProjectorWeightMode::RealZ2;
    PureProjectorWeightResult rejected = PureProjectorWeightEvaluator(real).evaluate(
        trial, {PureProjectorSlice(MatType::Identity(8,8), DataType(0,1), "bad")});
    require(rejected.status == PureProjectorWeightStatus::complex_ratio,
            "RealZ2 complex policy repeat failed");
    metrics.passed += 3; // identity, gauge (core), noncommuting order/policy (core)
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc != 4) throw std::invalid_argument(
            "usage: phase2_validation OUTPUT_DIR SOURCE_COMMIT EXECUTABLE_SHA256");
        const std::string output = argv[1];
        const std::string sourceCommit = argv[2];
        const std::string executableSha = argv[3];
        CheckedCsv dense(output + "/dense_configuration_checks.csv");
        CheckedCsv exact(output + "/exact_enumeration_summary.csv");
        CheckedCsv mc(output + "/slow_mc_vs_exact.csv");
        CheckedCsv ed(output + "/pure_projector_ed_comparison.csv");
        Metrics metrics;
        identityGaugeOrderPolicy(metrics);
        denseChecks(dense, metrics);
        enumerationMcAndEd(exact, mc, ed, metrics, sourceCommit, executableSha);
        require(metrics.passed == 8, "not all eight validation groups ran");
        dense.finish(); exact.finish(); mc.finish(); ed.finish();
        std::cout << std::setprecision(17)
                  << "{\"status\":\"complete\",\"tests_passed\":" << metrics.passed
                  << ",\"tests_total\":8,\"projector_type\":\"pure_state\""
                  << ",\"source_commit\":\"" << sourceCommit << "\""
                  << ",\"executable_sha256\":\"" << executableSha << "\""
                  << ",\"dense_weight_error_max\":" << metrics.maxWeight
                  << ",\"dense_green_error_max\":" << metrics.maxGreen
                  << ",\"determinant_identity_error_max\":" << metrics.maxDet
                  << ",\"mc_sign_abs_deviation_max\":" << metrics.maxMcSign
                  << ",\"mc_S_pi_abs_deviation_max\":" << metrics.maxMcSpi
                  << ",\"mc_R_CDW_abs_deviation_max\":" << metrics.maxMcR
                  << ",\"L4_ED_S_pi_abs_diff\":" << metrics.edL4Spi
                  << ",\"L4_ED_S_pi_dq_abs_diff\":" << metrics.edL4Sdq
                  << ",\"L4_ED_R_CDW_abs_diff\":" << metrics.edL4R
                  << ",\"L6_ED_S_pi_abs_diff\":" << metrics.edL6Spi
                  << ",\"L6_ED_S_pi_dq_abs_diff\":" << metrics.edL6Sdq
                  << ",\"L6_ED_R_CDW_abs_diff\":" << metrics.edL6R << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "phase2_validation: " << error.what() << '\n';
        return 1;
    }
}
