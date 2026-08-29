#ifndef PURE_PROJECTOR_STATIC_H
#define PURE_PROJECTOR_STATIC_H

#include "pure_projector_weight.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct PureStaticMeasurement {
    PureProjectorStatus status = PureProjectorStatus::invalid_dimension;
    MatType green;
    double overlap_rcond = 0.0;
    double solve_residual = std::numeric_limits<double>::infinity();
    bool ok() const { return status == PureProjectorStatus::success; }
};

class PureStaticProjectorContour {
public:
    PureStaticProjectorContour(const GaussianTrialState &trial, double theta, double dt,
                               int trialParity, double edgeSplitting)
        : trial_(trial), theta_(theta), dt_(dt), trial_parity_(trialParity),
          edge_splitting_(edgeSplitting) {
        if (!(theta > 0.0) || !(dt > 0.0) ||
            std::abs(theta / dt - std::round(theta / dt)) > 1e-10)
            throw std::invalid_argument("pure contour requires theta/dt positive integer");
        if (trialParity != -1 && trialParity != 1)
            throw std::invalid_argument("pure contour requires explicit trial parity");
    }

    void setKetActionOrder(std::vector<PureProjectorSlice> slices) {
        ket_ = std::move(slices);
    }
    void setBraProtocolOrder(std::vector<PureProjectorSlice> slices) {
        bra_protocol_ = std::move(slices);
    }

    std::vector<std::string> flattenedActionOrder() const {
        std::vector<std::string> result;
        for (const auto &slice : ket_) result.push_back("ket:" + slice.label);
        for (auto iterator = bra_protocol_.rbegin(); iterator != bra_protocol_.rend(); ++iterator)
            result.push_back("bra_adjoint:" + iterator->label);
        return result;
    }

    PureStaticMeasurement measurementGreen() const {
        MatType right = trial_.Phi;
        for (const PureProjectorSlice &slice : ket_) right = slice.matrix * right;
        MatType left = trial_.Phi;
        for (auto iterator = bra_protocol_.rbegin(); iterator != bra_protocol_.rend(); ++iterator)
            left = iterator->matrix.adjoint() * left;
        const PureProjectorGreenResult green = pureProjectorGreenThinQr(right, left);
        PureStaticMeasurement result;
        result.status = green.status;
        result.green = green.green;
        result.overlap_rcond = green.overlap_rcond;
        result.solve_residual = green.solve_residual;
        return result;
    }

    int measurementCut() const { return int(ket_.size()); }
    int trialParity() const { return trial_parity_; }
    double edgeSplitting() const { return edge_splitting_; }
    double theta() const { return theta_; }
    double dt() const { return dt_; }

private:
    GaussianTrialState trial_;
    double theta_;
    double dt_;
    int trial_parity_;
    double edge_splitting_;
    std::vector<PureProjectorSlice> ket_;
    std::vector<PureProjectorSlice> bra_protocol_;
};

struct PureProjectorConfiguration {
    std::vector<PureProjectorSlice> slices;
    std::vector<int> hs_fields;
};

struct PureProjectorProposal {
    int index = -1;
    int candidate_field = 0;
    double uniform = 0.0;
};

struct PureProjectorWalkerSnapshot { std::uint64_t state_hash = 0; };

struct PureProjectorProposalResult {
    int proposal_index = -1;
    double uniform = 0.0;
    bool accepted = false;
    bool zero_or_untrusted = false;
    std::uint64_t predecision_hash = 0;
};

class PureProjectorSlowWalker {
public:
    PureProjectorSlowWalker(const GaussianTrialState &trial,
                            PureProjectorConfiguration configuration,
                            unsigned long long seed,
                            PureProjectorWeightOptions options = PureProjectorWeightOptions())
        : trial_(trial), configuration_(std::move(configuration)), seed_(seed),
          evaluator_(options) {
        if (configuration_.slices.size() != configuration_.hs_fields.size())
            throw std::invalid_argument("slow walker requires one field per test slice");
        current_ = evaluator_.evaluate(trial_, configuration_.slices);
        if (!current_.ok()) throw std::runtime_error("invalid initial slow-walker weight");
    }

    PureProjectorWalkerSnapshot snapshot() const { return {stateHash(configuration_)}; }
    const PureProjectorWeightResult &currentWeight() const { return current_; }
    const PureProjectorConfiguration &configuration() const { return configuration_; }

    PureProjectorProposalResult propose(const PureProjectorProposal &proposal) {
        PureProjectorProposalResult result;
        result.proposal_index = proposal.index;
        result.uniform = proposal.uniform;
        result.predecision_hash = stateHash(configuration_);
        if (proposal.index < 0 || proposal.index >= int(configuration_.slices.size()) ||
            (proposal.candidate_field != -1 && proposal.candidate_field != 1) ||
            !(proposal.uniform >= 0.0 && proposal.uniform < 1.0)) {
            result.zero_or_untrusted = true;
            return result;
        }

        PureProjectorConfiguration candidate = configuration_;
        if (candidate.hs_fields[proposal.index] != proposal.candidate_field) {
            candidate.hs_fields[proposal.index] = proposal.candidate_field;
            candidate.slices[proposal.index].matrix =
                candidate.slices[proposal.index].matrix.adjoint().eval();
        }
        const PureProjectorWeightResult candidateWeight =
            evaluator_.evaluate(trial_, candidate.slices);
        if (!candidateWeight.ok()) {
            result.zero_or_untrusted = true;
            ++zero_or_untrusted_;
            return result;
        }
        const double delta = candidateWeight.log_abs_weight - current_.log_abs_weight;
        const double probability = delta >= 0.0 ? 1.0 : std::exp(delta);
        result.accepted = proposal.uniform < probability;
        ++attempted_;
        if (result.accepted) {
            configuration_ = std::move(candidate);
            current_ = candidateWeight;
            ++accepted_;
        }
        return result;
    }

    long long attempted() const { return attempted_; }
    long long accepted() const { return accepted_; }
    long long zeroOrUntrusted() const { return zero_or_untrusted_; }
    unsigned long long seed() const { return seed_; }

private:
    GaussianTrialState trial_;
    PureProjectorConfiguration configuration_;
    unsigned long long seed_;
    PureProjectorWeightEvaluator evaluator_;
    PureProjectorWeightResult current_;
    long long attempted_ = 0;
    long long accepted_ = 0;
    long long zero_or_untrusted_ = 0;

    static std::uint64_t stateHash(const PureProjectorConfiguration &configuration) {
        std::uint64_t hash = 1469598103934665603ULL;
        auto mix = [&hash](std::uint64_t value) {
            hash ^= value;
            hash *= 1099511628211ULL;
        };
        for (int field : configuration.hs_fields) mix(std::uint64_t(std::int64_t(field)));
        for (const auto &slice : configuration.slices) {
            for (int col = 0; col < slice.matrix.cols(); ++col)
                for (int row = 0; row < slice.matrix.rows(); ++row) {
                    std::uint64_t realBits = 0, imagBits = 0;
                    const double real = slice.matrix(row, col).real();
                    const double imag = slice.matrix(row, col).imag();
                    std::memcpy(&realBits, &real, sizeof(double));
                    std::memcpy(&imagBits, &imag, sizeof(double));
                    mix(realBits); mix(imagBits);
                }
        }
        return hash;
    }
};

#endif
