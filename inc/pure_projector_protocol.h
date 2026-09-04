#ifndef PURE_PROJECTOR_PROTOCOL_H
#define PURE_PROJECTOR_PROTOCOL_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

enum class PureImaginaryTimeProtocolKind { Constant, SuddenQuench, LinearRamp };

class PureImaginaryTimeProtocol {
public:
    static PureImaginaryTimeProtocol constant(double value, double tauFinal,
                                               double deltaTau) {
        return PureImaginaryTimeProtocol(
            PureImaginaryTimeProtocolKind::Constant, value, value, 0.0,
            tauFinal, deltaTau);
    }

    static PureImaginaryTimeProtocol suddenQuench(double initialValue,
                                                  double finalValue,
                                                  double tauFinal,
                                                  double deltaTau) {
        return PureImaginaryTimeProtocol(
            PureImaginaryTimeProtocolKind::SuddenQuench, initialValue,
            finalValue, 0.0, tauFinal, deltaTau);
    }

    static PureImaginaryTimeProtocol linearRamp(double initialValue,
                                                double rate,
                                                double tauFinal,
                                                double deltaTau) {
        return PureImaginaryTimeProtocol(
            PureImaginaryTimeProtocolKind::LinearRamp, initialValue,
            initialValue + rate * tauFinal, rate, tauFinal, deltaTau);
    }

    PureImaginaryTimeProtocolKind kind() const { return kind_; }
    double initialValue() const { return initial_value_; }
    double finalValue() const { return final_value_; }
    double rate() const { return rate_; }
    double tauFinal() const { return tau_final_; }
    double deltaTau() const { return delta_tau_; }
    int slices() const { return int(midpoint_values_.size()); }
    double midpointValue(int slice) const {
        if (slice < 0 || slice >= slices())
            throw std::out_of_range("imaginary-time protocol slice out of range");
        return midpoint_values_[slice];
    }
    const std::vector<double> &midpointValues() const { return midpoint_values_; }

    const char *name() const {
        switch (kind_) {
        case PureImaginaryTimeProtocolKind::Constant: return "constant";
        case PureImaginaryTimeProtocolKind::SuddenQuench: return "sudden_quench";
        case PureImaginaryTimeProtocolKind::LinearRamp: return "linear_ramp";
        }
        return "unknown";
    }

private:
    PureImaginaryTimeProtocolKind kind_;
    double initial_value_ = 0.0;
    double final_value_ = 0.0;
    double rate_ = 0.0;
    double tau_final_ = 0.0;
    double delta_tau_ = 0.0;
    std::vector<double> midpoint_values_;

    PureImaginaryTimeProtocol(PureImaginaryTimeProtocolKind kind,
                              double initialValue, double finalValue,
                              double rate, double tauFinal, double deltaTau)
        : kind_(kind), initial_value_(initialValue), final_value_(finalValue),
          rate_(rate), tau_final_(tauFinal), delta_tau_(deltaTau) {
        if (!(deltaTau > 0.0) || !(tauFinal > 0.0) || initialValue < 0.0 ||
            finalValue < 0.0 || !std::isfinite(initialValue) ||
            !std::isfinite(finalValue) || !std::isfinite(rate) ||
            !std::isfinite(tauFinal) || !std::isfinite(deltaTau))
            throw std::invalid_argument("invalid imaginary-time protocol");
        const double count = tauFinal / deltaTau;
        const long long rounded = std::llround(count);
        if (rounded <= 0 ||
            std::abs(count - rounded) > 1e-10 * std::max(1.0, std::abs(count)))
            throw std::invalid_argument("tau_f / Delta_tau must be a positive integer");
        midpoint_values_.reserve(std::size_t(rounded));
        for (long long slice = 0; slice < rounded; ++slice) {
            double value = finalValue;
            if (kind == PureImaginaryTimeProtocolKind::Constant)
                value = initialValue;
            else if (kind == PureImaginaryTimeProtocolKind::LinearRamp)
                value = initialValue + rate * (slice + 0.5) * deltaTau;
            if (value < 0.0 || !std::isfinite(value))
                throw std::invalid_argument("protocol has invalid midpoint V(tau)");
            midpoint_values_.push_back(value);
        }
    }
};

#endif
