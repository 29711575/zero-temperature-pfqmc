#ifndef PROJECTOR_KITAEV_CONTOUR_H
#define PROJECTOR_KITAEV_CONTOUR_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "kitaevChain.h"

inline void build_projector_static_contour(
    Spinless_tV &walker, const SpinlessTvChainUtils *config, rdGenerator *random,
    double theta, double beta, int &center, int &trialSlices, int &physicalSlices) {
    if (config == nullptr || random == nullptr) {
        throw std::invalid_argument("projector contour requires config and RNG");
    }
    if (!std::isfinite(config->dt) || !std::isfinite(theta) || !std::isfinite(beta) ||
        config->dt <= 0 || theta <= 0 || beta <= 0) {
        throw std::invalid_argument("projector contour requires dt, theta, beta_trial > 0");
    }
    if (!walker.op_array.empty()) {
        throw std::logic_error("projector contour requires an empty operator array");
    }
    const double thetaSlices = theta / config->dt;
    const long long roundedThetaSlices = std::llround(thetaSlices);
    if (roundedThetaSlices <= 0 ||
        std::abs(thetaSlices - roundedThetaSlices) > 1e-10) {
        throw std::invalid_argument("theta/dt must be a positive integer");
    }

    walker.nDim = config->nDim;
    center = -1;
    trialSlices = 0;
    physicalSlices = 2 * int(roundedThetaSlices);
    int bonds[2];
    if (config->boundaryType == 0) {
        bonds[0] = (config->Lx + 1) / 2;
        bonds[1] = config->Lx / 2;
    } else {
        bonds[0] = config->Lx / 2;
        bonds[1] = (config->Lx - 1) / 2;
    }

    MatType kinetic(config->nDim, config->nDim);
    config->KineticGenerator(kinetic);
    for (double remaining = beta; remaining > 1e-12;) {
        const double step = std::min(config->dt, remaining);
        MatType exponent = kinetic;
        MatType matrix = expm(exponent, -step);
        exponent = step * kinetic;
        walker.op_array.push_back(new DenseOperator(matrix, signOfHamiltonian(exponent)));
        ++trialSlices;
        remaining -= step;
    }
    if (trialSlices <= 0) throw std::logic_error("projector contour has no trial slices");

    MatType exponent = kinetic;
    MatType halfKinetic = expm(exponent, -config->dt / 2);
    exponent = (config->dt / 2) * kinetic;
    const DataType kineticSign = signOfHamiltonian(exponent);
    for (int slice = 0; slice < physicalSlices; ++slice) {
        walker.op_array.push_back(new DenseOperator(halfKinetic, kineticSign));
        for (int bond = 0; bond < 2; ++bond) {
            auto *aux = new iVecType(bonds[bond]);
            for (int index = 0; index < aux->size(); ++index) (*aux)(index) = random->rdZ2();
            walker.op_array.push_back(new SpinlessVOperator(config, aux, bond, random));
        }
        walker.op_array.push_back(new DenseOperator(halfKinetic, kineticSign));
        if (slice + 1 == physicalSlices / 2) center = int(walker.op_array.size());
    }

    const int expectedOperatorCount = trialSlices + 4 * physicalSlices;
    const int expectedCenter = trialSlices + 2 * physicalSlices;
    if (int(walker.op_array.size()) != expectedOperatorCount ||
        center != expectedCenter || center < 1 || center >= expectedOperatorCount) {
        throw std::logic_error("projector contour slice/operator count mismatch");
    }
}

#endif
