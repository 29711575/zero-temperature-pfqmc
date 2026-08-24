#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "../projector_kitaev/projector_contour.h"
#include "pfqmc.h"

struct StaticWalker : Spinless_tV {
    int center = 0;
    int nt = 0;
    int np = 0;

    StaticWalker(const SpinlessTvChainUtils *config, rdGenerator *random,
                 double theta, double beta) {
        build_projector_static_contour(*this, config, random, theta, beta, center, nt, np);
    }
};

MatType centerNoUpdate(PfQMC &qmc, int center) {
    MatType g = qmc.g;
    MatType temporary;
    for (int index = 0; index < center; ++index) qmc.op_array[index]->left_propagate(g, temporary);
    return g;
}

double operatorDifference(const Spinless_tV &left, const Spinless_tV &right) {
    if (left.op_array.size() != right.op_array.size()) return INFINITY;
    double difference = 0.0;
    for (size_t index = 0; index < left.op_array.size(); ++index) {
        auto *leftDense = dynamic_cast<DenseOperator *>(left.op_array[index]);
        auto *rightDense = dynamic_cast<DenseOperator *>(right.op_array[index]);
        auto *leftInteraction = dynamic_cast<SpinlessVOperator *>(left.op_array[index]);
        auto *rightInteraction = dynamic_cast<SpinlessVOperator *>(right.op_array[index]);
        if (bool(leftDense) != bool(rightDense) ||
            bool(leftInteraction) != bool(rightInteraction)) return INFINITY;
        if (leftDense) {
            difference = std::max(
                difference, (leftDense->mat - rightDense->mat).cwiseAbs().maxCoeff());
        }
        if (leftInteraction) {
            if (leftInteraction->bondType != rightInteraction->bondType ||
                leftInteraction->localV != rightInteraction->localV ||
                leftInteraction->s->size() != rightInteraction->s->size()) {
                return INFINITY;
            }
            for (int field = 0; field < leftInteraction->s->size(); ++field) {
                if ((*leftInteraction->s)(field) != (*rightInteraction->s)(field)) return INFINITY;
            }
            difference = std::max(
                difference, (leftInteraction->B - rightInteraction->B).cwiseAbs().maxCoeff());
        }
    }
    return difference;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: static_contour_compare V\n";
        return 2;
    }
    const double interaction = std::stod(argv[1]);
    const int L = 6;
    const double dt = .1;
    const double theta = 2;
    const double beta = 4;
    SpinlessTvChainUtils projectorConfig(L, dt, interaction, 40, 1, 1, 0, 0);
    SpinlessTvChainUtils drivenConfig(L, dt, interaction, 40, 1, 1, 0, 0);
    rdGenerator projectorRandom(777);
    rdGenerator drivenRandom(777);
    StaticWalker projector(&projectorConfig, &projectorRandom, theta, beta);
    StaticWalker drivenStatic(&drivenConfig, &drivenRandom, theta, beta);
    PfQMC projectorQmc(&projector, 10);
    PfQMC drivenQmc(&drivenStatic, 10);
    const MatType projectorGreen = centerNoUpdate(projectorQmc, projector.center);
    const MatType drivenGreen = centerNoUpdate(drivenQmc, drivenStatic.center);
    const double projectorSpi = -projectorConfig.StructureFactorCDW(projectorGreen).real();
    const double drivenSpi = -drivenConfig.StructureFactorCDW(drivenGreen).real();
    const double projectorSpidq =
        -projectorConfig.StructureFactorCDWOffset(projectorGreen).real();
    const double drivenSpidq = -drivenConfig.StructureFactorCDWOffset(drivenGreen).real();
    std::cout << std::setprecision(17)
              << "{\"V\":" << interaction
              << ",\"trial_slices\":" << projector.nt
              << ",\"physical_slices\":" << projector.np
              << ",\"operator_count\":" << projector.op_array.size()
              << ",\"center_boundary\":" << projector.center
              << ",\"operator_max_abs_diff\":"
              << operatorDifference(projector, drivenStatic)
              << ",\"weight_sign_abs_diff\":"
              << std::abs(projectorQmc.getSignRaw() - drivenQmc.getSignRaw())
              << ",\"green_max_abs_diff\":"
              << (projectorGreen - drivenGreen).cwiseAbs().maxCoeff()
              << ",\"S_pi_abs_diff\":" << std::abs(projectorSpi - drivenSpi)
              << ",\"S_pi_dq_abs_diff\":" << std::abs(projectorSpidq - drivenSpidq)
              << ",\"R_abs_diff\":"
              << std::abs((1.0 - projectorSpidq / projectorSpi) -
                          (1.0 - drivenSpidq / drivenSpi))
              << "}\n";
}
