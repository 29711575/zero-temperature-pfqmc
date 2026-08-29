#include "test_common.h"

#include <iostream>

namespace {

DataType directWickInteractionBond(const SpinlessTvChainUtils &config,
                                   const MatType &g, int site) {
    const int next = (site + 1) % config.Lx;
    const int i1 = config.majoranaCoord2Idx(site, 0);
    const int i2 = config.majoranaCoord2Idx(site, 1);
    const int j1 = config.majoranaCoord2Idx(next, 0);
    const int j2 = config.majoranaCoord2Idx(next, 1);
    return 0.25 * config.V *
        (g(i1, j1) * g(i2, j2) +
         g(i1, j2) * g(j1, i2) -
         g(i1, i2) * g(j1, j2));
}

DataType directWickInteraction(const SpinlessTvChainUtils &config,
                               const MatType &g) {
    DataType result = 0.0;
    const int bonds = config.boundaryType == 0 ? config.Lx : config.Lx - 1;
    for (int site = 0; site < bonds; ++site) {
        result += directWickInteractionBond(config, g, site);
    }
    return result;
}

MatType deterministicGreen(int dimension) {
    MatType g = MatType::Identity(dimension, dimension);
    for (int row = 0; row < dimension; ++row) {
        for (int col = row + 1; col < dimension; ++col) {
            const DataType value(0.013 * (row + 1) - 0.007 * (col + 1),
                                 0.005 * (row + col + 2));
            g(row, col) = value;
            g(col, row) = -value;
        }
    }
    return g;
}

double relativeScalarError(DataType actual, DataType expected) {
    return std::abs(actual - expected) /
        std::max(std::abs(expected), std::numeric_limits<double>::min());
}

}  // namespace

int main(int argc, char **argv) try {
    if (argc != 2) {
        throw std::runtime_error("usage: pbc_structural_guards_driver output.csv");
    }
    std::ofstream out(argv[1]);
    if (!out) throw std::runtime_error("cannot open output CSV");
    out << "test,L,boundary,status,value,tolerance,detail\n"
        << std::setprecision(17);
    int failures = 0;

    bool oddPbcRejected = false;
    try {
        SpinlessTvChainUtils invalid(3, 0.1, 2.0, 2, 0, 1.0, 0.0, 0);
        (void)invalid;
    } catch (const std::invalid_argument &error) {
        oddPbcRejected = std::string(error.what()).find("odd-L PBC") !=
                         std::string::npos;
    }
    failures += !oddPbcRejected;
    out << "odd_L_PBC_construction,3,PBC,"
        << (oddPbcRejected ? "PASS" : "FAIL") << ','
        << (oddPbcRejected ? 1 : 0) << ",1,shared_config_rejection\n";

    bool oddObcConstructed = false;
    try {
        SpinlessTvChainUtils config(3, 0.1, 2.0, 2, 1, 1.0, 0.0, 0);
        rdGenerator random(630031);
        Chain_tV walker(&config, &random);
        oddObcConstructed = walker.nBond[0] == 1 && walker.nBond[1] == 1 &&
                            !walker.op_array.empty();
    } catch (...) {
        oddObcConstructed = false;
    }
    failures += !oddObcConstructed;
    out << "odd_L_OBC_construction,3,OBC,"
        << (oddObcConstructed ? "PASS" : "FAIL") << ','
        << (oddObcConstructed ? 1 : 0) << ",1,finite_T_walker_constructed\n";

    double maxInverseResidual = 0.0;
    bool evenPbcConstructed = false;
    try {
        SpinlessTvChainUtils config(6, 0.1, 2.0, 2, 0, 1.0, 0.0, 0);
        rdGenerator random(630061);
        Chain_tV walker(&config, &random);
        evenPbcConstructed = true;
        const MatType identity = MatType::Identity(config.nDim, config.nDim);
        for (Operator *base : walker.op_array) {
            auto *interaction = dynamic_cast<SpinlessVOperator *>(base);
            if (!interaction) continue;
            interaction->reCalcInv();
            maxInverseResidual = std::max(
                maxInverseResidual,
                (interaction->B * interaction->B_inv - identity).norm());
        }
    } catch (...) {
        evenPbcConstructed = false;
    }
    const bool inversePass = evenPbcConstructed && maxInverseResidual < 1e-12;
    failures += !inversePass;
    out << "even_L_PBC_B_Binv,6,PBC,"
        << (inversePass ? "PASS" : "FAIL") << ',' << maxInverseResidual
        << ",1e-12,two_layer_interaction_operators\n";

    const int L = 4;
    const double V = 2.0;
    const MatType g = deterministicGreen(2 * L);
    for (int boundary = 0; boundary <= 1; ++boundary) {
        SpinlessTvChainUtils interacting(L, 0.1, V, 2, boundary, 0.3, 0.2, 0);
        SpinlessTvChainUtils gaussian(L, 0.1, 0.0, 2, boundary, 0.3, 0.2, 0);
        const DataType helperInteraction =
            interacting.energyFromGreensFunc(g) - gaussian.energyFromGreensFunc(g);
        const DataType directInteraction = directWickInteraction(interacting, g);
        const double error = relativeScalarError(helperInteraction, directInteraction);
        const bool pass = error < 1e-13;
        failures += !pass;
        out << "L4_direct_Wick_interaction_energy,4,"
            << (boundary == 0 ? "PBC" : "OBC") << ','
            << (pass ? "PASS" : "FAIL") << ',' << error
            << ",1e-13," << (boundary == 0 ? 4 : 3) << "_bonds\n";
    }

    std::cout << "{\"status\":\"" << (failures ? "failed" : "complete")
              << "\",\"failures\":" << failures
              << ",\"max_B_Binv_residual\":" << std::setprecision(17)
              << maxInverseResidual << "}\n";
    return failures ? 3 : 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
}
