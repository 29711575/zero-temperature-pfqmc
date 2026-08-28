#include <cmath>
#include <iomanip>
#include <iostream>

#include "kitaevChain.h"

int main(int argc, char **argv) {
    const int L = argc == 2 ? std::stoi(argv[1]) : 6;
    if (L < 2) {
        std::cerr << "usage: static_contact_test [L>=2]\n";
        return 2;
    }
    SpinlessTvChainUtils config(L, 0.1, 0.0, 40, 1, 1, 0, 0);
    const MatType identity = MatType::Identity(config.nDim, config.nDim);
    const double actual = -config.StructureFactorCDW(identity).real();
    const double expected = 1.0 / (4.0 * L);
    const double difference = std::abs(actual - expected);
    std::cout << std::setprecision(17)
              << "{\"L\":" << L
              << ",\"identity_S_pi\":" << actual
              << ",\"expected_S_pi\":" << expected
              << ",\"abs_diff\":" << difference
              << ",\"status\":\"" << (difference <= 1e-12 ? "PASS" : "FAIL") << "\"}\n";
    return difference <= 1e-12 ? 0 : 1;
}
