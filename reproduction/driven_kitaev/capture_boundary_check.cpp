#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../projector_kitaev/projector_contour.h"
#include "pfqmc.h"

struct CaptureWalker : Spinless_tV {
    int center = -1;
    int trialSlices = 0;
    int physicalSlices = 0;

    CaptureWalker(const SpinlessTvChainUtils *config, rdGenerator *random) {
        build_projector_static_contour(
            *this, config, random, 2.0, 4.0, center, trialSlices, physicalSlices);
    }
};

int main() {
    SpinlessTvChainUtils config(4, 0.1, 2.0, 40, 1, 1, 0, 0);
    rdGenerator random(771);
    CaptureWalker walker(&config, &random);
    PfQMC qmc(&walker, 10);
    MatType captured;
    DataType capturedSign;
    bool rejectedOutOfRange = false;
    bool rejectedPointerWithoutBoundary = false;
    bool validCapture = false;
    try {
        qmc.rightSweep(0, &captured, &capturedSign);
    } catch (const std::out_of_range &) {
        rejectedOutOfRange = true;
    }
    try {
        qmc.rightSweep(-1, &captured, &capturedSign);
    } catch (const std::invalid_argument &) {
        rejectedPointerWithoutBoundary = true;
    }
    try {
        qmc.rightSweep(walker.center, &captured, &capturedSign);
        validCapture = captured.rows() == config.nDim && captured.cols() == config.nDim;
    } catch (...) {
        validCapture = false;
    }
    const bool pass = rejectedOutOfRange && rejectedPointerWithoutBoundary && validCapture;
    std::cout << std::setprecision(17)
              << "{\"rejected_out_of_range\":" << (rejectedOutOfRange ? "true" : "false")
              << ",\"rejected_pointer_without_boundary\":"
              << (rejectedPointerWithoutBoundary ? "true" : "false")
              << ",\"valid_capture\":" << (validCapture ? "true" : "false")
              << ",\"status\":\"" << (pass ? "PASS" : "FAIL") << "\"}\n";
    return pass ? 0 : 1;
}
