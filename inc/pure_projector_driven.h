#ifndef PURE_PROJECTOR_DRIVEN_H
#define PURE_PROJECTOR_DRIVEN_H

#include "pure_projector_fast.h"

#include <complex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct PureDrivenContourResult {
    PureFastConfiguration configuration;
    int ket_factor_count = 0;
    int bra_factor_count = 0;
    bool bra_ket_hs_independent = true;
    std::string contour_order =
        "ket_canonical_action_then_bra_strict_reverse_adjoint";
};

inline PureDrivenContourResult pureProjectorIndependentBraKetContour(
        PureFastConfiguration ketCanonicalAction,
        PureFastConfiguration braCanonicalAction) {
    if (ketCanonicalAction.slices.empty() || braCanonicalAction.slices.empty() ||
        ketCanonicalAction.slices.size() != ketCanonicalAction.hs_fields.size() ||
        ketCanonicalAction.slices.size() != ketCanonicalAction.locations.size() ||
        braCanonicalAction.slices.size() != braCanonicalAction.hs_fields.size() ||
        braCanonicalAction.slices.size() != braCanonicalAction.locations.size() ||
        ketCanonicalAction.slices.size() != braCanonicalAction.slices.size())
        throw std::invalid_argument("invalid independent bra/ket branch factors");

    PureDrivenContourResult result;
    result.ket_factor_count = int(ketCanonicalAction.slices.size());
    result.bra_factor_count = int(braCanonicalAction.slices.size());
    result.configuration = std::move(ketCanonicalAction);
    for (int index = int(braCanonicalAction.slices.size()) - 1;
         index >= 0; --index) {
        const PureProjectorSlice &source = braCanonicalAction.slices[index];
        PureSliceLocation location = braCanonicalAction.locations[index];
        location.branch = PureBranch::Bra;
        result.configuration.slices.emplace_back(
            source.matrix.adjoint().eval(), std::conj(source.eta),
            "bra_adjoint:" + source.label);
        result.configuration.hs_fields.push_back(
            braCanonicalAction.hs_fields[index]);
        result.configuration.locations.push_back(location);
    }
    return result;
}

#endif
