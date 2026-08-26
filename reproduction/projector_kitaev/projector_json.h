#ifndef PFQMC_PROJECTOR_JSON_H
#define PFQMC_PROJECTOR_JSON_H

#include <cmath>
#include <ostream>
#include "pfqmc.h"

inline void projectorJsonNumber(std::ostream &out, double value,
                                bool resolved = true)
{
    if (resolved && std::isfinite(value)) out << value;
    else out << "null";
}

inline void projectorJsonBuildProvenance(std::ostream &out, const PfQMC &qmc)
{
#ifdef PFQMC_SCALE_SAFE_UDT
    out << ",\"scale_safe_udt\":true";
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
    out << ",\"udt_rank_loss_guard\":true"
        << ",\"udt_rank_loss_guard_bits\":"
        << scaleSafeUDTRankLossGuardBits
        << ",\"udt_orthogonality_gate\":"
        << scaleSafeUDTPartialQOrthogonalityLimit;
#else
    out << ",\"udt_rank_loss_guard\":false"
        << ",\"udt_rank_loss_guard_bits\":null"
        << ",\"udt_orthogonality_gate\":null";
#endif
#else
    out << ",\"scale_safe_udt\":false"
        << ",\"udt_rank_loss_guard\":false"
        << ",\"udt_rank_loss_guard_bits\":null"
        << ",\"udt_orthogonality_gate\":null";
#endif
    out << ",\"left_recovery_enabled\":"
        << (qmc.left_green_recovery ? "true" : "false")
        << ",\"condition_aware_ratio_enabled\":false";
#ifdef PFQMC_SOURCE_COMMIT
    out << ",\"source_commit\":\"" << PFQMC_SOURCE_COMMIT << "\"";
#else
    out << ",\"source_commit\":null";
#endif
#ifdef PFQMC_EXECUTABLE_SHA256
    out << ",\"executable_sha256\":\"" << PFQMC_EXECUTABLE_SHA256 << "\"";
#else
    out << ",\"executable_sha256\":null";
#endif
}

#endif
