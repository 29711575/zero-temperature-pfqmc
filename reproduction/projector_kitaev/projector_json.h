#ifndef PFQMC_PROJECTOR_JSON_H
#define PFQMC_PROJECTOR_JSON_H

#include <cmath>
#include <ostream>
#include <string>
#include "pfqmc.h"

struct ProjectorRawSignChecks
{
    long long comparisons = 0;
    long long agreement = 0;
    long long mismatch = 0;
    long long unavailable = 0;
    long long untrusted = 0;
    int last_lapack_info = 0;
    double min_pivot = std::numeric_limits<double>::infinity();

    void record(const DataType &transported, const PfaffianResult &raw,
                double phase_tolerance = 1e-2)
    {
        ++comparisons;
        min_pivot = std::min(min_pivot, raw.min_pivot);
        if (!raw.ok()) {
            if (raw.untrusted()) ++untrusted;
            else ++unavailable;
            last_lapack_info = raw.lapack_info;
            return;
        }
        const double magnitude = std::abs(raw.value);
        if (!std::isfinite(raw.value.real()) || !std::isfinite(raw.value.imag()) ||
            !std::isfinite(magnitude) || magnitude == 0.0 ||
            std::abs(magnitude-1.0) > phase_tolerance) {
            ++untrusted;
            return;
        }
        if ((transported.real() >= 0.0) == (raw.value.real() >= 0.0)) ++agreement;
        else ++mismatch;
    }

    const char *status() const
    {
        if (comparisons == 0) return "not_sampled";
        if (unavailable != 0) return "raw_check_unavailable";
        if (untrusted != 0) return "raw_check_untrusted";
        if (mismatch != 0) return "mismatch";
        return "agreement";
    }
};

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
        << ",\"condition_aware_ratio_enabled\":false"
        << ",\"sign_mode\":\""
        << (qmc.realZ2Mode() ? "real_z2" : "generic_complex") << "\""
        << ",\"raw_sign_trusted_count\":" << qmc.raw_sign_trusted_count
        << ",\"raw_sign_untrusted_count\":" << qmc.raw_sign_untrusted_count
        << ",\"raw_sign_mismatch_count\":" << qmc.raw_sign_mismatch_count
        << ",\"mp_oracle_adjudication_count\":" << qmc.mp_oracle_adjudication_count
        << ",\"mp_checkpoint_mutating\":false"
        << ",\"mp_trusted_count\":" << qmc.mp_trusted_count
        << ",\"mp_untrusted_count\":" << qmc.mp_untrusted_count
        << ",\"mp_precision_escalation_count\":" << qmc.mp_precision_escalation_count
        << ",\"mp_max_precision_digits\":" << qmc.mp_max_precision_digits
        << ",\"mp_candidate_mismatch_count\":" << qmc.mp_candidate_mismatch_count
        << ",\"mp_correction_count\":" << qmc.mp_correction_count
        << ",\"max_z2_ratio_reality_error\":" << qmc.max_z2_ratio_reality_error
        << ",\"mp_canonical_order\":true"
        << ",\"real_z2_policy\":\"transported_ratio_z2_mp_record_only\"";
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

inline void projectorJsonRawSignChecks(std::ostream &out,
                                       const ProjectorRawSignChecks &checks)
{
    out << ",\"raw_sign_check_status\":\"" << checks.status() << "\""
        << ",\"raw_sign_check_comparisons\":" << checks.comparisons
        << ",\"raw_sign_check_agreement\":" << checks.agreement
        << ",\"raw_sign_check_mismatch\":" << checks.mismatch
        << ",\"raw_sign_check_unavailable\":" << checks.unavailable
        << ",\"raw_sign_check_untrusted\":" << checks.untrusted
        << ",\"raw_sign_check_last_lapack_info\":"
        << checks.last_lapack_info
        << ",\"raw_sign_check_min_pivot\":";
    projectorJsonNumber(out, checks.min_pivot,
                        std::isfinite(checks.min_pivot));
}

#endif
