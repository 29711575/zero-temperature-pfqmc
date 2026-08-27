#ifndef PFQMC_H
#define PFQMC_H

#include "spinless_tV.h"
#include "qr_udt.h"
#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct LeftGreenRecoveryEvent
{
    std::string source;
    int boundary = -1;
    int aux = -1;
    double green_error_before = 0.0;
    double structure_pre_operation = 0.0;
    double structure_delta = 0.0;
    double structure_before = 0.0;
    double structure_after = 0.0;
};

enum class PfQMCSignMode { generic_complex, real_z2 };

class PfQMC
{
public:
    int stb;
    int nDim;
    MatType g;
    std::vector<Operator *> op_array;
    int op_length;
    std::vector<bool> need_stabilization;
    int checkpoints;
    std::vector<UDT> udtL;
    std::vector<UDT> udtR;

    DataType sign;
    PfQMCSignMode sign_mode = PfQMCSignMode::generic_complex;
    int z2_sign = 1;
    long long raw_sign_trusted_count = 0;
    long long raw_sign_untrusted_count = 0;
    long long raw_sign_mismatch_count = 0;
    long long mp_oracle_adjudication_count = 0;
    double max_complex_phase_imag = 0.0;
    bool last_z2_update_used_oracle = false;
    int last_mp_oracle_z2 = 0;
    std::function<int(const std::vector<Operator *> &)> initial_z2_oracle;

    // Off by default so legacy/static callers retain the original update()
    // path bit-for-bit.  The driven driver enables this explicitly.
    bool adaptive_guard = false;
    double guard_threshold = 0.8;
    double guard_ratio_upper = 100.0;
    long long proposal_attempt_count = 0;
    long long adaptive_rebuild_count = 0;
    long long pre_decision_rebuild_count = 0;
    long long post_accept_rebuild_count = 0;
    double min_update_denominator = std::numeric_limits<double>::infinity();

    // Experimental diagnostic fallback.  It is inert unless a caller
    // explicitly installs a rebuild callback and enables it.
    bool multiprecision_fallback = false;
    double multiprecision_core_condition_threshold = std::numeric_limits<double>::infinity();
    std::function<bool(int, MatType &)> multiprecision_rebuild_callback;
    long long multiprecision_fallback_count = 0;
    long long multiprecision_proxy_trigger_count = 0;
    std::vector<double> multiprecision_condition_samples;

    // Optional Green-only recovery for the legacy left sweep. It is inert by
    // default and never changes a proposal, uniform, acceptance, HS field, or
    // transported sign. A triggered operation is rebuilt at most once.
    bool left_green_recovery = false;
    double left_recovery_structure_threshold = 1e-10;
    double left_recovery_structure_delta_threshold = 1e-10;
    long long left_recovery_rank_update_count = 0;
    long long left_recovery_propagation_count = 0;
    std::function<void(const LeftGreenRecoveryEvent &)> left_recovery_event_hook;

    PfQMC(Spinless_tV *walker, int _stb = 10,
          PfQMCSignMode mode = PfQMCSignMode::generic_complex,
          std::function<int(const std::vector<Operator *> &)> z2_oracle = {});

    bool realZ2Mode() const { return sign_mode == PfQMCSignMode::real_z2; }
    int physicalZ2Sign() const { return realZ2Mode() ? z2_sign : (sign.real() >= 0.0 ? 1 : -1); }
    DataType diagnosticComplexPhase() const { return sign; }

    void configureAdaptiveGuard(bool enabled, double threshold,
                                double ratio_upper = 100.0)
    {
        adaptive_guard = enabled;
        guard_threshold = threshold;
        guard_ratio_upper = ratio_upper;
    }

    void configureMultiprecisionFallback(bool enabled, double core_condition_threshold,
                                         std::function<bool(int, MatType &)> callback = {})
    {
        multiprecision_fallback = enabled;
        multiprecision_core_condition_threshold = core_condition_threshold;
        multiprecision_rebuild_callback = std::move(callback);
    }

    void configureLeftSweepGreenRecovery(
        bool enabled, double structure_threshold = 1e-10,
        double structure_delta_threshold = 1e-10)
    {
        left_green_recovery = enabled;
        left_recovery_structure_threshold = structure_threshold;
        left_recovery_structure_delta_threshold = structure_delta_threshold;
    }

    void configureLeftSweepGreenRecoveryEventHook(
        std::function<void(const LeftGreenRecoveryEvent &)> hook)
    {
        left_recovery_event_hook = std::move(hook);
    }

    void rightInit()
    {
        MatType tmp = MatType::Identity(nDim, nDim);
        MatType Aseg = MatType::Identity(nDim, nDim);
        int curSeg = 0;
        for (int l = 0; l < op_length; l++)
        {
            op_array[l]->left_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);
            // %op_length is important, cannot be remove. or else the last segment will not be calculated
            if (need_stabilization[(l + 1) % op_length])
            {
                if (curSeg == 0)
                {
                    udtR[curSeg] = UDT(Aseg); // TODO: performance check
                }
                else
                {
                    udtR[curSeg] = Aseg * udtR[curSeg - 1];
                }
                Aseg = MatType::Identity(nDim, nDim);
                curSeg++;
            }
        }
        udtR[checkpoints - 1].onePlusInv(g);
    }

    void leftInit()
    {
        MatType tmp = MatType::Identity(nDim, nDim);
        MatType Aseg = MatType::Identity(nDim, nDim);
        int curSeg = checkpoints - 1;
        for (int l = op_length - 1; l > -1; l--)
        {
            op_array[l]->right_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);
            if (need_stabilization[l])
            {
                Aseg.adjointInPlace();
                if (curSeg == (checkpoints - 1))
                {
                    udtL[curSeg] = UDT(Aseg); // TODO: performance check
                }
                else
                {
                    udtL[curSeg] = Aseg * udtL[curSeg + 1];
                }
                Aseg = MatType::Identity(nDim, nDim);
                curSeg--;
            }
        }
        udtL[0].onePlusInv(g);
        g.adjointInPlace();
    }

    // after each sweep
    // the greens function g
    // is automatically updated
    // If capture_boundary is in [1, op_length], copy the stabilized/evolved
    // equal-time Green matrix and sign immediately after that many operators.
    // This is used by projector calculations to measure at a fixed midpoint;
    // the default preserves the finite-temperature API and behavior.
    void rightSweep(int capture_boundary = -1, MatType *captured_g = nullptr,
                    DataType *captured_sign = nullptr,
                    int *captured_z2_sign = nullptr,
                    bool *captured_z2_oracle_used = nullptr,
                    int *captured_oracle_z2 = nullptr);
    // Deprecated/debug-only: this cache splice is valid only at the matching
    // checkpoint.  At an arbitrary boundary it omits the unvisited suffix of
    // the current stabilization segment.  Do not use it as a rebuild backend.
    void rebuildRightSweepGreenAtCurrentBoundary(int curSeg, MatType Aseg,
                                                 MatType &out)
    {
        UDT right;
        if (curSeg == 0) right = UDT(Aseg);
        else right = Aseg * udtR[curSeg-1];
        if (curSeg == checkpoints-1) right.onePlusInv(out);
        else out = onePlusInv(udtL[curSeg+1], right);
    }

    // Read-only stabilized rebuild at the boundary immediately BEFORE
    // op_array[boundary_index].  rightSweep's g after operator l is therefore
    // boundary_index=(l+1)%op_length.  Since stabilizedLeftMultiply applies
    // F <- B_i F, iterating l,l+1,... cyclically constructs exactly
    //
    //   B_{l-1} ... B_0 B_{N-1} ... B_l .
    //
    // This is the cyclic product used by rightSweep at that boundary.  Build
    // it from the live operator B matrices in bounded blocks and refactor the
    // accumulated UDT after every block.  No sweep cache, field, proposal, or
    // RNG state is read or modified except for reading op_array's matrices.
    void rebuildGreenFromFullContourAtBoundary(int boundary_index, MatType &out,
                                               int block_size = -1) const
    {
        if (op_length <= 0) {
            out = MatType::Identity(nDim, nDim);
            return;
        }
        int boundary = boundary_index % op_length;
        if (boundary < 0) boundary += op_length;
        // The reference/safety backend refactorizes after every operator.
        // A larger block is an explicit opt-in only: at the frozen flip-254
        // configuration, accumulating stb=10 raw B matrices before the next
        // UDT already loses enough information to defeat the guard.
        const int blockLength = block_size > 0 ? block_size : 1;
        UDT product(nDim);
        MatType block = MatType::Identity(nDim, nDim);
        MatType tmp;
        int inBlock = 0;
        for (int offset = 0; offset < op_length; ++offset) {
            const int index = (boundary + offset) % op_length;
            op_array[index]->left_multiply(block, tmp);
            block.swap(tmp);
            ++inBlock;
            if (inBlock == blockLength || offset + 1 == op_length) {
                product = block * product;
                block.setIdentity();
                inBlock = 0;
            }
        }
        product.onePlusInv(out);
    }

    // Read-only condition proxy for the scaled solve used by UDT::onePlusInv.
    // It is deliberately separate from the normal rebuild path.
    double fullContourCoreConditionAtBoundary(int boundary_index) const
    {
        if (op_length <= 0) return 1.0;
        int boundary = boundary_index % op_length;
        if (boundary < 0) boundary += op_length;
        UDT product(nDim);
        MatType block = MatType::Identity(nDim, nDim), tmp;
        for (int offset = 0; offset < op_length; ++offset) {
            const int index = (boundary + offset) % op_length;
            op_array[index]->left_multiply(block, tmp);
            block.swap(tmp);
            product = block * product;
            block.setIdentity();
        }
        Eigen::FullPivLU<MatType> tFactor(product.T);
        if (!tFactor.isInvertible() || !std::isfinite(tFactor.rcond()) ||
            tFactor.rcond() <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        const MatType identity = MatType::Identity(nDim, nDim);
        MatType lhs = tFactor.solve(identity);
        const double solveResidual =
            (product.T*lhs-identity).norm()/std::max(
                identity.norm(), std::numeric_limits<double>::min());
        if (!lhs.allFinite() || !std::isfinite(solveResidual) ||
            solveResidual > 1e-8) {
            return std::numeric_limits<double>::infinity();
        }
        dVecType dplus(nDim), dminus(nDim);
        for (int i=0;i<nDim;++i) {
#ifdef PFQMC_SCALE_SAFE_UDT
            dplus(i) = product.dLargeInverse(i);
            dminus(i) = product.dSmallPart(i);
#else
            dplus(i) = 1.0 / std::max(product.D(i), 1.0);
            dminus(i) = std::min(product.D(i), 1.0);
#endif
        }
        MatType core = lhs * dplus.asDiagonal() + product.U * dminus.asDiagonal();
        Eigen::JacobiSVD<MatType> svd(core);
        const auto s = svd.singularValues();
        return s(0) / std::max(s(s.size()-1), std::numeric_limits<double>::min());
    }
    DataType leftRecoveryUpdateAtBoundary(Operator *op, int boundary);
    DataType realZ2UpdateAtBoundary(Operator *op);
    bool recoverLeftGreenAfterOperation(const char *source, int boundary,
                                        int aux,
                                        double structure_pre_operation);
    void leftSweep();

    // get sign by computing the pfaffian of
    // a 4N * 4N matrix
    DataType getSignRaw();
    PfaffianResult getSignRawWithStatus();
    PfaffianResult checkRawSignReadOnly();
    double rawContourExponentSpan() const;

private:
    void updatePhysicalZ2(const DataType &phase_factor);

    // should provide same result as getSignRaw
    // but by computing the pfaffian of a 2N * 2N matrix
    // TODO: this method currently has fundamental flaws
    // therefore should not be used
    // DataType getSign();
    
    // ~PfQMC()
    // {
        // for (int i = 0; i < udtR.size(); i++)
        // {
        //     delete udtR[i];
        // }
        // delete Al;
    // }
};

#endif
