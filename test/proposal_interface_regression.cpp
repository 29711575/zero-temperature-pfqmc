#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "honeycomb.h"
#include "kitaevChain.h"
#include "pfqmc.h"

namespace {

class TestSpinlessVOperator : public SpinlessVOperator {
  public:
    using SpinlessVOperator::SpinlessVOperator;
    using SpinlessVOperator::config;
};

struct FlipResult {
    DataType ratio{1.0};
    bool accepted{false};
    DataType signMultiplier{1.0};
};

double matrixMaxError(const MatType &a, const MatType &b) {
    return (a - b).cwiseAbs().maxCoeff();
}

// Test-only transcription of SpinlessVOperator::singleFlip from immediately
// before the proposal API refactor.  In particular, this is the post-local-V
// implementation: all coupling-dependent quantities come from the operator,
// not from its shared config object.
FlipResult legacySingleFlipReference(TestSpinlessVOperator &op, MatType &g,
                                     int idxAux, double rand) {
    FlipResult out;
    DataType r;
    const int auxCur = (*op.s)(idxAux);
    int idx1, idx2, idx3, idx4;
    DataType tmp[2];
    const int inc = 1;
    DataType alpha;

    op.config->aux2MajoranaIdx(idxAux, 0, op.bondType, idx1, idx2);
    op.config->aux2MajoranaIdx(idxAux, 1, op.bondType, idx3, idx4);

    if (op.hsScheme == 0) {
        tmp[0] = 1.0 - DataType(0.0, 1.0) * op.thlV * double(auxCur) * g(idx1, idx2);
        tmp[1] = 1.0 - DataType(0.0, 1.0) * op.thlV * double(auxCur) * g(idx3, idx4);
        r = tmp[0] * tmp[1];
        r += op.thlV * op.thlV *
             (g(idx1, idx3) * g(idx2, idx4) - g(idx2, idx3) * g(idx1, idx4));
        r *= op.etaM;
    } else {
        tmp[0] = 1.0 - DataType(0.0, 1.0) * op.thlV * double(auxCur) * g(idx1, idx3);
        tmp[1] = 1.0 + DataType(0.0, 1.0) * op.thlV * double(auxCur) * g(idx2, idx4);
        r = tmp[0] * tmp[1];
        r -= op.thlV * op.thlV *
             (g(idx1, idx2) * g(idx3, idx4) - g(idx3, idx2) * g(idx1, idx4));
        r *= op.etaM;
    }
    out.ratio = r;
    out.accepted = rand < std::abs(r);
    if (!out.accepted) return out;

    out.signMultiplier = r / std::abs(r);
    if (op.hsScheme == 0) {
        for (int imaj = 0; imaj < 2; ++imaj) {
            op.config->aux2MajoranaIdx(idxAux, imaj, op.bondType, idx1, idx2);
            if (imaj == 1) {
                tmp[1] = 1.0 - DataType(0.0, 1.0) * op.thlV * double(auxCur) * g(idx1, idx2);
            }
            (*op.s)(idxAux) = -auxCur;
            op.B(idx1, idx2) = -op.B(idx1, idx2);
            op.B(idx2, idx1) = -op.B(idx2, idx1);

            cVecType x1 = -g.col(idx1);
            cVecType x2 = -g.col(idx2);
            x1(idx1) += 2;
            x2(idx2) += 2;
            alpha = DataType(0.0, 1.0) * double(auxCur) * op.thlV / tmp[imaj];
            zgeru(&op.nDim, &op.nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                  g.data(), &op.nDim);
            alpha = -alpha;
            zgeru(&op.nDim, &op.nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                  g.data(), &op.nDim);
        }
    } else {
        int idxj1, idxk1, idxj2, idxk2;
        op.config->aux2MajoranaIdx(idxAux, 0, op.bondType, idxj1, idxk1);
        op.config->aux2MajoranaIdx(idxAux, 1, op.bondType, idxj2, idxk2);
        for (int iaux = 0; iaux < 2; ++iaux) {
            if (iaux == 0) {
                idx1 = idxj1;
                idx2 = idxj2;
            } else {
                idx1 = idxk1;
                idx2 = idxk2;
                tmp[1] = 1.0 + DataType(0.0, 1.0) * op.thlV * double(auxCur) * g(idx1, idx2);
            }
            (*op.s)(idxAux) = -auxCur;
            op.B(idx1, idx2) = -op.B(idx1, idx2);
            op.B(idx2, idx1) = -op.B(idx2, idx1);

            cVecType x1 = -g.col(idx1);
            cVecType x2 = -g.col(idx2);
            x1(idx1) += 2;
            x2(idx2) += 2;
            alpha = DataType(0.0, 1.0) * double(auxCur) * op.thlV / tmp[iaux];
            if (iaux == 1) alpha = -alpha;
            zgeru(&op.nDim, &op.nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                  g.data(), &op.nDim);
            alpha = -alpha;
            zgeru(&op.nDim, &op.nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                  g.data(), &op.nDim);
        }
    }
    return out;
}

DataType legacyOperatorUpdate(TestSpinlessVOperator &op, MatType &g) {
    DataType sign = 1.0;
    for (int aux = 0; aux < op.s->size(); ++aux) {
        const double u = op.rd->rdUniform01();
        const FlipResult result = legacySingleFlipReference(op, g, aux, u);
        sign *= result.signMultiplier;
    }
    return sign;
}

iVecType initialHs(int n) {
    iVecType hs(n);
    for (int i = 0; i < n; ++i) hs(i) = ((i * 7 + 3) % 5 < 2) ? -1 : 1;
    return hs;
}

MatType initialGreen(const TestSpinlessVOperator &op) {
    const MatType identity = MatType::Identity(op.nDim, op.nDim);
    return 2.0 * (identity + op.B).inverse() - identity;
}

bool exactNextRng(rdGenerator rngA, rdGenerator rngB, int count = 8) {
    for (int i = 0; i < count; ++i) {
        if (rngA.rdUniform01() != rngB.rdUniform01()) return false;
    }
    return true;
}

bool pendingRecomputeRegression(const SpinlessTvHoneycombUtils &config,
                                const iVecType &hs, double localV) {
    rdGenerator rng(481516);
    TestSpinlessVOperator op(&config, new iVecType(hs), 1, &rng, localV);
    MatType g = initialGreen(op);
    double uniform = -1.0;
    if (!op.prepareSingleFlip(g, &uniform)) return false;

    const int aux = op.pendingAux;
    const int oldSigma = op.pendingOldSigma;
    const int newSigma = op.pendingNewSigma;
    const DataType ratio = op.pendingRatio;
    const DataType denom1 = op.pendingDenom1;
    const DataType denom2 = op.pendingDenom2;
    const double savedUniform = op.pendingUniform;
    const iVecType hsSaved = *op.s;
    const MatType bSaved = op.B;
    const MatType gSaved = g;
    const rdGenerator rngSaved = rng;

    for (int i = 0; i < 3; ++i) {
        const DataType recomputed = op.recomputePreparedRatio(g);
        if (recomputed != ratio || op.pendingRatio != ratio ||
            op.pendingDenom1 != denom1 || op.pendingDenom2 != denom2 ||
            op.pendingAux != aux || op.pendingOldSigma != oldSigma ||
            op.pendingNewSigma != newSigma || op.pendingUniform != savedUniform ||
            uniform != savedUniform || !op.pendingValid ||
            !((*op.s).array() == hsSaved.array()).all() ||
            matrixMaxError(op.B, bSaved) != 0.0 || matrixMaxError(g, gSaved) != 0.0 ||
            !exactNextRng(rng, rngSaved)) return false;
    }
    return true;
}

struct DrivenWalker : Spinless_tV {
    std::vector<std::string> region;

    DrivenWalker(const SpinlessTvChainUtils *config, rdGenerator *rng,
                 double v0, double vf, double rate, double theta, double beta) {
        nDim = config->nDim;
        const int ninit = std::llround(theta / config->dt);
        const int nramp = std::llround((vf - v0) / (rate * config->dt));
        const int nb[2] = {(config->Lx + 1) / 2, config->Lx / 2};
        MatType h(nDim, nDim);
        config->KineticGenerator(h);
        for (double remaining = beta; remaining > 1e-12;) {
            const double step = std::min(config->dt, remaining);
            MatType x = h;
            const MatType b = expm(x, -step);
            x = step * h;
            op_array.push_back(new DenseOperator(b, signOfHamiltonian(x)));
            region.push_back("trial");
            remaining -= step;
        }
        MatType x = h;
        const MatType kh = expm(x, -config->dt / 2.0);
        x = config->dt * h / 2.0;
        const DataType sk = signOfHamiltonian(x);
        auto addSlice = [&](double v, bool dagger, const std::string &label) {
            op_array.push_back(new DenseOperator(kh, sk));
            region.push_back(label);
            for (int z = 0; z < 2; ++z) {
                const int bond = dagger ? 1 - z : z;
                auto *s = new iVecType(nb[bond]);
                for (int j = 0; j < s->size(); ++j) (*s)(j) = rng->rdZ2();
                op_array.push_back(new SpinlessVOperator(config, s, bond, rng, v));
                region.push_back(label);
            }
            op_array.push_back(new DenseOperator(kh, sk));
            region.push_back(label);
        };
        for (int l = 0; l < ninit; ++l) addSlice(v0, false, "ket_init");
        for (int l = 0; l < nramp; ++l) {
            const char *label = l < nramp / 3 ? "ket_early" :
                                l < 2 * nramp / 3 ? "ket_middle" : "ket_late";
            addSlice(v0 + rate * (l + 0.5) * config->dt, false, label);
        }
        for (int l = nramp - 1; l >= 0; --l) {
            const char *label = l >= 2 * nramp / 3 ? "bra_late" :
                                l >= nramp / 3 ? "bra_middle" : "bra_early";
            addSlice(v0 + rate * (l + 0.5) * config->dt, true, label);
        }
        for (int l = 0; l < ninit; ++l) addSlice(v0, true, "bra_init");
    }
};

MatType stabilizedGreenAt(DrivenWalker &walker, int boundary) {
    UDT product(walker.nDim);
    const int n = walker.op_array.size();
    for (int j = 0; j < n; ++j)
        walker.op_array[(boundary + j) % n]->stabilizedLeftMultiply(product);
    MatType result;
    product.onePlusInv(result);
    return result;
}

struct SweepBoundaryState {
    int curSeg = 0;
    MatType aseg;
    explicit SweepBoundaryState(int n) : aseg(MatType::Identity(n, n)) {}
};

void advanceRightBoundary(PfQMC &q, SweepBoundaryState &state, int opIndex) {
    MatType tmp;
    q.op_array[opIndex]->left_multiply(state.aseg, tmp);
    state.aseg.swap(tmp);
    if (q.need_stabilization[(opIndex + 1) % q.op_length]) {
        if (state.curSeg == 0) q.udtR[state.curSeg] = UDT(state.aseg);
        else q.udtR[state.curSeg] = state.aseg * q.udtR[state.curSeg - 1];
        state.aseg.setIdentity();
        ++state.curSeg;
    }
}

bool proposalIdentityPreserved(const SpinlessVOperator &op, int aux, int oldSigma,
                               int newSigma, double uniform, const iVecType &hs,
                               const MatType &b, const rdGenerator &rngSaved) {
    return op.pendingValid && op.pendingAux == aux &&
           op.pendingOldSigma == oldSigma && op.pendingNewSigma == newSigma &&
           op.pendingUniform == uniform &&
           ((*op.s).array() == hs.array()).all() && matrixMaxError(op.B, b) == 0.0 &&
           exactNextRng(*op.rd, rngSaved);
}

struct NormalBoundaryResult {
    double maxGreenError = 0.0, maxCheckpointError = 0.0;
    double maxNonCheckpointError = 0.0;
    int early = 0, middle = 0, late = 0;
    bool proposalPreserved = true;
    bool helperPassed = true;
    int checkpointsTested = 0;
};

NormalBoundaryResult normalBoundaryRegression() {
    SpinlessTvChainUtils config(6, 0.05, 0.0, 1, 1, 1.0, 0.0, 0);
    rdGenerator rng(424242);
    DrivenWalker walker(&config, &rng, 0.0, 4.0, 1.0, 6.0, 8.0);
    PfQMC q(&walker, 10);
    SweepBoundaryState state(q.nDim);
    NormalBoundaryResult result;
    bool checkpointSelected = false;
    for (int oi = 0; oi < q.op_length; ++oi) {
        const std::string &label = walker.region[oi];
        int *count = label.find("early") != std::string::npos ? &result.early :
                     label.find("middle") != std::string::npos ? &result.middle :
                     label.find("late") != std::string::npos ? &result.late : nullptr;
        const bool nonCheckpointSample = count && *count < 1 &&
            !q.need_stabilization[oi] &&
            dynamic_cast<SpinlessVOperator *>(q.op_array[oi]);
        const bool checkpointSample = q.need_stabilization[oi] && oi > 0 &&
                                      !checkpointSelected;
        MatType rebuiltAll;
        double allError = 0.0;
        if (nonCheckpointSample || checkpointSample) {
            q.rebuildGreenFromFullContourAtBoundary(oi, rebuiltAll);
            const MatType fullReference = stabilizedGreenAt(walker, oi);
            allError = matrixMaxError(rebuiltAll, fullReference);
            result.maxGreenError = std::max(result.maxGreenError, allError);
        }
        if (checkpointSample) {
            const int completedSegment = state.curSeg - 1;
            MatType checkpointGreen;
            if (completedSegment == q.checkpoints - 1)
                q.udtR[completedSegment].onePlusInv(checkpointGreen);
            else
                checkpointGreen = onePlusInv(q.udtL[completedSegment + 1],
                                             q.udtR[completedSegment]);
            result.maxCheckpointError = std::max(
                result.maxCheckpointError,
                std::max(allError, matrixMaxError(rebuiltAll, checkpointGreen)));
            checkpointSelected = true;
            ++result.checkpointsTested;
        } else if (nonCheckpointSample) {
            result.maxNonCheckpointError = std::max(result.maxNonCheckpointError,
                                                     allError);
        }
        if (nonCheckpointSample) {
            auto *op = dynamic_cast<SpinlessVOperator *>(q.op_array[oi]);
            double uniform = -1.0;
            op->prepareSingleFlip(q.g, &uniform);
            const int aux = op->pendingAux;
            const int oldSigma = op->pendingOldSigma;
            const int newSigma = op->pendingNewSigma;
            const double savedUniform = op->pendingUniform;
            const iVecType hs = *op->s;
            const MatType b = op->B;
            const rdGenerator rngSaved = rng;
            MatType rebuilt;
            q.rebuildGreenFromFullContourAtBoundary(oi, rebuilt);
            op->recomputePreparedRatio(rebuilt);
            result.proposalPreserved = result.proposalPreserved &&
                proposalIdentityPreserved(*op, aux, oldSigma, newSigma,
                                          savedUniform, hs, b, rngSaved) &&
                uniform == savedUniform;
            ++*count;
        }
        advanceRightBoundary(q, state, oi);
    }
    result.helperPassed = result.checkpointsTested == 1 &&
                          result.early >= 1 && result.middle >= 1 && result.late >= 1 &&
                          result.proposalPreserved &&
                          result.maxGreenError <= 1e-10;
    return result;
}

struct FrozenResult {
    double ratioBefore = 0.0, ratioAfter = 0.0;
    double ratioBeforeError = 0.0, ratioAfterError = 0.0;
    double greenBeforeError = 0.0, greenAfterError = 0.0;
    bool proposalPreserved = false;
    bool found = false;
};

FrozenResult frozenRegression(int targetFlip, int stb, double exactRatio) {
    SpinlessTvChainUtils config(6, 0.1, 0.0, 1, 1, 1.0, 0.0, 0);
    rdGenerator rng(424242);
    DrivenWalker walker(&config, &rng, 0.0, 5.0, 1.0, 6.0, 8.0);
    PfQMC q(&walker, stb);
    SweepBoundaryState state(q.nDim);
    MatType g = stabilizedGreenAt(walker, 0), tmp;
    FrozenResult result;
    int flip = 0;
    int opsSinceReset = 0;
    for (int oi = 0; oi < q.op_length && !result.found; ++oi) {
        if (auto *op = dynamic_cast<SpinlessVOperator *>(q.op_array[oi])) {
            const int aux = flip % op->s->size();
            op->pendingCursor = aux;
            double uniform = -1.0;
            if (!op->prepareSingleFlip(g, &uniform)) break;
            if (flip == targetFlip) {
                const int savedAux = op->pendingAux;
                const int oldSigma = op->pendingOldSigma;
                const int newSigma = op->pendingNewSigma;
                const double savedUniform = op->pendingUniform;
                const iVecType hs = *op->s;
                const MatType b = op->B;
                const rdGenerator rngSaved = rng;
                const MatType exactGreen = stabilizedGreenAt(walker, oi);
                result.ratioBefore = std::abs(op->pendingRatio);
                result.greenBeforeError = matrixMaxError(g, exactGreen);
                MatType rebuilt;
                q.rebuildGreenFromFullContourAtBoundary(oi, rebuilt);
                op->recomputePreparedRatio(rebuilt);
                result.ratioAfter = std::abs(op->pendingRatio);
                result.greenAfterError = matrixMaxError(rebuilt, exactGreen);
                result.ratioBeforeError = std::abs(result.ratioBefore - exactRatio);
                result.ratioAfterError = std::abs(result.ratioAfter - exactRatio);
                result.proposalPreserved = proposalIdentityPreserved(
                    *op, savedAux, oldSigma, newSigma, savedUniform, hs, b, rngSaved) &&
                    savedAux == aux && uniform == savedUniform;
                result.found = true;
            }
            op->finishSingleFlip(g, true, true);
            ++flip;
        }
        q.op_array[oi]->left_propagate(g, tmp);
        advanceRightBoundary(q, state, oi);
        ++opsSinceReset;
        if (opsSinceReset >= stb && oi + 1 < q.op_length) {
            g = stabilizedGreenAt(walker, oi + 1);
            opsSinceReset = 0;
        }
    }
    return result;
}

}  // namespace

int main() {
    mkl_set_num_threads(1);
    constexpr int kFields = 100;
    constexpr int kOperatorUpdates = 12;
    constexpr double kTol = 2e-12;
    constexpr double kLocalV = 1.37;
    const int seed = 20260814;

    SpinlessTvHoneycombUtils config(4, 4, 0.08, 0.61, 4);
    const iVecType hs0 = initialHs(config.nUnitcell);
    const bool pendingPass = pendingRecomputeRegression(config, hs0, kLocalV);
    bool ok = pendingPass;

    rdGenerator rngA(seed), rngB(seed);
    TestSpinlessVOperator opA(&config, new iVecType(hs0), 2, &rngA, kLocalV);
    TestSpinlessVOperator opB(&config, new iVecType(hs0), 2, &rngB, kLocalV);
    MatType gA = initialGreen(opA);
    MatType gB = gA;
    DataType signA = 1.0, signB = 1.0;
    double maxRatioError = 0.0, maxBError = 0.0, maxGreenError = 0.0;
    double maxSignError = 0.0;
    bool uniformExact = true, acceptExact = true, hsExact = true;
    bool proposalSequenceExact = true;

    for (int step = 0; step < kFields; ++step) {
        const int aux = step % opA.s->size();
        const double uA = rngA.rdUniform01();
        const FlipResult a = legacySingleFlipReference(opA, gA, aux, uA);

        double uB = -1.0;
        ok = opB.prepareSingleFlip(gB, &uB) && ok;
        proposalSequenceExact = proposalSequenceExact && opB.pendingAux == aux;
        const DataType ratioB = opB.pendingRatio;
        const bool acceptB = uB < std::abs(ratioB);
        const DataType signMultiplierB = opB.finishSingleFlip(gB, acceptB, true);
        signA *= a.signMultiplier;
        signB *= signMultiplierB;

        uniformExact = uniformExact && (uA == uB);
        acceptExact = acceptExact && (a.accepted == acceptB);
        hsExact = hsExact && ((*opA.s).array() == (*opB.s).array()).all();
        maxRatioError = std::max(maxRatioError, std::abs(a.ratio - ratioB));
        maxBError = std::max(maxBError, matrixMaxError(opA.B, opB.B));
        maxGreenError = std::max(maxGreenError, matrixMaxError(gA, gB));
        maxSignError = std::max(maxSignError, std::abs(a.signMultiplier - signMultiplierB));
        maxSignError = std::max(maxSignError, std::abs(signA - signB));
    }
    const bool fieldRngExact = exactNextRng(rngA, rngB);
    const bool fieldPass = uniformExact && proposalSequenceExact && acceptExact && hsExact &&
                           fieldRngExact && maxRatioError <= kTol && maxBError <= kTol &&
                           maxGreenError <= kTol && maxSignError <= kTol;
    ok = ok && fieldPass;

    rdGenerator wholeRngA(seed), wholeRngB(seed);
    TestSpinlessVOperator wholeA(&config, new iVecType(hs0), 2, &wholeRngA, kLocalV);
    TestSpinlessVOperator wholeB(&config, new iVecType(hs0), 2, &wholeRngB, kLocalV);
    MatType wholeGA = initialGreen(wholeA);
    MatType wholeGB = wholeGA;
    DataType wholeSignA = 1.0, wholeSignB = 1.0;
    double wholeBError = 0.0, wholeGreenError = 0.0, wholeSignError = 0.0;
    for (int i = 0; i < kOperatorUpdates; ++i) {
        wholeSignA *= legacyOperatorUpdate(wholeA, wholeGA);
        wholeSignB *= wholeB.update(wholeGB);
        wholeBError = std::max(wholeBError, matrixMaxError(wholeA.B, wholeB.B));
        wholeGreenError = std::max(wholeGreenError, matrixMaxError(wholeGA, wholeGB));
        wholeSignError = std::max(wholeSignError, std::abs(wholeSignA - wholeSignB));
    }
    const bool wholeHsExact = ((*wholeA.s).array() == (*wholeB.s).array()).all();
    const bool wholeRngExact = exactNextRng(wholeRngA, wholeRngB);
    const bool wholePass = wholeHsExact && wholeRngExact && wholeBError <= kTol &&
                           wholeGreenError <= kTol && wholeSignError <= kTol;
    ok = ok && wholePass;

    const NormalBoundaryResult normal = normalBoundaryRegression();
    const FrozenResult flip167 = frozenRegression(167, 10, 0.07652623318902);
    const FrozenResult flip254 = frozenRegression(254, 1, 0.001180054933443);
    const bool frozenIdentity = normal.proposalPreserved &&
                                flip167.proposalPreserved && flip254.proposalPreserved;
    const bool boundaryHelperPass = normal.helperPassed && flip167.found && flip254.found;
    ok = ok && boundaryHelperPass && frozenIdentity;

    std::cout << std::setprecision(17)
              << "pending_recompute_x3=" << (pendingPass ? "PASS" : "FAIL") << '\n'
              << "field_count=" << kFields << '\n'
              << "independent_same_seed_rng=true\n"
              << "proposal_sequence_exact=" << proposalSequenceExact << '\n'
              << "uniform_sequence_exact=" << uniformExact << '\n'
              << "max_ratio_error=" << maxRatioError << '\n'
              << "max_B_error=" << maxBError << '\n'
              << "max_Green_error=" << maxGreenError << '\n'
              << "max_sign_error=" << maxSignError << '\n'
              << "accept_sequence_exact=" << acceptExact << '\n'
              << "field_rng_exact=" << fieldRngExact << '\n'
              << "whole_max_B_error=" << wholeBError << '\n'
              << "whole_max_Green_error=" << wholeGreenError << '\n'
              << "whole_max_sign_error=" << wholeSignError << '\n'
              << "whole_HS_exact=" << wholeHsExact << '\n'
              << "whole_rng_exact=" << wholeRngExact << '\n'
              << "whole_operator=" << (wholePass ? "PASS" : "FAIL") << '\n'
              << "boundary_normal_max_Green_error=" << normal.maxGreenError << '\n'
              << "boundary_checkpoint_max_Green_error=" << normal.maxCheckpointError << '\n'
              << "boundary_noncheckpoint_max_Green_error=" << normal.maxNonCheckpointError << '\n'
              << "boundary_normal_counts=" << normal.early << ',' << normal.middle << ',' << normal.late << '\n'
              << "boundary_checkpoints_tested=" << normal.checkpointsTested << '\n'
              << "boundary_proposal_uniform_rng_exact=" << frozenIdentity << '\n'
              << "flip167_R_before=" << flip167.ratioBefore << '\n'
              << "flip167_R_before_error=" << flip167.ratioBeforeError << '\n'
              << "flip167_R_after=" << flip167.ratioAfter << '\n'
              << "flip167_R_after_error=" << flip167.ratioAfterError << '\n'
              << "flip167_Green_before_error=" << flip167.greenBeforeError << '\n'
              << "flip167_Green_after_error=" << flip167.greenAfterError << '\n'
              << "flip167_identity_exact=" << flip167.proposalPreserved << '\n'
              << "flip254_R_before=" << flip254.ratioBefore << '\n'
              << "flip254_R_before_error=" << flip254.ratioBeforeError << '\n'
              << "flip254_R_after=" << flip254.ratioAfter << '\n'
              << "flip254_R_after_error=" << flip254.ratioAfterError << '\n'
              << "flip254_Green_before_error=" << flip254.greenBeforeError << '\n'
              << "flip254_Green_after_error=" << flip254.greenAfterError << '\n'
              << "flip254_identity_exact=" << flip254.proposalPreserved << '\n'
              << "boundary_helper=" << (boundaryHelperPass ? "PASS" : "FAIL") << '\n'
              << "overall=" << (ok ? "PASS" : "FAIL") << '\n';
    return ok ? 0 : 1;
}
