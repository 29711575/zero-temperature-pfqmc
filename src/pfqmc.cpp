#include "pfqmc.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {
double leftRecoveryRelativeError(const MatType &a, const MatType &b)
{
    return (a-b).norm()/std::max(b.norm(), std::numeric_limits<double>::min());
}

double leftRecoveryStructureResidual(const MatType &green)
{
    const MatType residual = green + green.transpose()
                           - 2.0*MatType::Identity(green.rows(), green.cols());
    return residual.norm()/std::max(green.norm(),
                                    std::numeric_limits<double>::min());
}

bool leftRecoveryFinite(const MatType &green)
{
    for (int j=0; j<green.cols(); ++j)
        for (int i=0; i<green.rows(); ++i)
            if (!std::isfinite(green(i,j).real()) ||
                !std::isfinite(green(i,j).imag())) return false;
    return true;
}
}

PfQMC::PfQMC(Spinless_tV *walker, int _stb, PfQMCSignMode mode,
             std::function<MpZ2Result(const std::vector<Operator *> &)> z2_oracle)
{
    if (walker == nullptr) {
        throw std::invalid_argument("PfQMC requires a non-null walker");
    }
    if (_stb <= 0) {
        throw std::invalid_argument("PfQMC stabilization interval must be positive");
    }
    if (walker->nDim <= 0 || walker->op_array.empty()) {
        throw std::invalid_argument("PfQMC requires a non-empty contour");
    }
    stb = _stb;
    sign_mode = mode;
    initial_z2_oracle = std::move(z2_oracle);
    nDim = walker->nDim;
    g = MatType::Identity(nDim, nDim);
    op_array = walker->op_array;
    op_length = op_array.size();
    need_stabilization = std::vector<bool>(op_length);
    checkpoints = 0;
    // if a checkpoint is reached
    // stabilized Green's function is re-evaluated
    for (int l = 0; l < op_length; l++)
    {
        bool flag = ((l % stb) == 0);
        need_stabilization[l] = flag;
        if (flag)
        {
            checkpoints++;
        }
    }
    udtL = std::vector<UDT>(checkpoints);
    udtR = std::vector<UDT>(checkpoints);
    leftInit();
    rightInit();
    const PfaffianResult initialSign = getSignRawWithStatus();
    if (!initialSign.ok() && !initialSign.untrusted()) {
        throw std::runtime_error(
            std::string("initial raw sign unavailable: ") +
            pfaffianStatusName(initialSign.status));
    }
    sign = initialSign.value;
    max_complex_phase_imag = std::abs(sign.imag());
    if (realZ2Mode()) {
        if (initialSign.ok()) {
            ++raw_sign_trusted_count;
            z2_sign = initialSign.value.real() >= 0.0 ? 1 : -1;
        } else {
            ++raw_sign_untrusted_count;
            if (!initial_z2_oracle)
                throw std::runtime_error("real-Z2 initialization requires an oracle when raw sign is untrusted");
            const MpZ2Result oracleResult = initial_z2_oracle(op_array);
            ++mp_oracle_adjudication_count;
            recordMpZ2Result(oracleResult);
            if (!oracleResult.trusted())
                throw std::runtime_error(
                    std::string("real-Z2 initialization oracle is not trusted: ") +
                    mpZ2StatusName(oracleResult.status) + ": " + oracleResult.message);
            z2_sign = oracleResult.z2;
        }
    }
}

void PfQMC::updatePhysicalZ2(const DataType &phaseFactor)
{
    if (!realZ2Mode()) return;
    last_z2_update_used_oracle=false;
    last_mp_oracle_z2=0;
    const double magnitude = std::abs(phaseFactor);
    if (!std::isfinite(phaseFactor.real()) || !std::isfinite(phaseFactor.imag()) ||
        !std::isfinite(magnitude) || magnitude == 0.0)
        throw std::runtime_error("real-Z2 update received a nonfinite/zero phase factor");
    const double reality = std::abs(phaseFactor.imag()) / std::max(std::abs(phaseFactor.real()), 1e-300);
    if (std::isfinite(reality))
        max_z2_ratio_reality_error = std::max(max_z2_ratio_reality_error, reality);
    // In a real-weight policy the imaginary component is a conditioning
    // diagnostic, not a physical phase.  Require a clear real-axis sign
    // margin; machine-level phase accumulation must not drive Z2 or invoke a
    // mutating oracle.  Ratios outside this margin still fail closed.
    if (!std::isfinite(reality) || reality > 0.25 || std::abs(phaseFactor.real()) < 1e-12) {
        if (!initial_z2_oracle) {
            std::ostringstream message;
            message << std::setprecision(17)
                    << "real-Z2 update received a significantly complex/indeterminate ratio: real="
                    << phaseFactor.real() << "; imag=" << phaseFactor.imag()
                    << "; imag_over_real=" << reality
                    << "; no trusted ratio-adjudication oracle is installed";
            throw std::runtime_error(message.str());
        }
        const MpZ2Result oracleResult = initial_z2_oracle(op_array);
        ++mp_oracle_adjudication_count;
        ++mp_ratio_adjudication_count;
        recordMpZ2Result(oracleResult);
        last_z2_update_used_oracle = true;
        last_mp_oracle_z2 = oracleResult.z2;
        if (!oracleResult.trusted())
            throw std::runtime_error(
                std::string("real-Z2 ratio adjudication oracle is not trusted: ") +
                mpZ2StatusName(oracleResult.status) + ": " + oracleResult.message);
        // This is the primary update for an indeterminate accepted ratio, not
        // a periodic checkpoint correction.  The oracle sees the already
        // accepted post-proposal contour and returns its canonical absolute Z2.
        z2_sign = oracleResult.z2;
        return;
    }
    if (phaseFactor.real() < 0.0) z2_sign = -z2_sign;
}

void PfQMC::rightSweep(int capture_boundary, MatType *captured_g,
                       DataType *captured_sign, int *captured_z2_sign,
                       bool *captured_z2_oracle_used,
                       int *captured_oracle_z2)
{
    if (capture_boundary == -1) {
        if (captured_g != nullptr || captured_sign != nullptr || captured_z2_sign != nullptr ||
            captured_z2_oracle_used != nullptr || captured_oracle_z2 != nullptr) {
            throw std::invalid_argument(
                "capture outputs require a requested capture boundary");
        }
    } else {
        if (capture_boundary < 1 || capture_boundary > op_length) {
            throw std::out_of_range("capture boundary must be in [1, op_length]");
        }
        if (captured_g == nullptr || captured_sign == nullptr) {
            throw std::invalid_argument(
                "requested capture boundary requires Green and sign outputs");
        }
    }
    bool captureOccurred = false;
    const auto installMultiprecisionIfNeeded = [this](int boundary, MatType &current) {
        if (!multiprecision_rebuild_callback) return;
        const double condition = fullContourCoreConditionAtBoundary(boundary);
        multiprecision_condition_samples.push_back(condition);
        if (!multiprecision_fallback) return;
        if (!(condition > multiprecision_core_condition_threshold)) return;
        ++multiprecision_proxy_trigger_count;
        MatType rebuilt;
        if (multiprecision_rebuild_callback(boundary, rebuilt)) {
            current.swap(rebuilt);
            ++multiprecision_fallback_count;
        }
    };
    MatType tmp = MatType::Identity(nDim, nDim);
    MatType Aseg = MatType::Identity(nDim, nDim);
    int curSeg = 0;
    DataType signCur;
    for (int l = 0; l < op_length; l++)
    {
        Operator *op = op_array[l];
        const int proposalCount = (adaptive_guard || realZ2Mode())
            ? op->singleFlipProposalCount() : 0;
        if (proposalCount > 0)
        {
            signCur = DataType(1);
            for (int aux = 0; aux < proposalCount; ++aux)
            {
                double uniform = 0.0;
                if (!op->prepareSingleFlip(g, aux, &uniform))
                    throw std::runtime_error("single-field proposal unexpectedly unavailable");
                ++proposal_attempt_count;
                const double minDenominator = op->preparedMinDenominator();
                min_update_denominator = std::min(min_update_denominator, minDenominator);
                DataType ratio = op->preparedRatio();
                const bool dangerous = adaptive_guard &&
                    (minDenominator < guard_threshold ||
                     std::abs(ratio) > guard_ratio_upper);
                if (!dangerous)
                {
                    const bool accept = uniform < std::abs(ratio);
                    const DataType delta=op->finishSingleFlip(g, accept, true);
                    signCur *= delta;
                    updatePhysicalZ2(delta);
                    continue;
                }

                // The proposal identity and uniform stay inside the operator.
                // Rebuilding reads only the current live contour and consumes
                // no RNG.  The boundary is immediately before op_array[l].
                ++adaptive_rebuild_count;
                ++pre_decision_rebuild_count;
                rebuildGreenFromFullContourAtBoundary(l, g);
                installMultiprecisionIfNeeded(l, g);
                ratio = op->recomputePreparedRatio(g);
                const bool accept = uniform < std::abs(ratio);
                if (!accept)
                {
                    // Keep the newly stabilized current-configuration Green.
                    const DataType delta=op->finishSingleFlip(g, false, false);
                    signCur *= delta;
                    updatePhysicalZ2(delta);
                }
                else
                {
                    // Mutate HS/B without the ill-conditioned rank update,
                    // then rebuild once more for the accepted configuration.
                    const DataType delta=op->finishSingleFlip(g, true, false);
                    signCur *= delta;
                    updatePhysicalZ2(delta);
                    rebuildGreenFromFullContourAtBoundary(l, g);
                    installMultiprecisionIfNeeded(l, g);
                    ++adaptive_rebuild_count;
                    ++post_accept_rebuild_count;
                }
            }
        }
        else
        {
            signCur = op->update(g);
            updatePhysicalZ2(signCur);
        }
        this->sign *= signCur;
        max_complex_phase_imag = std::max(max_complex_phase_imag, std::abs(sign.imag()));

        op_array[l]->left_multiply(Aseg, tmp);
        std::swap(Aseg, tmp);
        //%op_length is important, cannot be remove. or else the last segment will not be calculated
        if (need_stabilization[(l + 1) % op_length])
        {
            // auto g2 = g;
            // op_array[i]->left_propagate(g2, tmp);
            // re-evaluate the UDT of current segment
            if (curSeg == 0)
            {
                udtR[curSeg] = UDT(Aseg); // TODO: performance check
            }
            else
            {
                udtR[curSeg] = Aseg * udtR[curSeg - 1];
            }
            Aseg = MatType::Identity(nDim, nDim);
            // re-evaluate the Green's function of next time slice
            if (curSeg == (checkpoints - 1))
            {
                udtR[curSeg].onePlusInv(g);
            }
            else
            {
                g = onePlusInv(udtL[curSeg + 1], udtR[curSeg]);
            }
            // std::cout<<"right g recal "<<(g2-g).norm()<<std::endl;
            curSeg++;
        }
        else
        {
            // no need for stabilization
            // direct propagate
            op_array[l]->left_propagate(g, tmp);
        }
        if (l + 1 == capture_boundary)
        {
            *captured_g = g;
            *captured_sign = sign;
            if (captured_z2_sign != nullptr) *captured_z2_sign = physicalZ2Sign();
            if (captured_z2_oracle_used != nullptr || captured_oracle_z2 != nullptr) {
                if (!realZ2Mode() || !initial_z2_oracle)
                    throw std::runtime_error("capture Z2 adjudication requires real-Z2 mode and an oracle");
                const MpZ2Result oracleResult=initial_z2_oracle(op_array);
                ++mp_oracle_adjudication_count;
                recordMpZ2Result(oracleResult);
                const int oracleZ2=oracleResult.z2;
                if (oracleZ2!=0 && oracleZ2!=z2_sign)
                    ++mp_candidate_mismatch_count;
                if (captured_z2_oracle_used != nullptr) *captured_z2_oracle_used=true;
                if (captured_oracle_z2 != nullptr) *captured_oracle_z2=oracleZ2;
            }
            captureOccurred = true;
        }
    }
    if (capture_boundary != -1 && !captureOccurred) {
        throw std::logic_error("requested capture boundary was not reached");
    }
}

bool PfQMC::recoverLeftGreenAfterOperation(
    const char *source, int boundary, int aux, double structurePreOperation)
{
    const double structureBefore = leftRecoveryStructureResidual(g);
    const double structureDelta = std::abs(structureBefore-structurePreOperation);
    const bool alarm = !leftRecoveryFinite(g) ||
        (structureBefore > left_recovery_structure_threshold &&
         structureDelta > left_recovery_structure_delta_threshold);
    if (!alarm) return false;

    MatType rebuilt;
    rebuildGreenFromFullContourAtBoundary(boundary, rebuilt);
    LeftGreenRecoveryEvent event;
    event.source = source;
    event.boundary = boundary;
    event.aux = aux;
    event.green_error_before = leftRecoveryRelativeError(g, rebuilt);
    event.structure_pre_operation = structurePreOperation;
    event.structure_delta = structureDelta;
    event.structure_before = structureBefore;
    event.structure_after = leftRecoveryStructureResidual(rebuilt);
    g.swap(rebuilt);
    if (event.source == "PROPAGATION") ++left_recovery_propagation_count;
    else ++left_recovery_rank_update_count;
    if (left_recovery_event_hook) left_recovery_event_hook(event);
    return true;
}

DataType PfQMC::leftRecoveryUpdateAtBoundary(Operator *op, int boundary)
{
    const int proposalCount = op->singleFlipProposalCount();
    if (proposalCount <= 0) {
        const DataType delta = op->update(g);
        updatePhysicalZ2(delta);
        return delta;
    }

    DataType signCur(1);
    for (int aux=0; aux<proposalCount; ++aux) {
        const double structurePreOperation = leftRecoveryStructureResidual(g);
        double uniform = 0.0;
        if (!op->prepareSingleFlip(g, aux, &uniform))
            throw std::runtime_error(
                "single-field proposal unexpectedly unavailable");
        ++proposal_attempt_count;
        const DataType ratio = op->preparedRatio();
        const bool accept = uniform < std::abs(ratio);
        const DataType delta = op->finishSingleFlip(g, accept, true);
        signCur *= delta;
        updatePhysicalZ2(delta);
        if (accept)
            recoverLeftGreenAfterOperation(
                "RANK_UPDATE", boundary, aux, structurePreOperation);
    }
    return signCur;
}

DataType PfQMC::realZ2UpdateAtBoundary(Operator *op)
{
    const int proposalCount=op->singleFlipProposalCount();
    if(proposalCount<=0){const DataType delta=op->update(g);updatePhysicalZ2(delta);return delta;}
    DataType signCur(1);
    for(int aux=0;aux<proposalCount;++aux){double uniform=0;if(!op->prepareSingleFlip(g,aux,&uniform))throw std::runtime_error("single-field proposal unexpectedly unavailable");++proposal_attempt_count;const DataType ratio=op->preparedRatio();const bool accept=uniform<std::abs(ratio);const DataType delta=op->finishSingleFlip(g,accept,true);signCur*=delta;updatePhysicalZ2(delta);}
    return signCur;
}

void PfQMC::leftSweep()
{
    MatType tmp = MatType::Identity(nDim, nDim);
    MatType Aseg = MatType::Identity(nDim, nDim);
    int curSeg = checkpoints - 1;
    DataType signCur;
    for (int l = op_length - 1; l > -1; l--)
    {
        const double structurePrePropagation = left_green_recovery
            ? leftRecoveryStructureResidual(g) : 0.0;
        op_array[l]->right_propagate(g, tmp);
        if (left_green_recovery)
            recoverLeftGreenAfterOperation(
                "PROPAGATION", l, -1, structurePrePropagation);
        signCur = left_green_recovery
            ? leftRecoveryUpdateAtBoundary(op_array[l], l)
            : (realZ2Mode() ? realZ2UpdateAtBoundary(op_array[l])
                            : op_array[l]->update(g));
        sign *= signCur;
        max_complex_phase_imag = std::max(max_complex_phase_imag, std::abs(sign.imag()));

        const double structurePreCompletion = left_green_recovery
            ? leftRecoveryStructureResidual(g) : 0.0;
        op_array[l]->right_multiply(Aseg, tmp);
        std::swap(Aseg, tmp);
        if (need_stabilization[l])
        {
            // auto g2 = g;
            Aseg.adjointInPlace();
            // re-evaluate the UDT of current segment
            if (curSeg == (checkpoints - 1))
            {
                udtL[curSeg] = UDT(Aseg); // TODO: performance check
            }
            else
            {
                udtL[curSeg] = Aseg * udtL[curSeg + 1];
            }
            Aseg = MatType::Identity(nDim, nDim);
            // re-evaluate the Green's function of next time slice
            if (curSeg == 0)
            {
                udtL[curSeg].onePlusInv(g);
                g.adjointInPlace();
            }
            else
            {
                g = onePlusInv(udtL[curSeg], udtR[curSeg - 1]);
            }
            curSeg--;
            // std::cout<<"left g recal "<<(g2-g).norm()<<std::endl;
            if (left_green_recovery)
                recoverLeftGreenAfterOperation(
                    "PROPAGATION", l, -1, structurePreCompletion);
        }
    }
}

PfaffianResult PfQMC::getSignRawWithStatus()
{
    const MatType identity = MatType::Identity(nDim, nDim);
    const DataType extraSign = ((nDim / 2) % 2 == 0) ? 1.0 : -1.0;
    UDT A(nDim);
    op_array[0]->stabilizedLeftMultiply(A);
    MatType gNext, gCur;
    DataType signCur, signNext;
    PfaffianResult result;
    result.status = PfaffianStatus::success;
    result.value = DataType(1.0, 0.0);
    signCur = op_array[0]->getSignOfWeight();
    A.onePlusInv(gCur);
    gCur -= identity;
    for (int i = 1; i < op_length; i++)
    {
        op_array[i]->getGreensMat(gNext);
        signNext = op_array[i]->getSignOfWeight();
        // std::cout << gNext << "==== gnext ====\n \n";
        // std::cout << gCur << "==== gcur ====\n \n";
        const PfaffianResult pfaffian =
            pfaffianForSignOfProductWithStatus(gNext, gCur);
        result.min_pivot = std::min(result.min_pivot, pfaffian.min_pivot);
        if (!pfaffian.ok()) {
            result.status = pfaffian.status;
            result.lapack_info = pfaffian.lapack_info;
            result.value = DataType(0.0, 0.0);
            return result;
        }
        // std::cout << "Raw sCur, sNext, sPfaf=" << signCur << " " << signNext << " " << signPfaf << "\n"; 
        signCur = (signCur * signNext * pfaffian.value * extraSign);

        // std::cout << signCur << " sign cur\n";
        if (i == op_length - 1)
            break;
        op_array[i]->stabilizedLeftMultiply(A);
        A.onePlusInv(gCur);
        gCur -= identity;
        // if(std::abs(gCur(0, 0))>1e-4) {
        //     std::cout << " ===== begin writing ==== \n";
        //     int n = A.nDim;
        //     // std::cout << "pf(Acopy) = " << pfaf(2*n, Acopy) << "\n";
        //     std::fstream myfile;
        //     myfile.open("1.dat", std::fstream::out);
        //     // std::cout << " r = " << r << "\n";
        //     // myfile << " r = " << r << "\n";
        //     A.onePlusInv(gCur);
        //     int npre = 20;
        //     std::cout << std::setprecision(npre) << "A.U = \n" << A.U << "\n";
        //     std::cout << std::setprecision(npre) << "A.D = \n" << A.D << "\n";
        //     std::cout << std::setprecision(npre) << "gCur = \n" << gCur << "\n"; 
        //     for(int i=0; i<n; i++) {
        //         for(int j=0; j<n; j++) {
        //             myfile<< std::setprecision(npre) << A.U(i, j) << " ";
        //             std::cout << ".";
        //         }
        //         myfile << "\n";
        //     }

        //     myfile << "\n\n\n";

        //     for(int i=0; i<n; i++) {
        //         myfile<< std::setprecision(npre) << A.D(i) << " ";
        //     }

        //     myfile << "\n\n\n";

        //     for(int i=0; i<n; i++) {
        //         for(int j=0; j<n; j++) {
        //             myfile<< std::setprecision(npre) << A.T(i, j) << " ";
        //         }
        //         myfile << "\n";
        //     }
        //     std::cout << " ===== end writing ==== \n";
        //     myfile.close();
        // }
    }
    result.value = signCur;
    if (!std::isfinite(signCur.real()) || !std::isfinite(signCur.imag())) {
        result.status = PfaffianStatus::nonfinite_pivot;
        result.value = DataType(0.0, 0.0);
    }
    if (realZ2Mode() && result.ok()) {
        result.condition_proxy = rawContourExponentSpan();
        const double magnitude = std::abs(signCur);
        result.phase_reality_error = std::abs(signCur.imag()) /
            std::max(std::abs(signCur.real()), std::numeric_limits<double>::min());
        if (!std::isfinite(result.condition_proxy) || result.condition_proxy > 45.0) {
            result.status = PfaffianStatus::untrusted_condition;
        } else if (!std::isfinite(magnitude) || magnitude == 0.0 ||
                   std::abs(magnitude-1.0) > 1e-8 ||
                   !std::isfinite(result.phase_reality_error) ||
                   result.phase_reality_error > 1e-8 ||
                   std::abs(signCur.real()) < 1e-12) {
            result.status = PfaffianStatus::untrusted_phase;
        }
    }
    return result;
}

double PfQMC::rawContourExponentSpan() const
{
    UDT product(nDim);
    for (Operator *op : op_array) op->stabilizedLeftMultiply(product);
#ifdef PFQMC_SCALE_SAFE_UDT
    if (product.Dexp.size() == 0) return std::numeric_limits<double>::infinity();
    return double(product.Dexp.maxCoeff()) - double(product.Dexp.minCoeff());
#else
    if (product.D.size() == 0) return std::numeric_limits<double>::infinity();
    const double dmax = product.D.maxCoeff();
    const double dmin = product.D.minCoeff();
    if (!(dmin > 0.0) || !std::isfinite(dmax) || !std::isfinite(dmin))
        return std::numeric_limits<double>::infinity();
    return std::log2(dmax/dmin);
#endif
}

PfaffianResult PfQMC::checkRawSignReadOnly()
{
    PfaffianResult raw = getSignRawWithStatus();
    if (!realZ2Mode()) return raw;
    if (raw.untrusted()) {
        ++raw_sign_untrusted_count;
        return raw;
    }
    if (!raw.ok()) return raw;
    ++raw_sign_trusted_count;
    const int rawZ2 = raw.value.real() >= 0.0 ? 1 : -1;
    if (rawZ2 == z2_sign) return raw;
    if (!initial_z2_oracle)
        throw std::runtime_error("trusted raw/transported Z2 mismatch requires an adjudication oracle");
    const MpZ2Result oracleResult = initial_z2_oracle(op_array);
    ++mp_oracle_adjudication_count;
    recordMpZ2Result(oracleResult);
    const int oracleZ2 = oracleResult.z2;
    if (oracleZ2 != 0 && oracleZ2 != z2_sign)
        ++mp_candidate_mismatch_count;
    if (!oracleResult.trusted()) return raw;
    if (oracleZ2 == z2_sign) {
        ++raw_sign_mismatch_count;
        return raw;
    }
    if (oracleZ2 == rawZ2)
        throw std::runtime_error("hard diagnostic failure: MP oracle supports raw Z2 over transported Z2");
    throw std::runtime_error("hard diagnostic failure: MP oracle disagrees with both raw and transported Z2");
}

DataType PfQMC::getSignRaw()
{
    const PfaffianResult result = getSignRawWithStatus();
    if (!result.ok()) {
        throw std::runtime_error(
            std::string("raw Pfaffian sign unavailable: ") +
            pfaffianStatusName(result.status));
    }
    return result.value;
}

// DataType PfQMC::getSign() {
//     //TODO: this current method has fundamental flaws
//     const MatType identity = MatType::Identity(nDim, nDim);
//     UDT A(nDim);
//     op_array[0]->stabilizedLeftMultiply(A);
//     MatType gNext, gCur, gInv, gTemp, t;
//     DataType signCur, signNext, signPfaf;
//     signCur = op_array[0]->getSignOfWeight();
//     A.onePlusInv(gCur);
//     gCur -= identity;
//     for (int i = 1; i < op_length; i++)
//     {
//         op_array[i]->getGreensMat(gNext);
//         // std::cout << gNext << "==== gnext ====\n \n";
//         // std::cout << "det of g=" << gNext.determinant() << "\n";
//         op_array[i]->getGreensMatInv(gInv);
//         MatType t = gNext * gInv;
//         // std::cout << (t - identity).squaredNorm() << "gNext * gInv, bType" << op_array[i]->getType()<<  " \n";
//         // EXPECT_NEAR(, 0.0, 1e-10);
//         signNext = op_array[i]->getSignOfWeight();
//         // std::cout << gNext << "==== gnext ====\n \n";
//         // std::cout << gCur << "==== gcur ====\n \n";
//         // std::cout << gInv << "=gInv\n";
//         gTemp = gInv + gCur;
//         signPfaf = signOfPfaf(gTemp) / op_array[i]->getSignPfGInv();
//         // signPfaf = 1.0;
//         // std::cout << "sCur, sNext, sPfaf=" << signCur << " " << signNext << " " << signPfaf << "\n\n"; 
//         signCur = (signCur * signNext * signPfaf);
//         // std::cout << signCur << " sign cur\n";
//         if (i == op_length - 1)
//             break;
//         op_array[i]->stabilizedLeftMultiply(A);
//         A.onePlusInv(gCur);
//         gCur -= identity;
//     }
//     return signCur;
// }
