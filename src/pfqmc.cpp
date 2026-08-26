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

PfQMC::PfQMC(Spinless_tV *walker, int _stb)
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
    if (!initialSign.ok()) {
        throw std::runtime_error(
            std::string("initial raw sign unavailable: ") +
            pfaffianStatusName(initialSign.status));
    }
    sign = initialSign.value;
}

void PfQMC::rightSweep(int capture_boundary, MatType *captured_g,
                       DataType *captured_sign)
{
    if (capture_boundary == -1) {
        if (captured_g != nullptr || captured_sign != nullptr) {
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
        const int proposalCount = adaptive_guard ? op->singleFlipProposalCount() : 0;
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
                const bool dangerous = minDenominator < guard_threshold ||
                                       std::abs(ratio) > guard_ratio_upper;
                if (!dangerous)
                {
                    const bool accept = uniform < std::abs(ratio);
                    signCur *= op->finishSingleFlip(g, accept, true);
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
                    signCur *= op->finishSingleFlip(g, false, false);
                }
                else
                {
                    // Mutate HS/B without the ill-conditioned rank update,
                    // then rebuild once more for the accepted configuration.
                    signCur *= op->finishSingleFlip(g, true, false);
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
        }
        this->sign *= signCur;

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
    if (proposalCount <= 0) return op->update(g);

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
        signCur *= op->finishSingleFlip(g, accept, true);
        if (accept)
            recoverLeftGreenAfterOperation(
                "RANK_UPDATE", boundary, aux, structurePreOperation);
    }
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
            : op_array[l]->update(g);
        sign *= signCur;

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
    return result;
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
