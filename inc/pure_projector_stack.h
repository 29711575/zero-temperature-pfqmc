#ifndef PURE_PROJECTOR_STACK_H
#define PURE_PROJECTOR_STACK_H

#include "gaussian_trial_state.h"
#include "pure_projector_green.h"
#include "pure_projector_weight.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

enum class PureStackStatus {
    success,
    invalid_configuration,
    factor_solve_failure,
    subspace_failure,
    overlap_untrusted
};

struct PureStackCheckpoint {
    int cut = 0;
    MatType q_right;
    MatType q_left;
};

class PureProjectorStackManager {
public:
    PureProjectorStackManager(const GaussianTrialState &trial,
                              std::vector<PureProjectorSlice> slices,
                              int blockSize,
                              PureProjectorOptions options = PureProjectorOptions())
        : trial_(trial), slices_(std::move(slices)), block_size_(blockSize),
          options_(options) {
        if (block_size_ <= 0) {
            status_ = PureStackStatus::invalid_configuration;
            return;
        }
        for (const auto &slice : slices_) {
            if (slice.matrix.rows() != trial_.Phi.rows() ||
                slice.matrix.cols() != trial_.Phi.rows() ||
                !pureProjectorMatrixFinite(slice.matrix)) {
                status_ = PureStackStatus::invalid_configuration;
                return;
            }
        }
        buildCheckpoints();
        rebuildAtCut(0);
    }

    bool ok() const { return status_ == PureStackStatus::success; }
    PureStackStatus status() const { return status_; }
    int cut() const { return cut_; }
    int length() const { return int(slices_.size()); }
    int blockSize() const { return block_size_; }
    int checkpointCount() const { return int(checkpoints_.size()); }
    const MatType &phiRight() const { return phi_right_; }
    const MatType &phiLeft() const { return phi_left_; }
    const std::vector<PureProjectorSlice> &slices() const { return slices_; }
    const GaussianTrialState &trial() const { return trial_; }

    bool acceptFactorAtCut(int index,const PureProjectorSlice &replacement,
                           const MatType &rightBefore) {
        if(!ok()||index<0||index>=length()||cut_!=index+1||
           replacement.matrix.rows()!=trial_.Phi.rows()||
           replacement.matrix.cols()!=trial_.Phi.rows()||
           rightBefore.rows()!=trial_.Phi.rows()||rightBefore.cols()!=trial_.Phi.cols()||
           !pureProjectorMatrixFinite(replacement.matrix))return false;
        slices_[index]=replacement;
        phi_right_=replacement.matrix*rightBefore;
        if(!stabilizeCurrent())return false;
        const PureProjectorGreenResult checked=green();
        if(!checked.ok()){status_=PureStackStatus::overlap_untrusted;return false;}
        return true;
    }

    PureProjectorGreenResult green() const {
        return pureProjectorGreenThinQr(phi_right_, phi_left_, options_);
    }

    bool moveToCut(int target) {
        if (target < 0 || target > length()) return false;
        while (cut_ < target) if (!forwardOne()) return false;
        while (cut_ > target) if (!backwardOne()) return false;
        return true;
    }

    bool forwardOne() {
        if (!ok() || cut_ >= length()) return false;
        const MatType &factor = slices_[cut_].matrix;
        MatType nextRight = factor * phi_right_;
        Eigen::FullPivLU<MatType> solve(factor.adjoint());
        if (!solve.isInvertible()) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        MatType nextLeft = solve.solve(phi_left_);
        if (!pureProjectorMatrixFinite(nextRight) || !pureProjectorMatrixFinite(nextLeft)) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        ++cut_;
        phi_right_.swap(nextRight);
        phi_left_.swap(nextLeft);
        if ((cut_ % block_size_) == 0 || cut_ == length()) {
            if (!stabilizeCurrent()) return false;
        }
        const PureProjectorGreenResult checked = green();
        if (!checked.ok()) {
            status_ = PureStackStatus::overlap_untrusted;
            return false;
        }
        return true;
    }

    bool backwardOne() {
        if (!ok() || cut_ <= 0) return false;
        const MatType &factor = slices_[cut_ - 1].matrix;
        Eigen::FullPivLU<MatType> solve(factor);
        if (!solve.isInvertible()) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        MatType nextRight = solve.solve(phi_right_);
        MatType nextLeft = factor.adjoint() * phi_left_;
        if (!pureProjectorMatrixFinite(nextRight) || !pureProjectorMatrixFinite(nextLeft)) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        --cut_;
        phi_right_.swap(nextRight);
        phi_left_.swap(nextLeft);
        if (cut_ == 0) {
            // Re-anchor the endpoint gauge deterministically so the Phase 3A
            // round-trip invariant remains exact.
            return rebuildAtCut(0);
        }
        if ((cut_ % block_size_) == 0) {
            if (!stabilizeCurrent()) return false;
        }
        const PureProjectorGreenResult checked = green();
        if (!checked.ok()) {
            status_ = PureStackStatus::overlap_untrusted;
            return false;
        }
        return true;
    }

    bool forwardToEnd() { return moveToCut(length()); }
    bool backwardToStart() { return moveToCut(0); }
    bool rebuildCurrent() { return rebuildAtCut(cut_); }
    bool rebuildToCut(int target) {
        if(target<0||target>length())return false;
        return rebuildAtCut(target);
    }

    std::uint64_t configurationHash() const {
        std::uint64_t hash = 1469598103934665603ULL;
        auto mix = [&hash](std::uint64_t value) {
            hash ^= value; hash *= 1099511628211ULL;
        };
        for (const auto &slice : slices_) {
            for (int col=0; col<slice.matrix.cols(); ++col)
                for (int row=0; row<slice.matrix.rows(); ++row) {
                    std::uint64_t re=0,im=0;
                    const double r=slice.matrix(row,col).real();
                    const double i=slice.matrix(row,col).imag();
                    std::memcpy(&re,&r,sizeof(double));
                    std::memcpy(&im,&i,sizeof(double));
                    mix(re); mix(im);
                }
        }
        return hash;
    }

private:
    GaussianTrialState trial_;
    std::vector<PureProjectorSlice> slices_;
    int block_size_ = 1;
    PureProjectorOptions options_;
    PureStackStatus status_ = PureStackStatus::success;
    int cut_ = 0;
    MatType phi_right_;
    MatType phi_left_;
    std::vector<PureStackCheckpoint> checkpoints_;

    bool stabilizeCurrent() {
        const ThinQrResult right = thinQrSubspace(phi_right_, options_);
        const ThinQrResult left = thinQrSubspace(phi_left_, options_);
        if (!right.ok() || !left.ok()) {
            status_ = PureStackStatus::subspace_failure;
            return false;
        }
        phi_right_ = right.q;
        phi_left_ = left.q;
        return true;
    }

    bool rebuildAtCut(int target) {
        MatType right = trial_.Phi;
        for (int index=0; index<target; ++index) {
            right = slices_[index].matrix * right;
            if (((index+1)%block_size_)==0) {
                ThinQrResult qr=thinQrSubspace(right,options_);
                if(!qr.ok()){status_=PureStackStatus::subspace_failure;return false;}
                right=qr.q;
            }
        }
        MatType left = trial_.Phi;
        int propagated=0;
        for (int index=length()-1; index>=target; --index) {
            left = slices_[index].matrix.adjoint() * left;
            ++propagated;
            if ((propagated%block_size_)==0) {
                ThinQrResult qr=thinQrSubspace(left,options_);
                if(!qr.ok()){status_=PureStackStatus::subspace_failure;return false;}
                left=qr.q;
            }
        }
        phi_right_.swap(right);
        phi_left_.swap(left);
        cut_=target;
        status_=PureStackStatus::success;
        const PureProjectorGreenResult checked=green();
        if(!checked.ok()){status_=PureStackStatus::overlap_untrusted;return false;}
        return true;
    }

    void buildCheckpoints() {
        checkpoints_.clear();
        std::vector<int> cuts;
        for (int cut=0; cut<=length(); cut+=block_size_) cuts.push_back(cut);
        if (cuts.empty() || cuts.back()!=length()) cuts.push_back(length());
        std::vector<MatType> rights(cuts.size()),lefts(cuts.size());

        MatType right=trial_.Phi;
        std::size_t next=0;
        for(int cut=0;cut<=length();++cut){
            if(next<cuts.size()&&cut==cuts[next]){
                ThinQrResult qr=thinQrSubspace(right,options_);
                if(!qr.ok()){status_=PureStackStatus::subspace_failure;return;}
                rights[next]=qr.q;right=qr.q;++next;
            }
            if(cut<length())right=slices_[cut].matrix*right;
        }

        MatType left=trial_.Phi;
        int previous=int(cuts.size())-1;
        for(int cut=length();cut>=0;--cut){
            if(previous>=0&&cut==cuts[previous]){
                ThinQrResult ql=thinQrSubspace(left,options_);
                if(!ql.ok()){status_=PureStackStatus::subspace_failure;return;}
                lefts[previous]=ql.q;left=ql.q;--previous;
            }
            if(cut>0)left=slices_[cut-1].matrix.adjoint()*left;
        }
        for(std::size_t i=0;i<cuts.size();++i)
            checkpoints_.push_back({cuts[i],rights[i],lefts[i]});
    }
};

#endif
