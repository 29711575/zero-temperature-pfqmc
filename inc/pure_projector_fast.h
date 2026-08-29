#ifndef PURE_PROJECTOR_FAST_H
#define PURE_PROJECTOR_FAST_H

#include "pure_projector_stack.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

enum class PureBranch { Ket, Bra };
enum class PureFastRunMode { AuditLockstep, FastStrict };
enum class PureFastRatioStatus {
    success,
    invalid_proposal,
    rank_failure,
    overlap_untrusted,
    pfaffian_untrusted,
    nonfinite_ratio,
    complex_ratio,
    zero_ratio,
    reference_failure,
    terminated
};

struct PureSliceLocation {
    PureBranch branch = PureBranch::Ket;
    int slice = -1;
    int factor = -1;
    int bond = -1;
    int aux = -1;
};

struct PureFastConfiguration {
    std::vector<PureProjectorSlice> slices;
    std::vector<int> hs_fields;
    std::vector<PureSliceLocation> locations;
};

struct PureFastProposal {
    int index = -1;
    int new_hs = 0;
    MatType new_factor;
    DataType new_eta = DataType(1.0,0.0);
    double uniform = 0.0;
};

struct PureProposalSnapshot {
    PureBranch branch = PureBranch::Ket;
    int slice = -1;
    int factor = -1;
    int bond = -1;
    int aux = -1;
    int old_hs = 0;
    int new_hs = 0;
    std::uint64_t configuration_hash = 0;
    DataType fast_ratio = DataType(0.0,0.0);
    double uniform = 0.0;
};

struct PureFastRatioResult {
    PureFastRatioStatus status = PureFastRatioStatus::invalid_proposal;
    DataType ratio = DataType(0.0,0.0);
    DataType slow_ratio = DataType(0.0,0.0);
    int low_rank = 0;
    double overlap_rcond = 0.0;
    double solve_residual = std::numeric_limits<double>::infinity();
    PfaffianStatus pfaffian_status = PfaffianStatus::success;
    double relative_reference_error = std::numeric_limits<double>::infinity();
    double green_reference_error = std::numeric_limits<double>::infinity();
    MatType fast_updated_green;
    bool trust_alarm = false;
    bool used_reference = false;
    bool ok() const { return status == PureFastRatioStatus::success; }
};

struct PureFastProposalResult {
    PureProposalSnapshot snapshot;
    PureFastRatioResult ratio;
    bool accepted = false;
    bool terminated = false;
    bool reference_accepted = false;
    bool z2_reference_mismatch = false;
    double fast_seconds = 0.0;
    double slow_seconds = 0.0;
};

struct PureFastOptions {
    PureProjectorWeightMode weight_mode = PureProjectorWeightMode::RealZ2;
    double rank_tolerance = 1e-12;
    double minimum_overlap_rcond = 1e-12;
    double residual_tolerance = 1e-10;
    double reality_tolerance = 1e-10;
    double zero_tolerance = 1e-14;
    double decision_margin_tolerance = 1e-12;
};

inline std::uint64_t pureFastConfigurationHash(const PureFastConfiguration &configuration) {
    std::uint64_t hash=1469598103934665603ULL;
    auto mix=[&](std::uint64_t value){hash^=value;hash*=1099511628211ULL;};
    for(int field:configuration.hs_fields)mix(std::uint64_t(std::int64_t(field)));
    for(const auto &slice:configuration.slices)
        for(int c=0;c<slice.matrix.cols();++c)for(int r=0;r<slice.matrix.rows();++r){
            std::uint64_t re=0,im=0;double a=slice.matrix(r,c).real(),b=slice.matrix(r,c).imag();
            std::memcpy(&re,&a,sizeof(double));std::memcpy(&im,&b,sizeof(double));mix(re);mix(im);
        }
    return hash;
}

inline PureFastRatioResult pureFastLocalRatioAtCut(
    const GaussianTrialState &trial,
    const std::vector<PureProjectorSlice> &slices,
    int index,const MatType &newFactor,DataType newEta,
    const MatType &rightBefore,const MatType &leftAfter,
    const PureFastOptions &options) {
    PureFastRatioResult result;
    if(index<0||index>=int(slices.size())||newFactor.rows()!=trial.Phi.rows()||
       newFactor.cols()!=trial.Phi.rows()||!pureProjectorMatrixFinite(newFactor))return result;
    if(rightBefore.rows()!=trial.Phi.rows()||leftAfter.rows()!=trial.Phi.rows()||
       rightBefore.cols()!=trial.Phi.cols()||leftAfter.cols()!=trial.Phi.cols())return result;
    const MatType oldFactor=slices[index].matrix;
    const MatType oldRight=oldFactor*rightBefore;
    const MatType deltaRight=(newFactor-oldFactor)*rightBefore;
    Eigen::JacobiSVD<MatType> deltaSvd(deltaRight,Eigen::ComputeThinU|Eigen::ComputeThinV);
    deltaSvd.setThreshold(options.rank_tolerance);
    const int rank=deltaSvd.rank();
    result.low_rank=rank;
    if(rank==0){
        result.ratio=newEta/slices[index].eta;
        if(!std::isfinite(result.ratio.real())||!std::isfinite(result.ratio.imag())){
            result.status=PureFastRatioStatus::nonfinite_ratio;return result;
        }
        if(std::abs(result.ratio)<=options.zero_tolerance){
            result.status=PureFastRatioStatus::zero_ratio;return result;
        }
        if(options.weight_mode==PureProjectorWeightMode::RealZ2&&
           std::abs(result.ratio.imag())>options.reality_tolerance*
               std::max(1.0,std::abs(result.ratio.real()))){
            result.status=PureFastRatioStatus::complex_ratio;return result;
        }
        result.status=PureFastRatioStatus::success;
        result.fast_updated_green=pureProjectorGreenThinQr(oldRight,leftAfter).green;
        result.overlap_rcond=1.0;result.solve_residual=0.0;return result;
    }
    MatType u=deltaSvd.matrixU().leftCols(rank);
    MatType c=deltaSvd.singularValues().head(rank).asDiagonal()*
              deltaSvd.matrixV().leftCols(rank).adjoint();
    const MatType overlap=leftAfter.adjoint()*oldRight;
    Eigen::JacobiSVD<MatType> overlapSvd(overlap,Eigen::ComputeThinU|Eigen::ComputeThinV);
    overlapSvd.setThreshold(options.rank_tolerance);
    if(overlapSvd.rank()!=overlap.rows()){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    const double largest=overlapSvd.singularValues()(0);
    const double smallest=overlapSvd.singularValues()(overlap.rows()-1);
    result.overlap_rcond=largest>0?smallest/largest:0;
    if(result.overlap_rcond<options.minimum_overlap_rcond){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    const MatType a=leftAfter.adjoint()*u;
    const MatType x=overlapSvd.solve(a);
    const MatType y=overlapSvd.solve(leftAfter.adjoint());
    result.solve_residual=std::max(
        (overlap*x-a).norm()/std::max(1.0,a.norm()),
        (overlap*y-leftAfter.adjoint()).norm()/std::max(1.0,leftAfter.norm()));
    if(!std::isfinite(result.solve_residual)||result.solve_residual>options.residual_tolerance){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    const MatType k=MatType::Identity(rank,rank)+c*x;
    Eigen::FullPivLU<MatType> klu(k);
    if(!klu.isInvertible()){result.status=PureFastRatioStatus::zero_ratio;return result;}
    const DataType etaRatio=newEta/slices[index].eta;
    const DataType ratioSquared=etaRatio*etaRatio*k.determinant();
    if(!std::isfinite(ratioSquared.real())||!std::isfinite(ratioSquared.imag())){
        result.status=PureFastRatioStatus::nonfinite_ratio;return result;
    }
    if(std::abs(ratioSquared)<=options.zero_tolerance){
        result.status=PureFastRatioStatus::zero_ratio;return result;
    }
    DataType ratio=std::sqrt(ratioSquared);

    Eigen::FullPivLU<MatType> oldTranspose(oldFactor.transpose());
    if(!oldTranspose.isInvertible()){
        result.status=PureFastRatioStatus::rank_failure;return result;
    }
    MatType transform=oldTranspose.solve(newFactor.transpose()).transpose();
    MatType onePlus=MatType::Identity(transform.rows(),transform.cols())+transform;
    Eigen::FullPivLU<MatType> cayleySolve(onePlus);
    if(!cayleySolve.isInvertible()){
        result.status=PureFastRatioStatus::pfaffian_untrusted;return result;
    }
    MatType transformGreen=cayleySolve.solve(2.0*MatType::Identity(
        transform.rows(),transform.cols()))-MatType::Identity(transform.rows(),transform.cols());
    PureProjectorGreenResult oldGreen=pureProjectorGreenThinQr(oldRight,leftAfter);
    if(!oldGreen.ok()){result.status=PureFastRatioStatus::overlap_untrusted;return result;}
    PfaffianResult pf=pfaffianForSignOfProductWithStatus(transformGreen,oldGreen.green);
    result.pfaffian_status=pf.status;
    if(!pf.ok()){result.status=PureFastRatioStatus::pfaffian_untrusted;return result;}
    const DataType extraSign=((trial.Phi.cols()%2)==0)?DataType(1.0):DataType(-1.0);
    const DataType target=pureProjectorUnitPhase(etaRatio*pf.value*extraSign);
    if((std::conj(target)*ratio).real()<0)ratio=-ratio;
    if(options.weight_mode==PureProjectorWeightMode::RealZ2&&
       std::abs(ratio.imag())>options.reality_tolerance*std::max(1.0,std::abs(ratio.real()))){
        result.status=PureFastRatioStatus::complex_ratio;return result;
    }
    if(!std::isfinite(ratio.real())||!std::isfinite(ratio.imag())){
        result.status=PureFastRatioStatus::nonfinite_ratio;return result;
    }
    result.ratio=ratio;
    const MatType z=klu.solve(c*y);
    result.fast_updated_green=oldGreen.green-2.0*(u-oldRight*x)*z;
    if(!pureProjectorMatrixFinite(result.fast_updated_green)){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    result.status=PureFastRatioStatus::success;
    return result;
}

inline PureFastRatioResult pureFastLocalRatio(
    const GaussianTrialState &trial,
    const std::vector<PureProjectorSlice> &slices,
    int index,const MatType &newFactor,DataType newEta,
    const PureFastOptions &options) {
    PureFastRatioResult result;
    PureProjectorStackManager stack(trial,slices,std::max(1,int(slices.size())));
    if(!stack.moveToCut(index)){result.status=PureFastRatioStatus::overlap_untrusted;return result;}
    const MatType rightBefore=stack.phiRight();
    if(!stack.moveToCut(index+1)){result.status=PureFastRatioStatus::overlap_untrusted;return result;}
    return pureFastLocalRatioAtCut(
        trial,slices,index,newFactor,newEta,rightBefore,stack.phiLeft(),options);
}

class PureProjectorFastWalker {
public:
    PureProjectorFastWalker(const GaussianTrialState &trial,
                            PureFastConfiguration configuration,int blockSize,
                            PureFastRunMode mode,
                            PureFastOptions options=PureFastOptions())
        : trial_(trial),configuration_(std::move(configuration)),block_size_(blockSize),
          mode_(mode),options_(options) {
        if(configuration_.slices.size()!=configuration_.hs_fields.size()||
           configuration_.slices.size()!=configuration_.locations.size())
            throw std::invalid_argument("fast walker configuration size mismatch");
        PureProjectorWeightOptions weightOptions;
        weightOptions.mode=PureProjectorWeightMode::GenericComplex;
        evaluator_=PureProjectorWeightEvaluator(weightOptions);
        current_=evaluator_.evaluate(trial_,configuration_.slices);
        if(!current_.ok())throw std::runtime_error("fast walker initial weight failed");
        z2_sign_=current_.z2_sign==0?1:current_.z2_sign;
        manager_.reset(new PureProjectorStackManager(
            trial_,configuration_.slices,block_size_));
        if(!manager_->ok())throw std::runtime_error("fast walker initial stack failed");
    }

    std::uint64_t configurationHash() const {return pureFastConfigurationHash(configuration_);}
    PureProjectorGreenResult measurementGreen() {
        manager_->moveToCut(int(configuration_.slices.size()/2));return manager_->green();
    }
    int z2Sign() const {return z2_sign_;}
    const PureFastConfiguration &configuration() const{return configuration_;}
    const PureProjectorWeightResult &currentWeight() const{return current_;}

    PureFastProposalResult propose(const PureFastProposal &proposal) {
        PureFastProposalResult output;
        if(terminated_){output.terminated=true;output.ratio.status=PureFastRatioStatus::terminated;return output;}
        if(proposal.index<0||proposal.index>=int(configuration_.slices.size())||
           (proposal.new_hs!=-1&&proposal.new_hs!=1)||proposal.uniform<0||proposal.uniform>=1){
            output.ratio.status=PureFastRatioStatus::invalid_proposal;return output;
        }
        const auto &location=configuration_.locations[proposal.index];
        output.snapshot.branch=location.branch;output.snapshot.slice=location.slice;
        output.snapshot.factor=location.factor;output.snapshot.bond=location.bond;
        output.snapshot.aux=location.aux;output.snapshot.old_hs=configuration_.hs_fields[proposal.index];
        output.snapshot.new_hs=proposal.new_hs;output.snapshot.configuration_hash=configurationHash();
        output.snapshot.uniform=proposal.uniform;
        const auto fastStart=std::chrono::steady_clock::now();
        if(!manager_->moveToCut(proposal.index)){
            output.ratio.status=PureFastRatioStatus::overlap_untrusted;
        }else{
            const MatType rightBefore=manager_->phiRight();
            if(!manager_->moveToCut(proposal.index+1))
                output.ratio.status=PureFastRatioStatus::overlap_untrusted;
            else output.ratio=pureFastLocalRatioAtCut(
                trial_,configuration_.slices,proposal.index,proposal.new_factor,
                proposal.new_eta,rightBefore,manager_->phiLeft(),options_);
        }
        output.fast_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-fastStart).count();
        output.snapshot.fast_ratio=output.ratio.ratio;

        PureFastConfiguration candidate=configuration_;
        candidate.slices[proposal.index].matrix=proposal.new_factor;
        candidate.slices[proposal.index].eta=proposal.new_eta;
        candidate.hs_fields[proposal.index]=proposal.new_hs;
        const double fastProbability=output.ratio.ok()?std::min(1.0,std::abs(output.ratio.ratio)):0.0;
        output.ratio.trust_alarm=!output.ratio.ok()||
            std::abs(proposal.uniform-fastProbability)<options_.decision_margin_tolerance;
        const bool needReference=mode_==PureFastRunMode::AuditLockstep||output.ratio.trust_alarm;
        PureProjectorWeightResult candidateWeight;
        if(needReference){
            const auto slowStart=std::chrono::steady_clock::now();
            candidateWeight=evaluator_.evaluate(trial_,candidate.slices);
            output.slow_seconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-slowStart).count();
            output.ratio.used_reference=true;
            if(!candidateWeight.ok()||std::abs(current_.weight)==0){
                output.ratio.status=PureFastRatioStatus::reference_failure;
                terminated_=true;output.terminated=true;return output;
            }
            output.ratio.slow_ratio=candidateWeight.weight/current_.weight;
            if(output.ratio.ok())output.ratio.relative_reference_error=
                std::abs(output.ratio.ratio-output.ratio.slow_ratio)/
                std::max(1e-14,std::abs(output.ratio.slow_ratio));
            const double slowProbability=std::min(1.0,std::abs(output.ratio.slow_ratio));
            output.reference_accepted=proposal.uniform<slowProbability;
            if(output.ratio.fast_updated_green.size()!=0){
                PureProjectorStackManager candidateStack(
                    trial_,candidate.slices,block_size_);
                if(candidateStack.moveToCut(proposal.index+1)){
                    const PureProjectorGreenResult referenceGreen=candidateStack.green();
                    if(referenceGreen.ok())output.ratio.green_reference_error=
                        (output.ratio.fast_updated_green-referenceGreen.green).norm()/
                        std::max(1.0,referenceGreen.green.norm());
                }
            }
            if(output.ratio.trust_alarm){
                output.ratio.ratio=output.ratio.slow_ratio;
                output.ratio.status=PureFastRatioStatus::success;
            }
        }
        if(!output.ratio.ok())return output;
        const double probability=std::min(1.0,std::abs(output.ratio.ratio));
        output.accepted=proposal.uniform<probability;
        if(mode_==PureFastRunMode::AuditLockstep&&output.accepted!=output.reference_accepted){
            output.ratio.status=PureFastRatioStatus::reference_failure;
            terminated_=true;output.terminated=true;return output;
        }
        if(output.accepted){
            const int ratioSign=output.ratio.ratio.real()>=0?1:-1;
            if(options_.weight_mode==PureProjectorWeightMode::RealZ2)z2_sign_*=ratioSign;
            configuration_=std::move(candidate);
            if(needReference)current_=candidateWeight;
            else{
                current_.log_abs_weight+=std::log(std::abs(output.ratio.ratio));
                current_.complex_phase*=pureProjectorUnitPhase(output.ratio.ratio);
                current_.weight=std::exp(current_.log_abs_weight)*current_.complex_phase;
                current_.z2_sign=z2_sign_;
            }
            manager_.reset(new PureProjectorStackManager(
                trial_,configuration_.slices,block_size_));
            if(!manager_->ok()){terminated_=true;output.terminated=true;}
        }
        if(needReference&&candidateWeight.ok()){
            const int expected=output.accepted?
                (candidateWeight.z2_sign==0?z2_sign_:candidateWeight.z2_sign):
                (current_.z2_sign==0?z2_sign_:current_.z2_sign);
            output.z2_reference_mismatch=(expected!=z2_sign_);
        }
        return output;
    }

private:
    GaussianTrialState trial_;
    PureFastConfiguration configuration_;
    int block_size_;
    PureFastRunMode mode_;
    PureFastOptions options_;
    PureProjectorWeightEvaluator evaluator_;
    PureProjectorWeightResult current_;
    std::unique_ptr<PureProjectorStackManager> manager_;
    int z2_sign_=1;
    bool terminated_=false;
};

#endif
