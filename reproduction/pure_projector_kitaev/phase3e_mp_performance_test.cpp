#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "gaussian_trial_state.h"
#include "pure_projector_mp.h"

namespace {
void require(bool value,const char *message){if(!value)throw std::runtime_error(message);}

MatType canonicalPhi(int modes){MatType p=MatType::Zero(2*modes,modes);
    for(int i=0;i<modes;++i){p(i,i)=1/std::sqrt(2.0);
        p(modes+i,i)=DataType(0,1/std::sqrt(2.0));}return p;}

MatType boostFactor(double angle){MatType result=MatType::Identity(2,2);
    result(0,0)=std::cosh(angle);result(1,1)=std::cosh(angle);
    result(0,1)=DataType(0.0,std::sinh(angle));
    result(1,0)=DataType(0.0,-std::sinh(angle));return result;}

void cacheMatchesLegacy(){
    const GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(1));
    const MatType a=MatType::Identity(2,2),b=MatType::Identity(2,2);
    std::vector<PureProjectorSlice> slices;
    for(int i=0;i<24;++i)slices.emplace_back(i%2?a:b,1.0,"repeated");
    const int index=11;const MatType candidate=boostFactor(0.02);
    PureMpProposalOptions legacy;legacy.enable_operator_cache=false;
    PureMpProposalOptions cached=legacy;cached.enable_operator_cache=true;
    const PureMpProposalResult before=pureProjectorMpSameProposal(
        trial,slices,index,candidate,1.0,legacy);
    const PureMpProposalResult after=pureProjectorMpSameProposal(
        trial,slices,index,candidate,1.0,cached);
    if(!before.ok()||!after.ok())std::cerr<<"legacy="<<pureMpProposalStatusName(before.status)
        <<":"<<before.message<<" cached="<<pureMpProposalStatusName(after.status)
        <<":"<<after.message<<'\n';
    require(before.ok()&&after.ok(),"legacy/cached MP evaluation failed");
    require(before.ratio_z2==after.ratio_z2,"cache changed ratio Z2");
    require(std::abs(before.ratio-after.ratio)<=1e-13*std::max(1.0,std::abs(before.ratio)),
            "cache changed MP ratio");
    require(after.profile.operator_cache_hits>0,"operator cache recorded no hits");
    require(after.profile.operator_cache_misses<after.profile.operator_requests,
            "operator cache did not reduce conversions");
    require(after.profile.sparse_apply_count>0,"sparse local operator path was unused");
    require(after.profile.precision_160_seconds>0&&after.profile.precision_320_seconds>0,
            "MP precision-stage timers were not recorded");
}

void genericComplexStaysLegacy(){
    const GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(1));
    std::vector<PureProjectorSlice> slices(8,
        PureProjectorSlice(MatType::Identity(2,2),1.0,"generic"));
    PureMpProposalOptions expected;expected.real_z2=false;
    expected.enable_operator_cache=false;
    PureMpProposalOptions production=expected;production.enable_operator_cache=true;
    const MatType candidate=boostFactor(0.01);
    const auto legacy=pureProjectorMpSameProposal(trial,slices,3,candidate,1.0,expected);
    const auto actual=pureProjectorMpSameProposal(trial,slices,3,candidate,1.0,production);
    require(legacy.ok()&&actual.ok(),"generic-complex legacy evaluation failed");
    require(legacy.ratio==actual.ratio&&legacy.ratio_z2==actual.ratio_z2,
            "generic-complex result changed");
    require(actual.profile.canonical_input_builds==0&&
            actual.profile.operator_cache_hits==0,"generic-complex entered cache path");
}
}

int main(){try{cacheMatchesLegacy();genericComplexStaysLegacy();
    std::cout<<"PASS phase3e_mp_operator_cache\n";return 0;}
catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
