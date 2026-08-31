#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "gaussian_trial_state.h"
#include "pure_projector_mp.h"

namespace {
void require(bool value,const char*message){if(!value)throw std::runtime_error(message);}

MatType canonicalPhi(int modes){MatType p=MatType::Zero(2*modes,modes);
    for(int i=0;i<modes;++i){p(i,i)=1/std::sqrt(2.0);
        p(modes+i,i)=DataType(0,1/std::sqrt(2.0));}return p;}

MatType localFactor(int dimension,int a,int b,double angle){
    MatType result=MatType::Identity(dimension,dimension);
    result(a,a)=std::cosh(angle);result(b,b)=std::cosh(angle);
    result(a,b)=DataType(0,std::sinh(angle));
    result(b,a)=DataType(0,-std::sinh(angle));return result;}

void compare(const PureMpProposalResult&a,const PureMpProposalResult&b){
    require(a.ok()&&b.ok(),"new/legacy oracle failed");
    require(a.status==b.status&&a.precision_digits==b.precision_digits,
            "status or precision changed");
    require(a.ratio_z2==b.ratio_z2,"Z2 changed");
    require(std::abs(a.ratio-b.ratio)<=1e-12*std::max(1.0,std::abs(b.ratio)),
            "ratio changed");
    require(std::abs(a.pre_endpoint_residual-b.pre_endpoint_residual)<=1e-12&&
            std::abs(a.post_endpoint_residual-b.post_endpoint_residual)<=1e-12,
            "endpoint residual changed");
}

void lifecycleAndStaleGate(){
    const GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(2));
    std::vector<PureProjectorSlice>slices;
    for(int i=0;i<96;++i)slices.emplace_back(
        localFactor(4,i%2,2+i%2,(i%3-1)*0.003),1.0,"ordered");
    PureMpSubspaceCache cache;
    PureMpProposalOptions optimized,legacy;
    optimized.enable_subspace_cache=true;legacy.enable_subspace_cache=false;
    std::uint64_t hash=pureProjectorMpConfigurationHash(slices);
    std::uint64_t invalidated=0;
    for(int step=0;step<8;++step){const int index=8+step*10;
        const MatType candidate=localFactor(4,index%2,2+index%2,0.004*(step%2?1:-1));
        const auto fresh=pureProjectorMpSameProposal(
            trial,slices,index,candidate,1.0,optimized,&cache,hash);
        const auto reference=pureProjectorMpSameProposal(
            trial,slices,index,candidate,1.0,legacy,nullptr,hash);
        compare(fresh,reference);
        if(step%2==0){const std::uint64_t oldHash=hash;
            slices[index].matrix=candidate;hash=pureProjectorMpConfigurationHash(slices);
            const PureMpPerformanceProfile invalidation=
                cache.acceptedUpdate(oldHash,hash,index,slices);
            invalidated+=invalidation.subspace_cache_invalidations;}
    }
    require(invalidated>0,"accepted updates never invalidated affected checkpoints");
    require(cache.profile().subspace_restore_count>0,"no checkpoint restore occurred");
    require(cache.profile().subspace_saved_factor_count>0,"no propagation was saved");

    // Deliberately mutate without notifying the cache.  The key mismatch must
    // reject the stale entries and recompute safely from the boundary.
    slices[7].matrix=localFactor(4,1,3,0.011);
    hash=pureProjectorMpConfigurationHash(slices);
    const MatType candidate=localFactor(4,0,2,-0.007);
    const auto recovered=pureProjectorMpSameProposal(
        trial,slices,47,candidate,1.0,optimized,&cache,hash);
    const auto reference=pureProjectorMpSameProposal(
        trial,slices,47,candidate,1.0,legacy,nullptr,hash);
    compare(recovered,reference);
    require(cache.profile().subspace_stale_rejection_count>0,
            "stale cache was not rejected");
}
}

int main(){try{lifecycleAndStaleGate();std::cout<<"PASS phase3f_mp_subspace\n";return 0;}
catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
