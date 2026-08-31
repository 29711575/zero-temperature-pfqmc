#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "gaussian_trial_state.h"
#include "pure_projector_mp.h"

namespace {
void require(bool value,const std::string&message){if(!value)throw std::runtime_error(message);}
MatType canonicalPhi(int modes){MatType p=MatType::Zero(2*modes,modes);
    for(int i=0;i<modes;++i){p(i,i)=1/std::sqrt(2.0);
        p(modes+i,i)=DataType(0,1/std::sqrt(2.0));}return p;}
MatType factor(int dimension,int a,int b,double angle){MatType m=MatType::Identity(dimension,dimension);
    m(a,a)=std::cosh(angle);m(b,b)=std::cosh(angle);
    m(a,b)=DataType(0,std::sinh(angle));m(b,a)=DataType(0,-std::sinh(angle));return m;}
}

int main(int argc,char**argv){try{
    require(argc==3,"usage: phase3f_oracle_equivalence oracle.csv invalidation.csv");
    const int modes=3,dimension=6,count=128,midpoint=count/2;
    const GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(modes));
    std::vector<int>sign(count,1);std::vector<PureProjectorSlice>slices;
    for(int i=0;i<count;++i){sign[i]=(i%3)?1:-1;slices.emplace_back(
        factor(dimension,i%modes,modes+i%modes,sign[i]*0.004),1.0,"canonical_order");}
    PureMpProposalOptions optimized,phase3e,fullLegacy;
    optimized.enable_subspace_cache=true;phase3e.enable_subspace_cache=false;
    fullLegacy.enable_operator_cache=false;fullLegacy.enable_subspace_cache=false;
    PureMpSubspaceCache cache;std::mt19937_64 rng(930601);
    std::uniform_real_distribution<double>uniform(0,1);std::ofstream oracle(argv[1]),events(argv[2]);
    require(bool(oracle)&&bool(events),"output open failed");
    oracle<<"proposal,branch,index,block_boundary,stale_injected,uniform,new_ratio_real,new_ratio_imag,phase3e_ratio_real,phase3e_ratio_imag,legacy_ratio_real,legacy_ratio_imag,ratio_error_phase3e,ratio_error_legacy,z2_new,z2_phase3e,z2_legacy,status_match,precision_match,decision_new,decision_phase3e,decision_legacy,pre_residual_new,pre_residual_legacy,post_residual_new,post_residual_legacy,restore_count,saved_factors,legacy_fallbacks,status\n";
    events<<"proposal,event,index,old_hash,new_hash,accepted,invalidated_entries,restore_count_total,stale_rejections_total,checkpoint_entries_peak,status\n";
    std::uint64_t hash=pureProjectorMpConfigurationHash(slices);double maximumError=0,maximumResidualDelta=0;
    int mismatch=0,acceptedCount=0,rejectedCount=0;
    for(int proposal=0;proposal<512;++proposal){int index;
        if(proposal%4==0)index=(proposal*7)%midpoint;
        else if(proposal%4==1)index=midpoint+(proposal*11)%midpoint;
        else if(proposal%4==2)index=((proposal/4)*8)%count;
        else index=int(rng()%count);
        bool stale=false;if(proposal==257){const int staleIndex=5;sign[staleIndex]*=-1;
            slices[staleIndex].matrix=factor(dimension,staleIndex%modes,
                modes+staleIndex%modes,sign[staleIndex]*0.004);
            hash=pureProjectorMpConfigurationHash(slices);stale=true;}
        const MatType candidate=factor(dimension,index%modes,modes+index%modes,-sign[index]*0.004);
        // Deterministically exercise rejected-cache semantics as well as ordinary accepted moves.
        // The scalar-prefactor change is part of the same frozen proposal and is sent identically
        // to all three oracles; it is never committed because its Metropolis ratio is below 0.9.
        const bool forceReject=(proposal%31)==0;
        const DataType candidateEta=forceReject?DataType(0.1,0):DataType(1,0);
        const double u=forceReject?0.9:uniform(rng);const auto current=pureProjectorMpSameProposal(
            trial,slices,index,candidate,candidateEta,optimized,&cache,hash);
        const auto reference=pureProjectorMpSameProposal(
            trial,slices,index,candidate,candidateEta,phase3e,nullptr,hash);
        const auto legacy=pureProjectorMpSameProposal(
            trial,slices,index,candidate,candidateEta,fullLegacy,nullptr,hash);
        require(current.ok()&&reference.ok()&&legacy.ok(),"oracle returned untrusted result");
        const double errorPhase3e=std::abs(current.ratio-reference.ratio)/
            std::max(1e-14,std::abs(reference.ratio));
        const double errorLegacy=std::abs(current.ratio-legacy.ratio)/
            std::max(1e-14,std::abs(legacy.ratio));
        maximumError=std::max({maximumError,errorPhase3e,errorLegacy});
        const double residualDelta=std::max(
            std::abs(current.pre_endpoint_residual-legacy.pre_endpoint_residual),
            std::abs(current.post_endpoint_residual-legacy.post_endpoint_residual));
        maximumResidualDelta=std::max(maximumResidualDelta,residualDelta);
        const bool decision=u<std::min(1.0,std::abs(current.ratio));
        const bool phase3eDecision=u<std::min(1.0,std::abs(reference.ratio));
        const bool legacyDecision=u<std::min(1.0,std::abs(legacy.ratio));
        const bool statusMatch=current.status==reference.status&&current.status==legacy.status;
        const bool precisionMatch=current.precision_digits==reference.precision_digits&&
            current.precision_digits==legacy.precision_digits;
        const bool equal=errorPhase3e<=1e-12&&errorLegacy<=1e-12&&
            current.ratio_z2==reference.ratio_z2&&current.ratio_z2==legacy.ratio_z2&&
            statusMatch&&precisionMatch&&decision==phase3eDecision&&decision==legacyDecision&&
            residualDelta<=1e-12;if(!equal)++mismatch;
        oracle<<proposal<<','<<(index<midpoint?"ket":"bra")<<','<<index<<','<<(index%8==0)<<','<<stale<<','
            <<std::setprecision(17)<<u<<','<<current.ratio.real()<<','<<current.ratio.imag()<<','
            <<reference.ratio.real()<<','<<reference.ratio.imag()<<','<<legacy.ratio.real()<<','
            <<legacy.ratio.imag()<<','<<errorPhase3e<<','<<errorLegacy<<','<<current.ratio_z2<<','
            <<reference.ratio_z2<<','<<legacy.ratio_z2<<','<<statusMatch<<','<<precisionMatch<<','
            <<decision<<','<<phase3eDecision<<','<<legacyDecision<<','<<current.pre_endpoint_residual<<','
            <<legacy.pre_endpoint_residual<<','<<current.post_endpoint_residual<<','
            <<legacy.post_endpoint_residual<<','<<current.profile.subspace_restore_count<<','
            <<current.profile.subspace_saved_factor_count<<','
            <<current.profile.subspace_legacy_fallback_count<<','<<(equal?"PASS":"FAIL")<<'\n';
        const std::uint64_t oldHash=hash,beforeInvalidations=cache.profile().subspace_cache_invalidations;
        if(decision){++acceptedCount;sign[index]*=-1;slices[index].matrix=candidate;
            hash=pureProjectorMpConfigurationHash(slices);const auto delta=cache.acceptedUpdate(
                oldHash,hash,index,slices);events<<proposal<<",accepted,"<<index<<','<<oldHash<<','<<hash
                <<",1,"<<delta.subspace_cache_invalidations<<','<<cache.profile().subspace_restore_count
                <<','<<cache.profile().subspace_stale_rejection_count<<','
                <<cache.profile().subspace_checkpoint_entries_peak<<",PASS\n";}
        else{++rejectedCount;events<<proposal<<",rejected,"<<index<<','<<oldHash<<','<<hash
                <<",0,"<<(cache.profile().subspace_cache_invalidations-beforeInvalidations)<<','
                <<cache.profile().subspace_restore_count<<','<<cache.profile().subspace_stale_rejection_count
                <<','<<cache.profile().subspace_checkpoint_entries_peak<<",PASS\n";}
    }
    oracle.flush();events.flush();require(bool(oracle)&&bool(events),"output write failed");
    oracle.close();events.close();require(!oracle.fail()&&!events.fail(),"output close failed");
    require(mismatch==0,"new and legacy MP oracles disagreed");
    require(acceptedCount>0&&rejectedCount>0,"accepted/rejected cache paths not both covered");
    require(cache.profile().subspace_restore_count>0,"production checkpoint restore never occurred");
    require(cache.profile().subspace_stale_rejection_count>0,"stale cache test did not fire");
    std::cout<<"{\"status\":\"PASS\",\"proposals\":512,\"accepted\":"<<acceptedCount
        <<",\"rejected\":"<<rejectedCount<<",\"mismatch\":"<<mismatch
        <<",\"maximum_ratio_error\":"<<std::setprecision(17)<<maximumError
        <<",\"maximum_residual_delta\":"<<maximumResidualDelta
        <<",\"restore_count\":"<<cache.profile().subspace_restore_count
        <<",\"saved_factors\":"<<cache.profile().subspace_saved_factor_count
        <<",\"stale_rejections\":"<<cache.profile().subspace_stale_rejection_count<<"}\n";
    return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
