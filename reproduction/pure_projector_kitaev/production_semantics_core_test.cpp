#include <cmath>
#include <iostream>
#include <stdexcept>

#include "gaussian_trial_state.h"
#include "pure_projector_fast.h"
#include "pure_projector_observables.h"

namespace {
void require(bool value,const char*message){if(!value)throw std::runtime_error(message);}

MatType canonicalPhi(int modes){MatType p=MatType::Zero(2*modes,modes);
    for(int i=0;i<modes;++i){p(i,i)=1/std::sqrt(2.0);p(modes+i,i)=DataType(0,1/std::sqrt(2.0));}
    return p;}

void normalization(){for(int L:{2,4,6}){GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(L));
    DataType value=pureProjectorStructureFactor(trial.G_T,L,0.0);
    DataType raw=pureProjectorStructureFactorUnnormalized(trial.G_T,L,0.0);
    require(std::abs(value-raw/double(L*L))<1e-14,"S(q) is not sum/L^2");}}

void physicalParity(){for(int L:{2,4,6}){GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(L));
    PurePhysicalParityResult p=pureProjectorPhysicalParity(trial.G_T);
    require(p.ok(),"physical parity failed");
    int reorder=((L*(L-1)/2)&1)?-1:1;
    require(p.physical_parity==reorder*p.internal_pfaffian_sign,
            "block-order parity factor missing");
    require(trial.fermionParity()==p.physical_parity,"trial parity is not physical parity");}}

void bondCounts(){for(int L:{2,4,6}){
    require(pureProjectorCheckerboardBondCounts(L,0)==std::make_pair(L/2,L/2),
            "PBC checkerboard bond counts are wrong");
    require(pureProjectorCheckerboardBondCounts(L,1)==std::make_pair(L/2,(L-1)/2),
            "OBC checkerboard bond counts are wrong");}}

void greenStructureAndRankZero(){GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(2));
    auto g=pureProjectorGreenThinQr(trial.Phi,trial.Phi);require(g.ok(),"identity Green failed");
    require(g.green_skew_residual<1e-12&&g.green_diagonal_residual<1e-12,
            "Green skew/diagonal diagnostics missing");
    MatType invalidRight=trial.Phi;invalidRight.row(0)*=1.3;
    PureProjectorOptions strict;strict.fail_on_green_structure=true;
    auto invalid=pureProjectorGreenThinQr(invalidRight,trial.Phi,strict);
    require(invalid.status==PureProjectorStatus::green_structure_exceeded,
            "non-Gaussian propagated subspace did not fail the Green structure gate");
    std::vector<PureProjectorSlice>s{{MatType::Identity(4,4),1.0,"identity"}};
    PureFastOptions options;auto ratio=pureFastLocalRatio(trial,s,0,MatType::Identity(4,4),1.0,options);
    require(ratio.ok()&&ratio.low_rank==0&&ratio.fast_updated_green.size()!=0,
            "rank-zero Green status was not checked");}
}

int main(){try{normalization();physicalParity();bondCounts();greenStructureAndRankZero();
    std::cout<<"PASS production_semantics_core\n";return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
