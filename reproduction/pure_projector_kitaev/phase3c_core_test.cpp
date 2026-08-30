#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "gaussian_trial_state.h"
#include "pure_projector_endpoint.h"
#include "pure_projector_fast.h"

namespace {

void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}

MatType canonicalPhi(int modes) {
    MatType phi = MatType::Zero(2 * modes, modes);
    for (int mode = 0; mode < modes; ++mode) {
        phi(2 * mode, mode) = 1.0 / std::sqrt(2.0);
        phi(2 * mode + 1, mode) = DataType(0.0, 1.0) / std::sqrt(2.0);
    }
    return phi;
}

MatType rotation(int dimension, int a, int b, double angle) {
    MatType result = MatType::Identity(dimension, dimension);
    result(a,a)=std::cos(angle); result(b,b)=std::cos(angle);
    result(a,b)=std::sin(angle); result(b,a)=-std::sin(angle);
    return result;
}

MatType boostFactor(double angle) {
    MatType result=MatType::Identity(2,2);
    result(0,0)=std::cosh(angle);result(1,1)=std::cosh(angle);
    result(0,1)=DataType(0.0,std::sinh(angle));
    result(1,0)=DataType(0.0,-std::sinh(angle));
    return result;
}

void endpointDoesNotInspectPrefixes() {
    const GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(2));
    // A deliberately singular intermediate trial overlap followed by its exact
    // inverse.  The final contour is the identity and is therefore trustworthy.
    MatType crack = MatType::Identity(4,4);
    crack.col(0).swap(crack.col(1));
    crack.col(2).swap(crack.col(3));
    std::vector<PureProjectorSlice> slices{
        {crack,1.0,"prefix_crack"},{crack.adjoint(),1.0,"repair"}};
    PureFastOptions fast;
    const PureProjectorWeightResult sequential =
        pureProjectorStableReferenceWeight(trial,slices,fast);
    require(!sequential.ok(), "fixture did not trigger the sequential prefix gate");
    PureEndpointRebuildResult endpoint = pureProjectorEndpointRebuild(
        trial,slices,2,1,PureProjectorOptions());
    require(endpoint.ok(), "endpoint-only rebuild rejected a trustworthy endpoint");
    require(endpoint.green_residual < 1e-12, "endpoint Green residual is too large");
    require((endpoint.green-trial.G_T).norm() < 1e-12,
            "identity endpoint Green changed");
    std::cout << "PASS endpoint_ignores_prefix_green\n";
}

void mirroredTheoremInitializer() {
    const GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(3));
    PureFastConfiguration ket;
    for (int i=0;i<5;++i) {
        ket.slices.emplace_back(rotation(6,i%4,(i+2)%6,0.03*(i+1)),
                                DataType(1.0,0.0),"f"+std::to_string(i));
        ket.hs_fields.push_back(i%2?1:-1);
        ket.locations.push_back({PureBranch::Ket,i/2,i,i%2,i});
    }
    PureMirroredInitializationResult mirrored =
        pureProjectorMirroredConfiguration(ket,1e-12);
    require(mirrored.ok(), "theorem-backed mirrored construction failed");
    require(mirrored.initial_z2_sign==1, "mirrored initializer did not set Z2=+1");
    require(mirrored.initialization_policy=="mirrored_theorem_z2_plus",
            "mirrored policy provenance is missing");
    const int cut=int(ket.slices.size());
    PureEndpointRebuildResult endpoint=pureProjectorEndpointRebuild(
        trial,mirrored.configuration.slices,cut,2,PureProjectorOptions());
    require(endpoint.ok(), "mirrored center endpoint rebuild failed");
    require(endpoint.overlap_rcond>1.0-1e-12,
            "mirrored center overlap is not unit conditioned");
    require(endpoint.green_residual<1e-12,"mirrored center Green residual is too large");

    PureFastOptions options;
    PureProjectorFastWalker walker(trial,mirrored.configuration,2,
        PureFastRunMode::FastStrict,options,
        PureFastInitializationPolicy::MirroredTheoremZ2Plus);
    require(walker.z2Sign()==1,"production walker lost theorem-backed initial Z2");
    require(walker.initializationPolicy()=="mirrored_theorem_z2_plus",
            "walker policy provenance mismatch");
    std::cout << "PASS mirrored_theorem_initializer\n";
}

void scalarPrefactorMustBePositiveReal() {
    PureFastConfiguration ket;
    ket.slices.emplace_back(MatType::Identity(2,2),DataType(0.0,0.0),"zero_eta");
    ket.hs_fields.push_back(1);
    ket.locations.push_back({PureBranch::Ket,0,0,0,0});
    PureMirroredInitializationResult result=pureProjectorMirroredConfiguration(ket,1e-12);
    require(!result.ok(),"zero scalar prefactor did not fail closed");
    std::cout << "PASS mirrored_scalar_prefactor_gate\n";
}

PureFastConfiguration oneModeMirrored(double oldAngle) {
    PureFastConfiguration ket;
    ket.slices.emplace_back(boostFactor(oldAngle),1.0,"boost");
    ket.hs_fields.push_back(1);ket.locations.push_back({PureBranch::Ket,0,0,0,0});
    auto mirrored=pureProjectorMirroredConfiguration(ket,1e-12);
    require(mirrored.ok(),"one-mode mirrored fixture failed");
    return mirrored.configuration;
}

void mpSameProposalTransactions() {
    const GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(1));
    PureFastOptions options;options.decision_margin_tolerance=1.0;
    // Pure-Majorana overlap obeys W^2=det(Phi^dagger B Phi), hence changing
    // one of the two mirrored factors gives the square-root ratio exp(-0.2).
    const DataType expected(std::exp(-0.2),0.0);
    for (const auto &entry:std::vector<std::pair<double,bool>>{{0.5,true},{0.9,false}}) {
        PureFastConfiguration configuration=oneModeMirrored(-0.2);
        PureProjectorFastWalker walker(trial,configuration,1,PureFastRunMode::FastStrict,
            options,PureFastInitializationPolicy::MirroredTheoremZ2Plus);
        const std::uint64_t before=walker.configurationHash();
        PureFastProposal proposal;proposal.index=0;proposal.new_hs=-1;
        proposal.new_factor=boostFactor(0.2);proposal.new_eta=1.0;proposal.uniform=entry.first;
        PureFastProposalResult result=walker.propose(proposal);
        require(result.ratio.mp_reference.ok(),"MP160/320 proposal did not converge");
        require(result.ratio.mp_reference.precision_digits==320,
                "simple proposal unexpectedly required MP640");
        if(!(std::abs(result.ratio.ratio-expected)<1e-12)){
            std::cerr<<"MP ratio="<<result.ratio.ratio<<" expected="<<expected
                     <<" z2="<<result.ratio.mp_reference.ratio_z2
                     <<" reality="<<result.ratio.mp_reference.reality_error<<'\n';
            throw std::runtime_error("MP signed ratio mismatch");
        }
        require(result.accepted==entry.second,"same-uniform MP decision mismatch");
        require(result.snapshot.uniform==entry.first,"MP fallback consumed/replaced uniform");
        require(result.snapshot.configuration_hash==before,"MP fallback mutated before decision");
        require(!result.terminated&&result.ratio.ok(),"MP transaction failed closed");
        require((walker.configurationHash()!=before)==entry.second,
                "accepted/rejected MP mutation ordering mismatch");
        require(walker.z2Sign()==1,"positive MP ratio changed Z2");
        require(walker.measurementGreen().ok(),"endpoint stack recovery failed after MP decision");
    }

    PureFastOptions declined=options;declined.mp_same_proposal_fallback=false;
    PureProjectorFastWalker walker(trial,oneModeMirrored(-0.2),1,
        PureFastRunMode::FastStrict,declined,
        PureFastInitializationPolicy::MirroredTheoremZ2Plus);
    PureFastProposal proposal;proposal.index=0;proposal.new_hs=-1;
    proposal.new_factor=boostFactor(0.2);proposal.new_eta=1.0;proposal.uniform=0.5;
    PureFastProposalResult fatal=walker.propose(proposal);
    require(fatal.terminated&&fatal.ratio.status==PureFastRatioStatus::reference_failure,
            "declined MP fallback was not fatal");
    std::cout << "PASS mp_same_proposal_accept\nPASS mp_same_proposal_reject\n"
              << "PASS mp_fallback_declined_fatal\n";
}

void observerDoesNotChangeTrajectory() {
    const GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(1));
    PureFastOptions off,on;on.read_only_audit_interval=1;
    PureProjectorFastWalker a(trial,oneModeMirrored(-0.2),1,PureFastRunMode::FastStrict,
        off,PureFastInitializationPolicy::MirroredTheoremZ2Plus);
    PureProjectorFastWalker b(trial,oneModeMirrored(-0.2),1,PureFastRunMode::FastStrict,
        on,PureFastInitializationPolicy::MirroredTheoremZ2Plus);
    for(int step=0;step<8;++step){PureFastProposal pa;pa.index=step%2;pa.new_hs=-a.configuration().hs_fields[pa.index];
        pa.new_factor=a.configuration().slices[pa.index].matrix.adjoint().eval();pa.new_eta=1.0;
        pa.uniform=0.13+0.09*step;PureFastProposal pb=pa;
        pb.new_hs=-b.configuration().hs_fields[pb.index];
        pb.new_factor=b.configuration().slices[pb.index].matrix.adjoint().eval();
        auto ra=a.propose(pa),rb=b.propose(pb);require(!ra.terminated&&!rb.terminated,
            "observer trajectory terminated");require(ra.accepted==rb.accepted,
            "observer changed accept/reject decision");require(a.configurationHash()==b.configurationHash(),
            "observer changed HS trajectory hash");require(a.z2Sign()==b.z2Sign(),
            "observer changed Z2 trajectory");}
    require(b.diagnostics().ratio_slow_reference_count==8,
        "read-only endpoint observer did not run");
    std::cout<<"PASS observer_trajectory_hash_invariant\n";
}

void malformedMpReferencesAreFatal() {
    const GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(1));
    for(int kind=0;kind<2;++kind){PureFastOptions options;options.decision_margin_tolerance=1.0;
        PureProjectorFastWalker walker(trial,oneModeMirrored(-0.2),1,
            PureFastRunMode::FastStrict,options,
            PureFastInitializationPolicy::MirroredTheoremZ2Plus);
        PureFastProposal proposal;proposal.index=0;proposal.new_hs=-1;proposal.new_eta=1.0;
        proposal.uniform=0.4;proposal.new_factor=kind==0?MatType::Zero(2,2):boostFactor(0.2);
        if(kind==1)proposal.new_factor(0,0)=DataType(
            std::numeric_limits<double>::quiet_NaN(),0.0);
        const auto result=walker.propose(proposal);
        require(result.terminated&&result.ratio.status==PureFastRatioStatus::reference_failure,
            "malformed/nonfinite MP reference was not fatal");}
    std::cout<<"PASS mp_malformed_fatal\nPASS mp_nonfinite_fatal\n";
}

} // namespace

int main() {
    try {
        endpointDoesNotInspectPrefixes();
        mirroredTheoremInitializer();
        scalarPrefactorMustBePositiveReal();
        mpSameProposalTransactions();
        observerDoesNotChangeTrajectory();
        malformedMpReferencesAreFatal();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
