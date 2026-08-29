#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gaussian_trial_state.h"
#include "pure_projector_fast.h"
#include "pure_projector_stack.h"

namespace {

void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}

double relativeError(const MatType &a, const MatType &b) {
    return (a - b).norm() / std::max(1.0, b.norm());
}

MatType canonicalPhi(int modes) {
    MatType phi = MatType::Zero(2 * modes, modes);
    for (int mode = 0; mode < modes; ++mode) {
        phi(2 * mode, mode) = 1.0 / std::sqrt(2.0);
        phi(2 * mode + 1, mode) = DataType(0.0, 1.0) / std::sqrt(2.0);
    }
    return phi;
}

MatType planeRotation(int dimension, int a, int b, double angle) {
    MatType factor = MatType::Identity(dimension, dimension);
    factor(a,a)=std::cos(angle); factor(b,b)=std::cos(angle);
    factor(a,b)=std::sin(angle); factor(b,a)=-std::sin(angle);
    return factor;
}

PureProjectorSlice slice(int dimension, int a, int b, double angle,
                         const std::string &label) {
    return PureProjectorSlice(planeRotation(dimension,a,b,angle), 1.0, label);
}

std::vector<PureProjectorSlice> noncommutingSlices() {
    return {slice(6,0,2,0.08,"a"), slice(6,1,2,-0.06,"b"),
            slice(6,2,4,0.05,"c"), slice(6,0,3,-0.04,"d"),
            slice(6,3,5,0.07,"e")};
}

void cutRoundTripAndBlocks() {
    const GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(3));
    const auto slices = noncommutingSlices();
    PureProjectorStackManager manager(trial, slices, 2);
    const MatType initialRight = manager.phiRight();
    const MatType initialLeft = manager.phiLeft();
    const MatType initialGreen = manager.green().green;
    const std::uint64_t initialHash = manager.configurationHash();
    require(manager.forwardToEnd(), "forward sweep failed");
    require(manager.backwardToStart(), "backward sweep failed");
    require(manager.cut() == 0, "round trip did not return to left boundary");
    require(relativeError(manager.phiRight(), initialRight) < 1e-11,
            "right subspace changed after round trip");
    require(relativeError(manager.phiLeft(), initialLeft) < 1e-11,
            "left subspace changed after round trip");
    require(relativeError(manager.green().green, initialGreen) < 1e-11,
            "Green changed after cut round trip");
    require(manager.configurationHash() == initialHash,
            "configuration changed during cut movement");

    MatType right = trial.Phi;
    for (int i=0;i<3;++i) right=slices[i].matrix*right;
    MatType left = trial.Phi;
    for (int i=int(slices.size())-1;i>=3;--i) left=slices[i].matrix.adjoint()*left;
    const PureProjectorGreenResult expected = pureProjectorGreenThinQr(right,left);
    require(expected.ok(), "direct center rebuild failed");
    for (int block : {1,2,4}) {
        PureProjectorStackManager blocked(trial,slices,block);
        require(blocked.moveToCut(3), "blocked move failed");
        require(relativeError(blocked.green().green,expected.green)<1e-11,
                "block-size Green mismatch");
    }
    std::cout << "PASS cut_round_trip\nPASS stable_rebuild_blocks\n";
}

void orientationTest() {
    const GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(3));
    const auto slices = noncommutingSlices();
    PureProjectorStackManager manager(trial,slices,2);
    require(manager.moveToCut(2), "move to noncommuting cut failed");
    MatType right=slices[1].matrix*slices[0].matrix*trial.Phi;
    MatType left=trial.Phi;
    for(int i=int(slices.size())-1;i>=2;--i)left=slices[i].matrix.adjoint()*left;
    PureProjectorGreenResult correct=pureProjectorGreenThinQr(right,left);
    require(correct.ok() && relativeError(manager.green().green,correct.green)<1e-11,
            "left/right orientation mismatch");
    MatType wrong=slices[0].matrix*slices[1].matrix*trial.Phi;
    PureProjectorGreenResult wrongGreen=pureProjectorGreenThinQr(wrong,left);
    require(wrongGreen.ok() && relativeError(manager.green().green,wrongGreen.green)>1e-5,
            "orientation test is insensitive to wrong factor order");
    std::cout << "PASS left_right_orientation\n";
}

void proposalTransactionTest() {
    const GaussianTrialState trial = GaussianTrialState::fromPhi(canonicalPhi(3));
    PureFastConfiguration configuration;
    configuration.slices=noncommutingSlices();
    configuration.hs_fields={1,1,1,1,1};
    configuration.locations.resize(configuration.slices.size());
    for(int i=0;i<int(configuration.locations.size());++i) {
        configuration.locations[i].branch=i<3?PureBranch::Ket:PureBranch::Bra;
        configuration.locations[i].slice=i/3;
        configuration.locations[i].factor=i;
        configuration.locations[i].bond=i%2;
        configuration.locations[i].aux=0;
    }
    PureProjectorFastWalker walker(trial,configuration,2,PureFastRunMode::AuditLockstep);
    const std::uint64_t before=walker.configurationHash();
    const MatType oldGreen=walker.measurementGreen().green;
    PureFastProposal proposal;
    proposal.index=1; proposal.new_hs=-1; proposal.uniform=0.37;
    proposal.new_factor=configuration.slices[1].matrix.adjoint();
    PureFastProposalResult result=walker.propose(proposal);
    require(result.snapshot.uniform==proposal.uniform,"proposal uniform changed");
    require(result.snapshot.configuration_hash==before,"mutation occurred before decision");
    require(result.ratio.ok(),"trusted single proposal failed");
    require(result.ratio.relative_reference_error<1e-10,"fast/slow ratio mismatch");
    if(!result.accepted) {
        require(walker.configurationHash()==before,"rejected proposal mutated HS");
        require(relativeError(walker.measurementGreen().green,oldGreen)<1e-12,
                "rejected proposal mutated Green");
    }
    std::cout << "PASS proposal_transaction\n";
}

} // namespace

int main() {
    try {
        cutRoundTripAndBlocks();
        orientationTest();
        proposalTransactionTest();
        return 0;
    } catch(const std::exception &error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
