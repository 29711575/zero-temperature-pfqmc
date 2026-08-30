#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pure_projector_fast.h"

namespace {

struct Task {int id;double V;std::uint64_t seed;};
Task taskFor(int id){switch(id){case 16:return{16,4,712041};case 17:return{17,4,712042};
    case 20:return{20,6,712061};case 21:return{21,6,712062};default:throw std::invalid_argument("task must be 16,17,20,21");}}
MatType kinetic(double mu){auto matrix=[&](double d,double m){SpinlessTvChainUtils c(12,1,0,2,0,d,m,0);
    MatType h=MatType::Zero(24,24);c.KineticGenerator(h);return h;};MatType base=matrix(0,0);
    return base+(matrix(1,0)-base)+mu*(matrix(0,1)-base);}
MatType exponential(MatType h,double scale){return expm(h,scale);}
MatType hsFactor(const SpinlessTvChainUtils&c,int layer,int aux,int sigma){MatType h=MatType::Zero(24,24);
    const double lambda=std::acosh(std::exp(.5*c.V*c.dt));int a,b,d,e;
    c.aux2MajoranaIdx(aux,0,layer,a,b);c.aux2MajoranaIdx(aux,1,layer,d,e);
    const DataType z(0,lambda*sigma);h(a,b)=z;h(b,a)=-z;h(d,e)=z;h(e,d)=-z;
    return exponential(h,1.0);}
using HsCache=std::array<std::array<std::array<MatType,2>,6>,2>;
HsCache makeHsCache(const SpinlessTvChainUtils&model){HsCache cache;
    for(int layer=0;layer<2;++layer)for(int aux=0;aux<6;++aux){
        cache[layer][aux][0]=hsFactor(model,layer,aux,-1);
        cache[layer][aux][1]=hsFactor(model,layer,aux,1);}return cache;}
PureFastConfiguration ketProtocol(const SpinlessTvChainUtils&model,const HsCache&cache,
                                  std::uint64_t seed){
    std::mt19937_64 rng(seed);std::uniform_int_distribution<int>bit(0,1);
    const MatType half=exponential(kinetic(0),-.05);PureFastConfiguration out;int ordinal=0;
    auto push=[&](const MatType&m,int field,int slice,int bond,int aux,const std::string&label){
        out.slices.emplace_back(m,1.0,label);out.hs_fields.push_back(field);
        out.locations.push_back({PureBranch::Ket,slice,ordinal++,bond,aux});};
    for(int slice=0;slice<120;++slice){push(half,0,slice,-1,-1,"K/2");
        for(int layer=0;layer<2;++layer)for(int aux=0;aux<6;++aux){int sigma=bit(rng)?1:-1;
            push(cache[layer][aux][sigma>0],sigma,slice,layer,aux,"V");}
        push(half,0,slice,-1,-1,"K/2");}return out;}
std::vector<int> centerBranchProposalIndices(const PureFastConfiguration&c){
    std::vector<int> ket,bra,result;
    for(int i=0;i<int(c.locations.size());++i)if(c.locations[i].aux>=0)
        (c.locations[i].branch==PureBranch::Ket?ket:bra).push_back(i);
    if(ket.empty()||bra.empty())throw std::runtime_error("missing branch proposal indices");
    for(int step=0;step<32;++step){result.push_back(ket[ket.size()-1-step%ket.size()]);
        result.push_back(bra[step%bra.size()]);}return result;}
PureFastProposal flip(const HsCache&cache,const PureFastConfiguration&c,int index,double u){
    const auto&loc=c.locations[index];PureFastProposal p;p.index=index;p.new_hs=-c.hs_fields[index];
    p.new_factor=cache[loc.bond][loc.aux][p.new_hs>0];p.new_eta=1.0;p.uniform=u;return p;}
int mirrorMismatches(const PureFastConfiguration&c){const int half=int(c.slices.size()/2);int mismatch=0;
    for(int i=0;i<half;++i)if(c.locations[i].aux>=0&&c.hs_fields[i]!=c.hs_fields[2*half-1-i])++mismatch;
    return mismatch;}

} // namespace

int main(int argc,char**argv){try{if(argc!=3&&argc!=5)
        throw std::invalid_argument("usage: task output.csv [begin end]");
    const int begin=argc==5?std::stoi(argv[3]):0,end=argc==5?std::stoi(argv[4]):128;
    if(begin<0||end>128||begin>=end)throw std::invalid_argument("invalid attempt range");
    const Task task=taskFor(std::stoi(argv[1]));SpinlessTvChainUtils model(12,.1,task.V,2,0,1,0,0);
    const HsCache cache=makeHsCache(model);
    GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(kinetic(.3));
    if(trial.fermionParity()!=-1)throw std::runtime_error("L12 trial parity mismatch");
    std::ofstream out(argv[2]);if(!out)throw std::runtime_error("CSV open failed");
    out<<"task,L,V,seed,attempt,initialization_policy,initialization_success,initial_z2,center_overlap_rcond,center_solve_residual,center_green_residual,proposals,accepted,independent_field_mismatches,mp_fallback_count,mp_fallback_failures,endpoint_rebuild_residual_max,final_hash,status\n";
    for(int attempt=begin;attempt<end;++attempt){const std::uint64_t ketSeed=task.seed+0x9e3779b97f4a7c15ULL*std::uint64_t(attempt);
        auto mirrored=pureProjectorMirroredConfiguration(ketProtocol(model,cache,ketSeed),1e-12);
        if(!mirrored.ok())throw std::runtime_error("mirrored construction failed");const int center=int(mirrored.configuration.slices.size()/2);
        PureEndpointRebuildResult endpoint=pureProjectorEndpointRebuild(trial,mirrored.configuration.slices,center,8);
        if(!endpoint.ok()||endpoint.overlap_rcond<1-1e-10||endpoint.green_residual>1e-10)
            throw std::runtime_error("mirrored endpoint initialization failed");
        PureFastOptions options;PureProjectorFastWalker walker(trial,std::move(mirrored.configuration),8,
            PureFastRunMode::FastStrict,options,PureFastInitializationPolicy::MirroredTheoremZ2Plus);
        std::mt19937_64 proposalRng(ketSeed^0xd1b54a32d192ed03ULL);std::uniform_real_distribution<double>uniform(0,1);
        auto indices=centerBranchProposalIndices(walker.configuration());int accepted=0;
        for(int step=0;step<64;++step){const int index=indices[step%indices.size()];
            auto result=walker.propose(flip(cache,walker.configuration(),index,uniform(proposalRng)));
            if(result.terminated||!result.ratio.ok())throw std::runtime_error("ordinary independent proposal failed");
            accepted+=result.accepted;}
        const auto&d=walker.diagnostics();const int mismatch=mirrorMismatches(walker.configuration());
        if(mismatch==0)throw std::runtime_error("bra/ket fields remained mirror constrained");
        out<<task.id<<",12,"<<task.V<<','<<task.seed<<','<<attempt
           <<",mirrored_theorem_z2_plus,1,1,"<<std::setprecision(17)<<endpoint.overlap_rcond<<','
           <<endpoint.solve_residual<<','<<endpoint.green_residual<<",64,"<<accepted<<','<<mismatch<<','
           <<d.mp_fallback_count<<','<<d.mp_fallback_failure_count<<','
           <<d.maximum_endpoint_rebuild_green_residual<<','<<walker.configurationHash()<<",PASS\n";
        out.flush();if(!out)throw std::runtime_error("CSV row flush failed");}
    out.flush();if(!out)throw std::runtime_error("CSV flush failed");out.close();if(out.fail())
        throw std::runtime_error("CSV close failed");std::cout<<"PASS mirrored_task_"<<task.id
            <<' '<<begin<<':'<<end<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
