#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pure_projector_fast.h"

namespace {

struct Task {int id,L;double V;std::uint64_t seed;int parity,attempted,index;
    double uniform,ratio;std::uint64_t pre_hash;bool accepted;};
const std::vector<Task> tasks={
    {4,6,4.0,706041,1,110,146,0.52928038756354223,
     0.707419140654256205736,7807527216050905631ULL,true},
    {8,6,6.0,706061,1,415,553,0.56577673777691728,
     0.002657388953068520511,5178673164861548103ULL,false}};

MatType kineticGenerator(int L,double delta,double mu){
    auto matrix=[&](double d,double m){SpinlessTvChainUtils c(L,1,0,2,0,d,m,0);
        MatType h=MatType::Zero(2*L,2*L);c.KineticGenerator(h);return h;};
    MatType base=matrix(0,0);return base+delta*(matrix(1,0)-base)+mu*(matrix(0,1)-base);
}
MatType exponential(MatType generator,double scale){return expm(generator,scale);}
MatType localHsFactor(const SpinlessTvChainUtils&c,int bond,int aux,int sigma){
    MatType generator=MatType::Zero(c.nDim,c.nDim);
    const double lambda=std::acosh(std::exp(.5*c.V*c.dt));int a,b,d,e;
    c.aux2MajoranaIdx(aux,0,bond,a,b);c.aux2MajoranaIdx(aux,1,bond,d,e);
    const DataType z(0,lambda*sigma);generator(a,b)=z;generator(b,a)=-z;
    generator(d,e)=z;generator(e,d)=-z;return exponential(generator,1.0);
}
GaussianTrialState trialFor(const Task&t){
    GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(
        kineticGenerator(t.L,1.0,0.3));
    if(trial.fermionParity()!=t.parity)throw std::runtime_error("trial parity mismatch");
    return trial;
}
PureFastConfiguration independentContour(const Task&t,const SpinlessTvChainUtils&model,
                                         std::mt19937_64&rng){
    const MatType half=exponential(kineticGenerator(t.L,1.0,0.0),-.05);
    const int perLayer=t.L/2,timeSlices=10*t.L;std::uniform_int_distribution<int>bit(0,1);
    auto branch=[&](PureBranch side){PureFastConfiguration out;int ordinal=0;
        auto push=[&](const MatType&m,int field,int slice,int bond,int aux,const std::string&label){
            out.slices.emplace_back(m,1.0,label);out.hs_fields.push_back(field);
            out.locations.push_back({side,slice,ordinal++,bond,aux});};
        for(int slice=0;slice<timeSlices;++slice){push(half,0,slice,-1,-1,"K/2");
            for(int layer=0;layer<2;++layer)for(int aux=0;aux<perLayer;++aux){
                const int sigma=bit(rng)?1:-1;push(localHsFactor(model,layer,aux,sigma),
                    sigma,slice,layer,aux,"V");}push(half,0,slice,-1,-1,"K/2");}return out;};
    PureFastConfiguration ket=branch(PureBranch::Ket),bra=branch(PureBranch::Bra),full=ket;
    for(int i=int(bra.slices.size())-1;i>=0;--i){full.slices.push_back(bra.slices[i]);
        full.hs_fields.push_back(bra.hs_fields[i]);full.locations.push_back(bra.locations[i]);}
    return full;
}
std::vector<int> proposalIndices(const PureFastConfiguration&c){std::vector<int> out;
    for(int i=0;i<int(c.locations.size());++i)if(c.locations[i].aux>=0)out.push_back(i);return out;}
PureFastProposal flip(const SpinlessTvChainUtils&model,const PureFastConfiguration&c,
                      int index,double uniform){const auto&loc=c.locations[index];PureFastProposal p;
    p.index=index;p.new_hs=-c.hs_fields[index];p.new_factor=localHsFactor(model,loc.bond,loc.aux,p.new_hs);
    p.new_eta=1.0;p.uniform=uniform;return p;}

struct Replay {PureFastProposalResult frozen;std::uint64_t post_hash=0;int continued=0;
    std::uint64_t mp_count=0,mp_failures=0;};
Replay replay(const Task&t){
    GaussianTrialState trial=trialFor(t);SpinlessTvChainUtils model(t.L,.1,t.V,2,0,1,0,0);
    std::mt19937_64 rng(t.seed);PureFastOptions options;options.read_only_audit_interval=0;
    std::unique_ptr<PureProjectorFastWalker> walker;
    while(!walker){PureFastConfiguration initial=independentContour(t,model,rng);
        try{walker.reset(new PureProjectorFastWalker(trial,std::move(initial),8,
            PureFastRunMode::FastStrict,options,PureFastInitializationPolicy::SequentialAudit));}
        catch(const std::runtime_error&){}}
    std::uniform_real_distribution<double> uniform(0,1);std::size_t cursor=0;int direction=1;
    Replay replay;
    const int total=t.attempted+64;
    for(int attempted=1;attempted<=total;++attempted){auto indices=proposalIndices(walker->configuration());
        const int index=indices[cursor];if(indices.size()>1){if(direction>0&&cursor+1==indices.size())direction=-1;
            else if(direction<0&&cursor==0)direction=1;cursor=std::size_t(int(cursor)+direction);}
        const double u=uniform(rng);const std::uint64_t pre=walker->configurationHash();
        PureFastProposalResult result=walker->propose(flip(model,walker->configuration(),index,u));
        if(attempted<t.attempted&&(result.terminated||!result.ratio.ok())){
            std::cerr<<"task="<<t.id<<" early_attempt="<<attempted<<" index="<<index
                <<" status="<<int(result.ratio.status)<<" mp_status="
                <<pureMpProposalStatusName(result.ratio.mp_reference.status)<<" message="
                <<result.ratio.mp_reference.message<<'\n';
            throw std::runtime_error("trajectory failed before frozen proposal");}
        if(attempted==t.attempted){
            if(index!=t.index||pre!=t.pre_hash||std::abs(u-t.uniform)>1e-15)
                throw std::runtime_error("frozen proposal identity/hash/uniform mismatch");
            if(result.terminated||!result.ratio.ok()||!result.ratio.mp_reference.ok()){
                std::cerr<<"task="<<t.id<<" ratio_status="<<int(result.ratio.status)
                    <<" mp_status="<<pureMpProposalStatusName(result.ratio.mp_reference.status)
                    <<" precision="<<result.ratio.mp_reference.precision_digits
                    <<" pre_rcond="<<result.ratio.mp_reference.pre_endpoint_rcond
                    <<" post_rcond="<<result.ratio.mp_reference.post_endpoint_rcond
                    <<" message="<<result.ratio.mp_reference.message<<'\n';
                throw std::runtime_error("frozen proposal MP fallback failed");}
            if(std::abs(result.ratio.ratio.real()-t.ratio)>1e-11||
               std::abs(result.ratio.ratio.imag())>1e-11||result.accepted!=t.accepted)
                throw std::runtime_error("frozen ratio/decision mismatch");
            replay.frozen=result;
        }else if(attempted>t.attempted){
            if(result.terminated||!result.ratio.ok()){
                std::cerr<<"task="<<t.id<<" continuation_attempt="<<attempted
                    <<" index="<<index<<" status="<<int(result.ratio.status)
                    <<" mp_status="<<pureMpProposalStatusName(result.ratio.mp_reference.status)
                    <<" pre_rcond="<<result.ratio.mp_reference.pre_endpoint_rcond
                    <<" post_rcond="<<result.ratio.mp_reference.post_endpoint_rcond
                    <<" message="<<result.ratio.mp_reference.message<<'\n';
                throw std::runtime_error("trajectory did not continue after frozen proposal");
            }
            ++replay.continued;
        }
    }
    replay.post_hash=walker->configurationHash();replay.mp_count=walker->diagnostics().mp_fallback_count;
    replay.mp_failures=walker->diagnostics().mp_fallback_failure_count;return replay;
}

} // namespace

int main(int argc,char**argv){try{const std::string path=argc>1?argv[1]:"mp_same_proposal_fallback.csv";
    std::ofstream out(path);if(!out)throw std::runtime_error("CSV open failed");
    out<<"task,proposal_index,uniform,pre_hash,ratio_real,ratio_imag,ratio_z2,decision,expected_decision,precision_digits,converged,pre_endpoint_rcond,post_endpoint_rcond,reference_seconds,continued_proposals,final_hash,mp_fallback_count,mp_fallback_failures,status\n";
    for(const Task&t:tasks){Replay r=replay(t);const auto&mp=r.frozen.ratio.mp_reference;
        out<<t.id<<','<<t.index<<','<<std::setprecision(17)<<t.uniform<<','<<t.pre_hash<<','
           <<r.frozen.ratio.ratio.real()<<','<<r.frozen.ratio.ratio.imag()<<','<<mp.ratio_z2<<','
           <<(r.frozen.accepted?"accept":"reject")<<','<<(t.accepted?"accept":"reject")<<','
           <<mp.precision_digits<<','<<mp.converged<<','<<mp.pre_endpoint_rcond<<','
           <<mp.post_endpoint_rcond<<','<<r.frozen.slow_seconds<<','<<r.continued<<','<<r.post_hash<<','<<r.mp_count<<','
           <<r.mp_failures<<",PASS\n";}
    out.flush();if(!out)throw std::runtime_error("CSV flush failed");out.close();if(out.fail())
        throw std::runtime_error("CSV close failed");std::cout<<"PASS frozen_task_4\nPASS frozen_task_8\n";return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
