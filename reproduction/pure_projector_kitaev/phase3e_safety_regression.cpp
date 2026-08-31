#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pure_projector_fast.h"

namespace {

void require(bool ok,const std::string&message){if(!ok)throw std::runtime_error(message);}

MatType kineticGenerator(int L,double delta,double mu){
    auto make=[&](double d,double m){SpinlessTvChainUtils c(L,1,0,2,0,d,m,0);
        MatType h=MatType::Zero(2*L,2*L);c.KineticGenerator(h);return h;};
    MatType base=make(0,0);return base+delta*(make(1,0)-base)+mu*(make(0,1)-base);
}
MatType exponential(MatType generator,double scale=1){return expm(generator,scale);}
MatType localHsFactor(const SpinlessTvChainUtils&c,int bond,int aux,int sigma){
    MatType generator=MatType::Zero(c.nDim,c.nDim);
    const double lambda=std::acosh(std::exp(.5*c.V*c.dt));int a,b,d,e;
    c.aux2MajoranaIdx(aux,0,bond,a,b);c.aux2MajoranaIdx(aux,1,bond,d,e);
    const DataType z(0,lambda*sigma);
    if(c.hsScheme==0){generator(a,b)=z;generator(b,a)=-z;generator(d,e)=z;generator(e,d)=-z;}
    else{generator(a,d)=z;generator(d,a)=-z;generator(b,e)=-z;generator(e,b)=z;}
    return exponential(generator);
}

struct FrozenTask {int id,L;double V;std::uint64_t seed;int attempt,index;
    double uniform,ratio;std::uint64_t pre_hash;bool accepted;};
const FrozenTask frozenTasks[]={
    {4,6,4,706041,110,146,.52928038756354223,.707419140654256205736,7807527216050905631ULL,true},
    {8,6,6,706061,415,553,.56577673777691728,.002657388953068520511,5178673164861548103ULL,false}};

GaussianTrialState frozenTrial(const FrozenTask&t){auto trial=GaussianTrialState::fromMajoranaHamiltonian(
    kineticGenerator(t.L,1,.3));require(trial.fermionParity()==-1,
    "frozen physical trial parity");return trial;}

PureFastConfiguration independentContour(int L,double V,std::mt19937_64&rng,int timeSlices=-1){
    SpinlessTvChainUtils model(L,.1,V,2,0,1,0,0);
    if(timeSlices<0)timeSlices=10*L;const MatType half=exponential(kineticGenerator(L,1,0),-.05);
    const int perLayer=L/2;std::uniform_int_distribution<int>bit(0,1);
    auto branch=[&](PureBranch side){PureFastConfiguration out;int ordinal=0;
        auto push=[&](const MatType&m,int field,int slice,int bond,int aux,const std::string&label){
            out.slices.emplace_back(m,1.0,label);out.hs_fields.push_back(field);
            out.locations.push_back({side,slice,ordinal++,bond,aux});};
        for(int slice=0;slice<timeSlices;++slice){push(half,0,slice,-1,-1,"K/2");
            for(int layer=0;layer<2;++layer)for(int aux=0;aux<perLayer;++aux){int s=bit(rng)?1:-1;
                push(localHsFactor(model,layer,aux,s),s,slice,layer,aux,"V");}
            push(half,0,slice,-1,-1,"K/2");}return out;};
    auto ket=branch(PureBranch::Ket),bra=branch(PureBranch::Bra),full=ket;
    for(int i=int(bra.slices.size())-1;i>=0;--i){full.slices.push_back(bra.slices[i]);
        full.hs_fields.push_back(bra.hs_fields[i]);full.locations.push_back(bra.locations[i]);}
    return full;
}

std::vector<int> proposalIndices(const PureFastConfiguration&c){std::vector<int> out;
    for(int i=0;i<int(c.locations.size());++i)if(c.locations[i].aux>=0)out.push_back(i);return out;}
PureFastProposal flip(const SpinlessTvChainUtils&model,const PureFastConfiguration&c,int index,double u){
    const auto&loc=c.locations[index];PureFastProposal p;p.index=index;p.new_hs=-c.hs_fields[index];
    p.new_factor=localHsFactor(model,loc.bond,loc.aux,p.new_hs);p.new_eta=1;p.uniform=u;return p;}

struct Checkpoint {PureFastConfiguration configuration;std::mt19937_64 rng;std::size_t cursor=0;
    int direction=1,z2=1,attempted=0;double log_abs_weight=0;DataType phase=DataType(1,0);};
void writeCheckpoint(const std::string&path,const Checkpoint&c){std::ofstream o(path,std::ios::binary);
    require(bool(o),"checkpoint open");o<<std::setprecision(17)<<c.cursor<<' '<<c.direction<<' '<<c.z2<<' '<<c.attempted<<' '
        <<c.log_abs_weight<<' '<<c.phase.real()<<' '<<c.phase.imag()<<'\n'<<c.rng<<'\n';
    o<<c.configuration.slices.size()<<'\n'<<std::setprecision(17);
    for(std::size_t k=0;k<c.configuration.slices.size();++k){const auto&s=c.configuration.slices[k];const auto&l=c.configuration.locations[k];
        o<<c.configuration.hs_fields[k]<<' '<<int(l.branch)<<' '<<l.slice<<' '<<l.factor<<' '<<l.bond<<' '<<l.aux<<' '
         <<s.eta.real()<<' '<<s.eta.imag()<<' '<<s.matrix.rows()<<' '<<s.matrix.cols()<<' '<<s.label.size()<<'\n'<<s.label<<'\n';
        for(int i=0;i<s.matrix.rows();++i)for(int j=0;j<s.matrix.cols();++j)o<<s.matrix(i,j).real()<<' '<<s.matrix(i,j).imag()<<'\n';}
    o.flush();require(bool(o),"checkpoint write");o.close();require(!o.fail(),"checkpoint close");}
Checkpoint readCheckpoint(const std::string&path){std::ifstream in(path,std::ios::binary);require(bool(in),"checkpoint read open");
    Checkpoint c;double pr,pi;in>>c.cursor>>c.direction>>c.z2>>c.attempted>>c.log_abs_weight>>pr>>pi;
    c.phase=DataType(pr,pi);in>>c.rng;std::size_t n;in>>n;
    for(std::size_t k=0;k<n;++k){int hs,branch,slice,factor,bond,aux,rows,cols;double er,ei;std::size_t labelSize;
        in>>hs>>branch>>slice>>factor>>bond>>aux>>er>>ei>>rows>>cols>>labelSize;in.get();std::string label;std::getline(in,label);
        require(label.size()==labelSize,"checkpoint label");MatType m(rows,cols);for(int i=0;i<rows;++i)for(int j=0;j<cols;++j){double r,x;in>>r>>x;m(i,j)=DataType(r,x);}
        c.configuration.slices.emplace_back(m,DataType(er,ei),label);c.configuration.hs_fields.push_back(hs);
        c.configuration.locations.push_back({branch?PureBranch::Bra:PureBranch::Ket,slice,factor,bond,aux});}
    require(bool(in),"checkpoint parse");return c;}

std::string greenLine(int step,const PureFastProposalResult&r,const PureProjectorFastWalker&w){
    auto g=const_cast<PureProjectorFastWalker&>(w).measurementGreen();require(g.ok(),"measurement Green");
    std::ostringstream o;o<<step<<','<<r.snapshot.slice<<','<<r.snapshot.factor<<','<<std::setprecision(17)
      <<r.snapshot.uniform<<','<<r.accepted<<','<<w.configurationHash()<<','<<w.z2Sign();
    for(int i=0;i<g.green.rows();++i)for(int j=0;j<g.green.cols();++j)o<<','<<g.green(i,j).real()<<','<<g.green(i,j).imag();
    o<<'\n';return o.str();}

struct TailResult {std::string csv;PureFastProposalResult first;std::uint64_t final_hash=0,rng_hash=0;int z2=0;};
std::uint64_t hashRng(const std::mt19937_64&rng){std::ostringstream s;s<<rng;std::uint64_t h=1469598103934665603ULL;
    for(unsigned char c:s.str()){h^=c;h*=1099511628211ULL;}return h;}
TailResult runTail(const FrozenTask&t,const GaussianTrialState&trial,const SpinlessTvChainUtils&model,
                   Checkpoint checkpoint,int count){
    PureFastOptions options;options.read_only_audit_interval=0;PureProjectorFastWalker w(trial,checkpoint.configuration,8,
        PureFastRunMode::FastStrict,options,PureFastInitializationPolicy::MirroredTheoremZ2Plus);
    w.restoreTrustedCheckpointForTest(checkpoint.z2,checkpoint.log_abs_weight,checkpoint.phase);
    require(w.z2Sign()==checkpoint.z2,"restart Z2 mismatch");std::uniform_real_distribution<double>u(0,1);TailResult out;
    for(int k=0;k<count;++k){auto ids=proposalIndices(w.configuration());int index=ids[checkpoint.cursor];
        if(ids.size()>1){if(checkpoint.direction>0&&checkpoint.cursor+1==ids.size())checkpoint.direction=-1;
            else if(checkpoint.direction<0&&checkpoint.cursor==0)checkpoint.direction=1;
            checkpoint.cursor=std::size_t(int(checkpoint.cursor)+checkpoint.direction);}double uniform=u(checkpoint.rng);
        auto r=w.propose(flip(model,w.configuration(),index,uniform));require(!r.terminated&&r.ratio.ok(),"restart tail proposal");
        if(k==0)out.first=r;out.csv+=greenLine(k,r,w);}
    out.final_hash=w.configurationHash();out.rng_hash=hashRng(checkpoint.rng);out.z2=w.z2Sign();return out;}

void restartTest(const std::string&outdir){std::ofstream summary(outdir+"/restart_determinism.csv");
    summary<<"task,checkpoint_attempt,proposal_index,uniform,ratio,decision,proposal_identity_match,uniform_match,decision_match,hs_hash_match,rng_hash_match,z2_match,green_csv_match,measurement_csv_byte_identical,status\n";
    for(const auto&t:frozenTasks){auto trial=frozenTrial(t);SpinlessTvChainUtils model(t.L,.1,t.V,2,0,1,0,0);std::mt19937_64 rng(t.seed);
        std::unique_ptr<PureProjectorFastWalker>w;PureFastConfiguration initial;
        while(!w){initial=independentContour(t.L,t.V,rng);try{w.reset(new PureProjectorFastWalker(trial,initial,8,PureFastRunMode::FastStrict));}catch(const std::runtime_error&){}}
        std::uniform_real_distribution<double>u(0,1);std::size_t cursor=0;int direction=1;
        for(int attempted=1;attempted<t.attempt;++attempted){auto ids=proposalIndices(w->configuration());int index=ids[cursor];
            if(ids.size()>1){if(direction>0&&cursor+1==ids.size())direction=-1;else if(direction<0&&cursor==0)direction=1;cursor=std::size_t(int(cursor)+direction);}
            auto r=w->propose(flip(model,w->configuration(),index,u(rng)));require(!r.terminated&&r.ratio.ok(),"checkpoint prefix");}
        Checkpoint cp{w->configuration(),rng,cursor,direction,w->z2Sign(),t.attempt-1,
            w->currentWeight().log_abs_weight,w->currentWeight().complex_phase};
        std::string path=outdir+"/restart_task"+std::to_string(t.id)+".checkpoint";
        writeCheckpoint(path,cp);Checkpoint disk=readCheckpoint(path);TailResult continuous=runTail(t,trial,model,cp,65),restart=runTail(t,trial,model,disk,65);
        require(continuous.first.ratio.mp_reference.ok()&&restart.first.ratio.mp_reference.ok(),"frozen first MP missing");
        bool identity=continuous.first.snapshot.configuration_hash==restart.first.snapshot.configuration_hash&&
            continuous.first.snapshot.slice==restart.first.snapshot.slice&&continuous.first.snapshot.factor==restart.first.snapshot.factor;
        bool uniformMatch=continuous.first.snapshot.uniform==restart.first.snapshot.uniform;
        bool decision=continuous.first.accepted==restart.first.accepted;bool hs=continuous.final_hash==restart.final_hash;
        bool rngMatch=continuous.rng_hash==restart.rng_hash,z2=continuous.z2==restart.z2,csv=continuous.csv==restart.csv;
        std::ofstream a(outdir+"/restart_task"+std::to_string(t.id)+"_continuous_measurements.csv");a<<continuous.csv;a.close();
        std::ofstream b(outdir+"/restart_task"+std::to_string(t.id)+"_restarted_measurements.csv");b<<restart.csv;b.close();
        require(identity&&uniformMatch&&decision&&hs&&rngMatch&&z2&&csv,"restart determinism");
        summary<<t.id<<','<<t.attempt-1<<','<<t.index<<','<<std::setprecision(17)<<t.uniform<<','<<continuous.first.ratio.ratio.real()<<','
          <<continuous.first.accepted<<','<<identity<<','<<uniformMatch<<','<<decision<<','<<hs<<','<<rngMatch<<','<<z2<<','<<csv<<','<<csv<<",PASS\n";}
    summary.flush();require(bool(summary),"restart summary write");}

void ratioTest(const std::string&outdir){std::ofstream out(outdir+"/ratio_reciprocity.csv");
    out<<"L,proposal,mode,index,forward_real,forward_imag,reverse_real,reverse_imag,product_real,product_imag,reciprocity_error,z2_product,detailed_balance_error,decision_failure,reference_trusted,status\n";
    std::mt19937_64 rng(930501);std::uniform_real_distribution<double>u(0,1);PureFastOptions options;
    for(int L:{4,6}){auto trial=GaussianTrialState::fromMajoranaHamiltonian(kineticGenerator(L,1,.3));auto config=independentContour(L,1.2,rng,2);
        SpinlessTvChainUtils model(L,.1,1.2,2,0,1,0,0);auto ids=proposalIndices(config);
        for(int p=0;p<500;++p){int index=ids[rng()%ids.size()];auto proposal=flip(model,config,index,u(rng));PureFastConfiguration candidate=config;
            candidate.slices[index].matrix=proposal.new_factor;candidate.slices[index].eta=proposal.new_eta;candidate.hs_fields[index]=proposal.new_hs;
            PureFastRatioResult ff=pureFastLocalRatio(trial,config.slices,index,proposal.new_factor,proposal.new_eta,options);
            PureFastRatioResult fr=pureFastLocalRatio(trial,candidate.slices,index,config.slices[index].matrix,config.slices[index].eta,options);
            require(ff.ok()&&fr.ok(),"fast reciprocity untrusted");DataType f=ff.ratio,r=fr.ratio;std::string mode="fast";bool trusted=true;
            auto oldW=pureProjectorStableReferenceWeight(trial,config.slices,options),newW=pureProjectorStableReferenceWeight(trial,candidate.slices,options);
            require(oldW.ok()&&newW.ok(),"slow reciprocity untrusted");DataType sf=std::exp(newW.log_abs_weight-oldW.log_abs_weight)*newW.complex_phase/oldW.complex_phase;
            DataType sr=std::exp(oldW.log_abs_weight-newW.log_abs_weight)*oldW.complex_phase/newW.complex_phase;
            if(p==0){auto mf=pureProjectorMpSameProposal(trial,config.slices,index,proposal.new_factor,proposal.new_eta);
                auto mr=pureProjectorMpSameProposal(trial,candidate.slices,index,config.slices[index].matrix,config.slices[index].eta);
                require(mf.ok()&&mr.ok(),"forced MP reciprocity untrusted");f=mf.ratio;r=mr.ratio;mode="forced_mp";}
            else if(p%2){f=sf;r=sr;mode="slow_reference";}
            DataType product=f*r;double error=std::abs(product-DataType(1,0));int zprod=(f.real()>=0?1:-1)*(r.real()>=0?1:-1);
            double pf=std::min(1.0,std::abs(f)),pr=std::min(1.0,std::abs(r));double db=std::abs(pf/pr-std::abs(f));
            double uniform=u(rng);bool af=uniform<pf,ar=uniform<pr;bool decisionFailure=(af!=(uniform<pf))||(ar!=(uniform<pr));
            require(error<2e-9&&zprod==1&&db<2e-9&&!decisionFailure,"reciprocity/detailed balance");
            out<<L<<','<<p<<','<<mode<<','<<index<<','<<std::setprecision(17)<<f.real()<<','<<f.imag()<<','<<r.real()<<','<<r.imag()<<','
               <<product.real()<<','<<product.imag()<<','<<error<<','<<zprod<<','<<db<<','<<decisionFailure<<','<<trusted<<",PASS\n";
            if(p%7==0)config=std::move(candidate);}
        const int badIndex=ids.front();MatType zero=MatType::Zero(2*L,2*L);
        auto badMp=pureProjectorMpSameProposal(trial,config.slices,badIndex,zero,DataType(1,0));
        PureFastOptions strict;strict.decision_margin_tolerance=1.0;
        PureProjectorFastWalker failClosed(trial,config,2,PureFastRunMode::FastStrict,strict,
            PureFastInitializationPolicy::SequentialAudit);
        PureFastProposal bad;bad.index=badIndex;bad.new_hs=-config.hs_fields[badIndex];
        bad.new_factor=zero;bad.new_eta=1;bad.uniform=.5;auto badResult=failClosed.propose(bad);
        require(!badMp.ok()&&badResult.terminated&&
            badResult.ratio.status==PureFastRatioStatus::reference_failure,
            "untrusted MP reference did not fail closed");
        out<<L<<",500,untrusted_reference,"<<badIndex
           <<",0,0,0,0,0,0,0,0,0,0,0,PASS\n";
    }out.flush();require(bool(out),"ratio summary write");}

} // namespace

int main(int argc,char**argv){try{require(argc==3,"usage: phase3c_local_extra MODE OUTDIR");std::string mode=argv[1],out=argv[2];
    if(mode=="restart")restartTest(out);else if(mode=="ratio")ratioTest(out);else throw std::runtime_error("unknown mode");
    std::cout<<"PASS "<<mode<<'\n';return 0;}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
