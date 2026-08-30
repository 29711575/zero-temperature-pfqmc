#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pure_projector_fast.h"

namespace {

constexpr double kPi=3.141592653589793238462643383279502884;

void require(bool value,const std::string &message){if(!value)throw std::runtime_error(message);}
double matrixError(const MatType&a,const MatType&b){return(a-b).norm()/std::max(1.0,b.norm());}

struct CheckedCsv{
    std::string path;std::ofstream out;bool finished=false;
    explicit CheckedCsv(const std::string&p):path(p),out(p){if(!out)throw std::runtime_error("open failed: "+p);}
    void finish(){out.flush();if(!out)throw std::runtime_error("flush failed: "+path);out.close();if(out.fail())throw std::runtime_error("close failed: "+path);finished=true;}
    ~CheckedCsv(){if(!finished&&out.is_open())out.close();}
};

MatType trialHamiltonian(const SpinlessTvChainUtils&config,double split){
    MatType h=MatType::Zero(config.nDim,config.nDim);config.KineticGenerator(h);
    if(split!=0){int a=config.majoranaCoord2Idx(0,1),b=config.majoranaCoord2Idx(config.Lx-1,1);
        h(a,b)+=DataType(0,split);h(b,a)-=DataType(0,split);}return h;
}

MatType exponential(const MatType&generator){MatType copy=generator;return expm(copy,1.0);}

MatType localHsFactor(const SpinlessTvChainUtils&config,int bond,int aux,int sigma){
    MatType a=MatType::Zero(config.nDim,config.nDim);
    const double lambda=std::acosh(std::exp(0.5*config.V*config.dt));
    int i1,i2,i3,i4;config.aux2MajoranaIdx(aux,0,bond,i1,i2);config.aux2MajoranaIdx(aux,1,bond,i3,i4);
    const DataType z(0,lambda*sigma);
    if(config.hsScheme==0){a(i1,i2)=z;a(i2,i1)=-z;a(i3,i4)=z;a(i4,i3)=-z;}
    else{a(i1,i3)=z;a(i3,i1)=-z;a(i2,i4)=-z;a(i4,i2)=z;}
    return exponential(a);
}

std::pair<int,int> bondCounts(const SpinlessTvChainUtils&c){
    return pureProjectorCheckerboardBondCounts(c.Lx,c.boundaryType);
}

PureFastConfiguration makeConfiguration(const SpinlessTvChainUtils&config,unsigned seed){
    std::mt19937 random(seed);std::uniform_int_distribution<int>bit(0,1);
    MatType kinetic=MatType::Zero(config.nDim,config.nDim);config.KineticGenerator(kinetic);
    const MatType half=exponential(-0.5*config.dt*kinetic);
    const auto counts=bondCounts(config);
    auto protocol=[&](PureBranch branch,int protocolSlice){
        PureFastConfiguration part;
        auto push=[&](const MatType&m,int field,int factor,int bond,int aux,const std::string&label){
            part.slices.emplace_back(m,1.0,label);part.hs_fields.push_back(field);
            part.locations.push_back({branch,protocolSlice,factor,bond,aux});};
        int factor=0;push(half,0,factor++,-1,-1,"K/2");
        for(int bond=0;bond<2;++bond){int count=bond==0?counts.first:counts.second;
            for(int aux=0;aux<count;++aux){int sigma=bit(random)?1:-1;
                push(localHsFactor(config,bond,aux,sigma),sigma,factor++,bond,aux,
                     std::string("V")+std::to_string(bond)+":"+std::to_string(aux));}}
        push(half,0,factor++,-1,-1,"K/2");return part;};
    PureFastConfiguration ket=protocol(PureBranch::Ket,0),braProtocol=protocol(PureBranch::Bra,0),full=ket;
    for(int i=int(braProtocol.slices.size())-1;i>=0;--i){
        full.slices.push_back(braProtocol.slices[i]);full.hs_fields.push_back(braProtocol.hs_fields[i]);
        full.locations.push_back(braProtocol.locations[i]);}
    return full;
}

std::vector<int> proposalIndices(const PureFastConfiguration&c){std::vector<int>x;
    for(int i=0;i<int(c.locations.size());++i)if(c.locations[i].aux>=0)x.push_back(i);return x;}

PureFastProposal flipProposal(const SpinlessTvChainUtils&model,const PureFastConfiguration&c,
                              int index,double uniform){
    const auto &loc=c.locations[index];const int next=-c.hs_fields[index];
    PureFastProposal p;p.index=index;p.new_hs=next;p.uniform=uniform;p.new_eta=1.0;
    p.new_factor=localHsFactor(model,loc.bond,loc.aux,next);return p;
}

struct Observables{double spi=0,sdq=0,rcdw=0;};
DataType structure(const MatType&g,int L,double q){
    return pureProjectorStructureFactor(g,L,q);}
Observables observe(const MatType&g,int L){DataType p=structure(g,L,kPi),d=structure(g,L,kPi-2*kPi/L);
    return{p.real(),d.real(),(1.0-d/p).real()};}

struct Metrics{
    int passed=3;double ratioMax=0,greenMax=0;long long z2Mismatch=0,trajectoryMismatch=0;
    bool obcAlwaysPlus=true;double fastSeconds=0,slowSeconds=0;double stressMax=0;int firstAlarm=-1;
};

void singleProposalChecks(CheckedCsv&csv,Metrics&m){
    csv.out<<"L,boundary,hs_scheme,proposal,index,branch,bond,aux,low_rank,fast_real,fast_imag,slow_real,slow_imag,ratio_relative_error,green_relative_error,z2_mismatch,accepted,uniform,trust_alarm,used_reference,fast_seconds,slow_seconds\n";
    std::mt19937_64 random(31001);std::uniform_real_distribution<double>uniform(0,1);
    for(int L:{2,4}){
        SpinlessTvChainUtils model(L,.08,1.1,2,L==2?1:0,.7,.31,L==2?1:0);
        MatType h=trialHamiltonian(model,0);GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(h);
        PureFastConfiguration config=makeConfiguration(model,32000+L);
        PureProjectorFastWalker walker(trial,config,2,PureFastRunMode::AuditLockstep);
        for(int step=0;step<500;++step){auto indices=proposalIndices(walker.configuration());int index=indices[random()%indices.size()];
            PureFastProposal p=flipProposal(model,walker.configuration(),index,uniform(random));
            PureFastProposalResult r=walker.propose(p);require(!r.terminated&&r.ratio.ok(),"single proposal failed");
            require(r.ratio.relative_reference_error<1e-10,"single proposal ratio mismatch");
            require(r.ratio.green_reference_error<1e-10,"single proposal Green mismatch");
            m.ratioMax=std::max(m.ratioMax,r.ratio.relative_reference_error);
            m.greenMax=std::max(m.greenMax,r.ratio.green_reference_error);m.z2Mismatch+=r.z2_reference_mismatch;
            m.fastSeconds+=r.fast_seconds;m.slowSeconds+=r.slow_seconds;
            const auto &loc=walker.configuration().locations[index];
            csv.out<<L<<','<<model.boundaryType<<','<<model.hsScheme<<','<<step<<','<<index<<','
                   <<(loc.branch==PureBranch::Ket?"ket":"bra")<<','<<loc.bond<<','<<loc.aux<<','
                   <<r.ratio.low_rank<<','<<std::setprecision(17)<<r.ratio.ratio.real()<<','<<r.ratio.ratio.imag()<<','
                   <<r.ratio.slow_ratio.real()<<','<<r.ratio.slow_ratio.imag()<<','
                   <<r.ratio.relative_reference_error<<','<<r.ratio.green_reference_error<<','
                   <<r.z2_reference_mismatch<<','<<r.accepted<<','<<p.uniform<<','<<r.ratio.trust_alarm<<','
                   <<r.ratio.used_reference<<','<<r.fast_seconds<<','<<r.slow_seconds<<'\n';
        }
    }
    require(m.z2Mismatch==0,"single proposal Z2 mismatch");m.passed++;
}

void decisionCases(Metrics&m){
    SpinlessTvChainUtils model(2,.08,.7,2,1,.6,.2,0);MatType h=trialHamiltonian(model,0);
    GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(h);PureFastConfiguration base=makeConfiguration(model,33001);
    auto idx=proposalIndices(base);require(!idx.empty(),"no decision proposal");
    {PureProjectorFastWalker w(trial,base,1,PureFastRunMode::AuditLockstep);auto p=flipProposal(model,w.configuration(),idx[0],0.0);auto r=w.propose(p);require(r.accepted,"positive case not accepted");}
    {PureProjectorFastWalker w(trial,base,1,PureFastRunMode::AuditLockstep);PureFastProposal p;
      p.index=idx[0];p.new_hs=base.hs_fields[idx[0]];p.new_factor=base.slices[idx[0]].matrix;p.new_eta=-1;p.uniform=.2;
      auto r=w.propose(p);require(r.accepted&&w.z2Sign()==-1,"negative ratio Z2 transport failed");}
    {PureProjectorFastWalker w(trial,base,1,PureFastRunMode::AuditLockstep);auto p=flipProposal(model,w.configuration(),idx[0],.999999999);auto r=w.propose(p);require(!r.accepted,"rejection case accepted");}
    {PureFastOptions o;o.decision_margin_tolerance=1e-8;PureProjectorFastWalker w(trial,base,1,PureFastRunMode::FastStrict,o);
      auto p=flipProposal(model,w.configuration(),idx[0],0);PureFastRatioResult q=pureFastLocalRatio(trial,base.slices,p.index,p.new_factor,1.0,o);
      p.uniform=std::min(1.0,std::abs(q.ratio));if(p.uniform>=1)p.uniform=std::nextafter(1.0,0.0);auto r=w.propose(p);
      require(r.ratio.used_reference&&r.ratio.trust_alarm,"decision-margin fallback missing");}
    {PureProjectorFastWalker w(trial,base,1,PureFastRunMode::FastStrict);PureFastProposal p;
      p.index=idx[0];p.new_hs=-base.hs_fields[idx[0]];p.uniform=.2;p.new_factor=MatType::Zero(model.nDim,model.nDim);
      auto r=w.propose(p);require(r.terminated&&r.ratio.status==PureFastRatioStatus::reference_failure,"zero proposal did not fail closed");}
    {PureProjectorFastWalker w(trial,base,1,PureFastRunMode::FastStrict);PureFastProposal p;
      p.index=idx[0];p.new_hs=base.hs_fields[idx[0]];p.uniform=.2;p.new_factor=base.slices[idx[0]].matrix;p.new_eta=DataType(0,1);
      auto r=w.propose(p);require(r.ratio.used_reference,"complex ratio did not invoke reference");}
    m.passed++;
}

void trajectoryChecks(CheckedCsv&csv,Metrics&m){
    csv.out<<"L,step,index,uniform,audit_accept,fast_accept,audit_hash,fast_hash,audit_z2,fast_z2,green_error,S_pi_error,S_pi_dq_error,R_CDW_error,trust_alarm,trajectory_mismatch\n";
    for(int L:{2,4}){SpinlessTvChainUtils model(L,.07,1.0,2,L==2?1:0,.65,.27,L==2?1:0);
        MatType h=trialHamiltonian(model,0);GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(h);
        PureFastConfiguration base=makeConfiguration(model,34000+L);
        PureProjectorFastWalker audit(trial,base,2,PureFastRunMode::AuditLockstep),fast(trial,base,2,PureFastRunMode::FastStrict);
        std::mt19937_64 rng(35000+L);std::uniform_real_distribution<double>u(0,1);
        for(int step=0;step<1000;++step){auto ids=proposalIndices(audit.configuration());int index=ids[rng()%ids.size()];double uniform=u(rng);
            auto pa=flipProposal(model,audit.configuration(),index,uniform),pf=flipProposal(model,fast.configuration(),index,uniform);
            auto ra=audit.propose(pa),rf=fast.propose(pf);require(!ra.terminated&&!rf.terminated,"trajectory terminated");
            auto ga=audit.measurementGreen(),gf=fast.measurementGreen();require(ga.ok()&&gf.ok(),"trajectory Green failed");
            double ge=matrixError(gf.green,ga.green);Observables oa=observe(ga.green,L),of=observe(gf.green,L);
            bool mismatch=ra.accepted!=rf.accepted||audit.configurationHash()!=fast.configurationHash()||
                audit.z2Sign()!=fast.z2Sign()||ge>1e-10;m.trajectoryMismatch+=mismatch;m.greenMax=std::max(m.greenMax,ge);
            csv.out<<L<<','<<step<<','<<index<<','<<std::setprecision(17)<<uniform<<','<<ra.accepted<<','<<rf.accepted<<','
                   <<audit.configurationHash()<<','<<fast.configurationHash()<<','<<audit.z2Sign()<<','<<fast.z2Sign()<<','<<ge<<','
                   <<std::abs(of.spi-oa.spi)<<','<<std::abs(of.sdq-oa.sdq)<<','<<std::abs(of.rcdw-oa.rcdw)<<','
                   <<rf.ratio.trust_alarm<<','<<mismatch<<'\n';
        }
    }require(m.trajectoryMismatch==0,"full trajectory mismatch");m.passed++;
}

void l6Regression(CheckedCsv&csv,Metrics&m){
    csv.out<<"case,seed,proposals,audit_average_sign,fast_average_sign,audit_S_pi,fast_S_pi,audit_S_pi_dq,fast_S_pi_dq,audit_R_CDW,fast_R_CDW,audit_acceptance,fast_acceptance,z2_always_plus,trajectory_mismatch,fast_seconds,slow_seconds\n";
    for(int kind=0;kind<2;++kind)for(int seed=0;seed<4;++seed){int boundary=kind==0?0:1,hs=kind==0?0:1;
        SpinlessTvChainUtils model(6,.06,.9,2,boundary,kind==0?.7:1.0,kind==0?.25:0.0,hs);
        double split=kind==0?0:1e-8;MatType h=trialHamiltonian(model,split);
        GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(h);PureFastConfiguration base=makeConfiguration(model,36000+100*kind+seed);
        PureProjectorFastWalker audit(trial,base,2,PureFastRunMode::AuditLockstep),fast(trial,base,2,PureFastRunMode::FastStrict);
        std::mt19937_64 rng(37000+100*kind+seed);std::uniform_real_distribution<double>u(0,1);
        double za=0,zf=0,spa=0,spf=0,sda=0,sdf=0,ra=0,rf=0;int aa=0,af=0,n=200;bool plus=true;long long mismatch=0;double ft=0,st=0;
        for(int step=0;step<n;++step){auto ids=proposalIndices(audit.configuration());int index=ids[rng()%ids.size()];double uniform=u(rng);
            auto a=audit.propose(flipProposal(model,audit.configuration(),index,uniform));
            auto f=fast.propose(flipProposal(model,fast.configuration(),index,uniform));
            require(!a.terminated&&!f.terminated,"L6 trajectory terminated");aa+=a.accepted;af+=f.accepted;ft+=f.fast_seconds;st+=a.slow_seconds;
            mismatch+=a.accepted!=f.accepted||audit.configurationHash()!=fast.configurationHash()||audit.z2Sign()!=fast.z2Sign();
            auto ga=audit.measurementGreen(),gf=fast.measurementGreen();Observables oa=observe(ga.green,6),of=observe(gf.green,6);
            za+=audit.z2Sign();zf+=fast.z2Sign();spa+=oa.spi;spf+=of.spi;sda+=oa.sdq;sdf+=of.sdq;ra+=oa.rcdw;rf+=of.rcdw;
            if(kind==1&&(audit.z2Sign()!=1||fast.z2Sign()!=1))plus=false;}
        if(kind==1)m.obcAlwaysPlus=m.obcAlwaysPlus&&plus;m.trajectoryMismatch+=mismatch;m.fastSeconds+=ft;m.slowSeconds+=st;
        csv.out<<(kind==0?"PBC_hs0":"OBC_hs1")<<','<<seed<<','<<n<<','<<za/n<<','<<zf/n<<','<<spa/n<<','<<spf/n<<','
               <<sda/n<<','<<sdf/n<<','<<ra/n<<','<<rf/n<<','<<double(aa)/n<<','<<double(af)/n<<','<<plus<<','<<mismatch<<','<<ft<<','<<st<<'\n';
    }require(m.obcAlwaysPlus,"OBC hs1 Z2 was not always +1");require(m.trajectoryMismatch==0,"L6 mismatch");m.passed++;
}

MatType rotation(int n,int a,int b,double x){MatType m=MatType::Identity(n,n);m(a,a)=std::cos(x);m(b,b)=std::cos(x);m(a,b)=std::sin(x);m(b,a)=-std::sin(x);return m;}
MatType canonicalPhi(int modes){MatType p=MatType::Zero(2*modes,modes);for(int i=0;i<modes;++i){p(2*i,i)=1/std::sqrt(2.);p(2*i+1,i)=DataType(0,1/std::sqrt(2.));}return p;}

void stabilization(CheckedCsv&csv,Metrics&m){csv.out<<"block_size,cut,green_relative_error,overlap_rcond,solve_residual,alarm\n";
    GaussianTrialState trial=GaussianTrialState::fromPhi(canonicalPhi(6));std::vector<PureProjectorSlice>s;
    for(int i=0;i<48;++i)s.emplace_back(rotation(12,i%12,(i*5+3)%12,.012+.001*(i%5)),1.0,"stress");
    for(int block:{1,2,4,8}){PureProjectorStackManager stack(trial,s,block);for(int cut=0;cut<=48;++cut){require(stack.moveToCut(cut),"stress cut failed");
        MatType r=trial.Phi,l=trial.Phi;for(int i=0;i<cut;++i)r=s[i].matrix*r;for(int i=47;i>=cut;--i)l=s[i].matrix.adjoint()*l;
        auto direct=pureProjectorGreenThinQr(r,l),got=stack.green();double e=matrixError(got.green,direct.green);bool alarm=e>1e-10||!got.ok();
        m.stressMax=std::max(m.stressMax,e);if(alarm&&m.firstAlarm<0)m.firstAlarm=cut;
        csv.out<<block<<','<<cut<<','<<std::setprecision(17)<<e<<','<<got.overlap_rcond<<','<<got.solve_residual<<','<<alarm<<'\n';}}
    require(m.firstAlarm<0,"normal stabilization stress alarmed");m.passed++;
}

} // namespace

int main(int argc,char**argv){try{if(argc!=5)throw std::invalid_argument("usage: phase3a_validation OUTPUT_DIR SOURCE_COMMIT EXE_SHA THREADS");
    std::string out=argv[1];CheckedCsv ratio(out+"/proposal_ratio_checks.csv"),traj(out+"/trajectory_lockstep.csv"),block(out+"/block_stabilization.csv"),mc(out+"/fast_vs_slow_mc.csv");Metrics metrics;
    singleProposalChecks(ratio,metrics);decisionCases(metrics);trajectoryChecks(traj,metrics);l6Regression(mc,metrics);stabilization(block,metrics);
    require(metrics.passed==8,"not all validation groups ran");ratio.finish();traj.finish();block.finish();mc.finish();
    std::cout<<std::setprecision(17)<<"{\"status\":\"complete\",\"tests_passed\":"<<metrics.passed<<",\"tests_total\":8,\"source_commit\":\""<<argv[2]
             <<"\",\"executable_sha256\":\""<<argv[3]<<"\",\"ratio_error_max\":"<<metrics.ratioMax<<",\"green_error_max\":"<<metrics.greenMax
             <<",\"z2_mismatch_count\":"<<metrics.z2Mismatch<<",\"trajectory_mismatch_count\":"<<metrics.trajectoryMismatch
             <<",\"obc_hs1_always_plus\":"<<(metrics.obcAlwaysPlus?"true":"false")<<",\"fast_seconds\":"<<metrics.fastSeconds
             <<",\"slow_seconds\":"<<metrics.slowSeconds<<",\"fast_over_slow\":"<<metrics.fastSeconds/metrics.slowSeconds
             <<",\"stabilization_error_max\":"<<metrics.stressMax<<",\"first_alarm\":"<<metrics.firstAlarm<<"}\n";return 0;
}catch(const std::exception&e){std::cerr<<"phase3a_validation: "<<e.what()<<'\n';return 1;}}
