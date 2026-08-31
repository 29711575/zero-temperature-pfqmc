#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pure_projector_fast.h"

namespace {

constexpr double kPi=3.141592653589793238462643383279502884;

struct Parameters {
    int L=0,boundary=-1,hs=-1,trial_parity=0,burn=0,measurements=0,block=0;
    int measurement_stride=1;
    double V=0,t=1,delta=0,mu=0,theta=0,dt=0;
    double trial_t=1,trial_delta=0,trial_mu=0,edge_splitting=0;
    std::uint64_t seed=0,audit_interval=0;
    std::string retained="off",source_commit="unknown",executable_sha256="unknown";
    std::string walker_mode="fast-strict",initialization_policy="mirrored-theorem";
    std::string run_units="sweeps";
};

std::string jsonEscape(const std::string&s){std::ostringstream o;for(unsigned char c:s){
    if(c=='"'||c=='\\')o<<'\\'<<c;else if(c=='\n')o<<"\\n";else if(c<32)o<<"?";else o<<c;}return o.str();}

Parameters parse(int argc,char**argv){
    std::map<std::string,std::string> values;
    for(int i=1;i<argc;++i){std::string key=argv[i];if(key.rfind("--",0)!=0||i+1>=argc)
        throw std::invalid_argument("arguments must be --name value pairs");values[key.substr(2)]=argv[++i];}
    auto get=[&](const char*k)->std::string{auto it=values.find(k);if(it==values.end())
        throw std::invalid_argument(std::string("missing --")+k);return it->second;};
    // Reject odd periodic rings before requiring the remaining configuration;
    // this is also a cheap fail-fast contract probe.
    if(values.count("L")&&values.count("boundary")&&std::stoi(values["L"])%2&&
       (values["boundary"]=="pbc"||values["boundary"]=="0"))
        throw std::invalid_argument("odd-L PBC is unsupported");
    Parameters p;p.L=std::stoi(get("L"));p.V=std::stod(get("V"));p.t=std::stod(get("t"));
    p.delta=std::stod(get("delta"));p.mu=std::stod(get("mu"));p.theta=std::stod(get("theta"));
    p.dt=std::stod(get("dt"));std::string boundary=get("boundary");
    p.boundary=boundary=="pbc"||boundary=="0"?0:boundary=="obc"||boundary=="1"?1:-1;
    std::string hs=get("hs-scheme");p.hs=hs=="hs0"||hs=="0"?0:hs=="hs1"||hs=="1"?1:-1;
    p.trial_t=std::stod(get("trial-t"));p.trial_delta=std::stod(get("trial-delta"));
    p.trial_mu=std::stod(get("trial-mu"));p.trial_parity=std::stoi(get("trial-parity"));
    p.edge_splitting=std::stod(get("edge-splitting"));p.burn=std::stoi(get("burn"));
    p.measurements=std::stoi(get("measurements"));p.seed=std::stoull(get("seed"));
    p.block=std::stoi(get("stabilization-block"));p.retained=get("retained");
    if(values.count("audit-interval"))p.audit_interval=std::stoull(values["audit-interval"]);
    if(values.count("walker-mode"))p.walker_mode=values["walker-mode"];
    if(values.count("initialization-policy"))
        p.initialization_policy=values["initialization-policy"];
    if(values.count("run-units"))p.run_units=values["run-units"];
    if(values.count("measurement-stride"))
        p.measurement_stride=std::stoi(values["measurement-stride"]);
    if(values.count("source-commit"))p.source_commit=values["source-commit"];
    if(values.count("executable-sha256"))p.executable_sha256=values["executable-sha256"];
    if(p.L<2||p.boundary<0||p.hs<0||p.dt<=0||p.theta<0||p.V<0||p.burn<0||
       p.measurements<=0||p.block<=0||(p.trial_parity!=1&&p.trial_parity!=-1)||
       (p.walker_mode!="fast-strict"&&p.walker_mode!="audit-lockstep")||
       (p.run_units!="sweeps"&&p.run_units!="proposals")||p.measurement_stride<=0||
       (p.initialization_policy!="mirrored-theorem"&&
        p.initialization_policy!="sequential-audit")||
       (p.initialization_policy=="sequential-audit"&&p.L>6))
        throw std::invalid_argument("invalid production parameter");
    double slices=p.theta/p.dt;if(std::abs(slices-std::round(slices))>1e-10*std::max(1.0,slices))
        throw std::invalid_argument("theta/dt must be integral");
    return p;
}

MatType kineticGenerator(int L,int boundary,double t,double delta,double mu){
    auto matrix=[&](double d,double m){SpinlessTvChainUtils c(L,1,0,2,boundary,d,m,0);
        MatType h=MatType::Zero(2*L,2*L);c.KineticGenerator(h);return h;};
    MatType base=matrix(0,0),pair=matrix(1,0)-base,chemical=matrix(0,1)-base;
    return t*base+delta*pair+mu*chemical;
}

MatType exponential(MatType generator,double scale){return expm(generator,scale);}

std::pair<int,int> bondCounts(const SpinlessTvChainUtils&c){
    return pureProjectorCheckerboardBondCounts(c.Lx,c.boundaryType);
}

MatType localHsFactor(const SpinlessTvChainUtils&c,int bond,int aux,int sigma){
    MatType generator=MatType::Zero(c.nDim,c.nDim);double lambda=std::acosh(std::exp(.5*c.V*c.dt));
    int a,b,d,e;c.aux2MajoranaIdx(aux,0,bond,a,b);c.aux2MajoranaIdx(aux,1,bond,d,e);
    DataType z(0,lambda*sigma);if(c.hsScheme==0){generator(a,b)=z;generator(b,a)=-z;
        generator(d,e)=z;generator(e,d)=-z;}else{generator(a,d)=z;generator(d,a)=-z;
        generator(b,e)=-z;generator(e,b)=z;}return exponential(generator,1.0);
}

GaussianTrialState makeTrial(const Parameters&p){
    MatType h=kineticGenerator(p.L,p.boundary,p.trial_t,p.trial_delta,p.trial_mu);
    if(p.edge_splitting!=0){SpinlessTvChainUtils coordinates(p.L,p.dt,p.V,2,p.boundary,
            p.trial_delta,p.trial_mu,p.hs);int left=coordinates.majoranaCoord2Idx(0,1);
        int right=coordinates.majoranaCoord2Idx(p.L-1,1);DataType z(0,p.edge_splitting);
        h(left,right)+=z;h(right,left)-=z;}
    GaussianTrialState trial=GaussianTrialState::fromMajoranaHamiltonian(h);
    int actual=trial.fermionParity();if(actual!=p.trial_parity)
        throw std::invalid_argument("trial parity does not match explicit trial-parity policy");
    return trial;
}

PureMirroredInitializationResult makeMirroredContour(
    const Parameters&p,const SpinlessTvChainUtils&model,std::mt19937_64&rng){
    MatType kinetic=kineticGenerator(p.L,p.boundary,p.t,p.delta,p.mu);
    MatType half=exponential(kinetic,-.5*p.dt);auto counts=bondCounts(model);
    int timeSlices=int(std::llround(p.theta/p.dt));std::uniform_int_distribution<int>bit(0,1);
    auto ketBranch=[&](){PureFastConfiguration out;int ordinal=0;
        auto push=[&](const MatType&m,int field,int slice,int bond,int aux,const std::string&label){
            out.slices.emplace_back(m,1.0,label);out.hs_fields.push_back(field);
            out.locations.push_back({PureBranch::Ket,slice,ordinal++,bond,aux});};
        for(int slice=0;slice<timeSlices;++slice){push(half,0,slice,-1,-1,"K/2");
            for(int layer=0;layer<2;++layer)for(int aux=0;aux<(layer?counts.second:counts.first);++aux){
                int sigma=bit(rng)?1:-1;push(localHsFactor(model,layer,aux,sigma),sigma,slice,layer,aux,
                    std::string("V")+std::to_string(layer)+":"+std::to_string(aux));}
            push(half,0,slice,-1,-1,"K/2");}return out;};
    // This is the sole production construction path.  The shared helper
    // appends the strict adjoint/reverse protocol already used by the Phase 2
    // noncommuting-order tests; no second independent bra RNG stream exists.
    return pureProjectorMirroredConfiguration(ketBranch(),1e-12);
}

// Retained only for small-system trajectory audits against Phase 3B.  It is
// never selected by the production default and is forbidden for L>6.
PureFastConfiguration makeSequentialAuditContour(
    const Parameters&p,const SpinlessTvChainUtils&model,std::mt19937_64&rng){
    MatType kinetic=kineticGenerator(p.L,p.boundary,p.t,p.delta,p.mu);
    MatType half=exponential(kinetic,-.5*p.dt);auto counts=bondCounts(model);
    int timeSlices=int(std::llround(p.theta/p.dt));std::uniform_int_distribution<int>bit(0,1);
    auto branch=[&](PureBranch side){PureFastConfiguration out;int ordinal=0;
        auto push=[&](const MatType&m,int field,int slice,int bond,int aux,const std::string&label){
            out.slices.emplace_back(m,1.0,label);out.hs_fields.push_back(field);
            out.locations.push_back({side,slice,ordinal++,bond,aux});};
        for(int slice=0;slice<timeSlices;++slice){push(half,0,slice,-1,-1,"K/2");
            for(int layer=0;layer<2;++layer)for(int aux=0;aux<(layer?counts.second:counts.first);++aux){
                int sigma=bit(rng)?1:-1;push(localHsFactor(model,layer,aux,sigma),sigma,slice,layer,aux,
                    std::string("V")+std::to_string(layer)+":"+std::to_string(aux));}
            push(half,0,slice,-1,-1,"K/2");}return out;};
    PureFastConfiguration ket=branch(PureBranch::Ket),bra=branch(PureBranch::Bra),full=ket;
    for(int i=int(bra.slices.size())-1;i>=0;--i){full.slices.push_back(bra.slices[i]);
        full.hs_fields.push_back(bra.hs_fields[i]);full.locations.push_back(bra.locations[i]);}
    return full;
}

std::vector<int> proposalIndices(const PureFastConfiguration&c){std::vector<int> result;
    for(int i=0;i<int(c.locations.size());++i)if(c.locations[i].aux>=0)result.push_back(i);return result;}

PureFastProposal flip(const SpinlessTvChainUtils&model,const PureFastConfiguration&c,
                      int index,double uniform){const auto&location=c.locations[index];PureFastProposal p;
    p.index=index;p.new_hs=-c.hs_fields[index];p.new_factor=localHsFactor(model,location.bond,
        location.aux,p.new_hs);p.new_eta=1.0;p.uniform=uniform;return p;}

DataType energy(const Parameters&p,const MatType&g){auto value=[&](double V,double d,double m){
    SpinlessTvChainUtils c(p.L,p.dt,V,2,p.boundary,d,m,p.hs);return c.energyFromGreensFunc(g);};
    DataType base=value(0,0,0);DataType result=p.t*base+p.delta*(value(0,1,0)-base)+
        p.mu*(value(0,0,1)-base)+(value(p.V,0,0)-base);
    if(!std::isfinite(result.real())||!std::isfinite(result.imag()))
        throw std::runtime_error("energy is nonfinite");return result;}

double parity(const SpinlessTvChainUtils&,const MatType&g){
    const PurePhysicalParityResult result=pureProjectorPhysicalParity(g);
    if(!result.ok())throw std::runtime_error("physical fermion parity is untrusted");
    return result.physical_parity;
}

std::uint64_t hashRng(const std::mt19937_64&rng){std::ostringstream state;state<<rng;
    std::uint64_t h=1469598103934665603ULL;for(unsigned char c:state.str()){h^=c;h*=1099511628211ULL;}return h;}

class RetainedCsv {public:explicit RetainedCsv(const std::string&path):path_(path){if(path!="off"&&path!="-"){
    enabled_=true;stream_.open(path);if(!stream_)throw std::runtime_error("retained CSV open failed: "+path);
    stream_<<"measurement,z2,S_pi,S_pi_dq,R_CDW,energy,parity,signed_S_pi,signed_S_pi_dq,signed_energy,signed_parity,acceptance,configuration_hash\n";
    if(!stream_)throw std::runtime_error("retained CSV header write failed: "+path);}}
    void row(int i,int z,double spi,double sdq,double r,double e,double parity,double acceptance,std::uint64_t hash){
        if(!enabled_)return;stream_<<i<<','<<z<<','<<std::setprecision(17)<<spi<<','<<sdq<<','<<r<<','<<e<<','<<parity<<','
            <<z*spi<<','<<z*sdq<<','<<z*e<<','<<z*parity<<','<<acceptance<<','<<hash<<'\n';
        if(!stream_)throw std::runtime_error("retained CSV write failed: "+path_);}
    void finish(){if(!enabled_){finished_=true;return;}stream_.flush();if(!stream_)
        throw std::runtime_error("retained CSV flush failed: "+path_);stream_.close();if(stream_.fail())
        throw std::runtime_error("retained CSV close failed: "+path_);finished_=true;}
    ~RetainedCsv(){if(enabled_&&!finished_&&stream_.is_open())stream_.close();}
private:std::string path_;std::ofstream stream_;bool enabled_=false,finished_=false;};

double binnedError(const std::vector<double>&x){if(x.size()<2)return 0;std::size_t width=std::max<std::size_t>(1,std::sqrt(x.size()));
    std::vector<double> bins;for(std::size_t begin=0;begin<x.size();begin+=width){std::size_t end=std::min(x.size(),begin+width);
        double sum=0;for(std::size_t i=begin;i<end;++i)sum+=x[i];bins.push_back(sum/double(end-begin));}
    if(bins.size()<2)return 0;double mean=0;for(double v:bins)mean+=v;mean/=bins.size();double variance=0;
    for(double v:bins)variance+=(v-mean)*(v-mean);return std::sqrt(variance/(bins.size()*(bins.size()-1)));}

void nullable(std::ostream&o,double value,bool resolved){if(resolved&&std::isfinite(value))o<<std::setprecision(17)<<value;else o<<"null";}

} // namespace

int main(int argc,char**argv){try{
    const auto runStarted=std::chrono::steady_clock::now();
    Parameters p=parse(argc,argv);GaussianTrialState trial=makeTrial(p);
    const PurePhysicalParityResult trialParity=pureProjectorPhysicalParity(trial.G_T);
    if(!trialParity.ok())throw std::runtime_error("trial physical parity diagnostic failed");
    SpinlessTvChainUtils model(p.L,p.dt,p.V,2,p.boundary,p.delta,p.mu,p.hs);
    std::mt19937_64 rng(p.seed);PureFastConfiguration initial;
    PureFastInitializationPolicy initializationPolicy=
        PureFastInitializationPolicy::MirroredTheoremZ2Plus;
    if(p.initialization_policy=="mirrored-theorem"){
        PureMirroredInitializationResult initialized=makeMirroredContour(p,model,rng);
        if(!initialized.ok())throw std::runtime_error("mirrored production initializer failed: "+
            initialized.message);
        initial=std::move(initialized.configuration);
    }else{
        initial=makeSequentialAuditContour(p,model,rng);
        initializationPolicy=PureFastInitializationPolicy::SequentialAudit;
    }
    std::vector<int> indices=proposalIndices(initial);if(indices.empty())throw std::runtime_error("contour has no HS proposals");
    PureFastOptions options;options.weight_mode=PureProjectorWeightMode::RealZ2;
    options.read_only_audit_interval=p.audit_interval;
    PureFastRunMode runMode=p.walker_mode=="audit-lockstep"?PureFastRunMode::AuditLockstep:
        PureFastRunMode::FastStrict;
    PureProjectorFastWalker walker(trial,std::move(initial),p.block,runMode,options,
        initializationPolicy);
    RetainedCsv retained(p.retained);std::uniform_real_distribution<double>uniform(0,1);
    long long accepted=0,attempted=0;std::size_t cursor=0;int direction=1,sweepDirection=1;
    auto proposalStep=[&](int index){
        double u=uniform(rng);PureFastProposalResult result=
            walker.propose(flip(model,walker.configuration(),index,u));++attempted;accepted+=result.accepted;
        if(result.terminated||!result.ratio.ok())throw std::runtime_error("proposal failed closed");};
    auto legacyProposalStep=[&](){indices=proposalIndices(walker.configuration());
        int index=indices[cursor];
        if(indices.size()>1){if(direction>0&&cursor+1==indices.size())direction=-1;
            else if(direction<0&&cursor==0)direction=1;cursor=std::size_t(int(cursor)+direction);}
        proposalStep(index);};
    auto completeSweep=[&](){indices=proposalIndices(walker.configuration());
        if(sweepDirection>0){for(int index:indices)proposalStep(index);}
        else{for(auto it=indices.rbegin();it!=indices.rend();++it)proposalStep(*it);}
        sweepDirection=-sweepDirection;};
    auto advance=[&](int count){for(int i=0;i<count;++i){
        if(p.run_units=="sweeps")completeSweep();else legacyProposalStep();}};
    advance(p.burn);const long long burnProposals=attempted;
    double signSum=0,spiNum=0,sdqNum=0,energyNum=0,energyImagNum=0,energyImagMax=0,parityNum=0;std::vector<double> signs;
    for(int measurement=0;measurement<p.measurements;++measurement){advance(p.measurement_stride);auto green=walker.measurementGreen();
        if(!green.ok())throw std::runtime_error("center Green rebuild failed");int z=walker.z2Sign();
        double spi=pureProjectorStructureFactor(green.green,p.L,kPi).real();
        double sdq=pureProjectorStructureFactor(green.green,p.L,kPi-2*kPi/p.L).real();
        double r=std::abs(spi)>1e-15?1-sdq/spi:std::numeric_limits<double>::quiet_NaN();
        DataType complexEnergy=energy(p,green.green);double e=complexEnergy.real(),fparity=parity(model,green.green);
        energyImagMax=std::max(energyImagMax,std::abs(complexEnergy.imag()));signs.push_back(z);signSum+=z;
        spiNum+=z*spi;sdqNum+=z*sdq;energyNum+=z*e;energyImagNum+=z*complexEnergy.imag();parityNum+=z*fparity;
        retained.row(measurement,z,spi,sdq,r,e,fparity,double(accepted)/attempted,walker.configurationHash());}
    if(const char*forced=std::getenv("PFQMC_TEST_FORCE_ZERO_AVERAGE_SIGN"))if(std::string(forced)=="1")signSum=0;
    retained.finish(); // completion is forbidden before the retained stream passes write/flush/close checks.
    const bool resolved=std::abs(signSum)>1e-12;double spi=resolved?spiNum/signSum:0;
    double sdq=resolved?sdqNum/signSum:0,rcdw=resolved&&std::abs(spi)>1e-15?1-sdq/spi:0;
    const PureFastDiagnostics&d=walker.diagnostics();
    std::cout<<std::setprecision(17)<<"{\"status\":\"complete\",\"projector_type\":\"pure_state\","
        "\"condition_aware_ratio\":false,\"left_recovery\":false,\"mutating_raw_checkpoint\":false,"
        "\"mutating_mp_checkpoint\":false,\"z2_oracle_correction_count\":0,\"observable_status\":\""
        <<(resolved?"resolved":"average_sign_too_small")<<"\",\"average_z2\":"<<signSum/p.measurements
        <<",\"average_z2_bin_error\":"<<binnedError(signs)<<",\"S_pi\":";nullable(std::cout,spi,resolved);
    std::cout<<",\"S_pi_dq\":";nullable(std::cout,sdq,resolved);std::cout<<",\"R_CDW\":";nullable(std::cout,rcdw,resolved);
    std::cout<<",\"energy\":";nullable(std::cout,resolved?energyNum/signSum:0,resolved);
    std::cout<<",\"energy_imaginary\":";nullable(std::cout,resolved?energyImagNum/signSum:0,resolved);
    std::cout<<",\"fermion_parity\":";nullable(std::cout,resolved?parityNum/signSum:0,resolved);
    std::cout<<",\"signed_S_pi_numerator\":"<<spiNum<<",\"signed_S_pi_dq_numerator\":"<<sdqNum
        <<",\"signed_energy_numerator\":"<<energyNum<<",\"signed_parity_numerator\":"<<parityNum
        <<",\"energy_imaginary_signed_numerator\":"<<energyImagNum
        <<",\"energy_imaginary_estimator_max\":"<<energyImagMax
        <<",\"sign_denominator\":"<<signSum<<",\"acceptance\":"<<double(accepted)/attempted
        <<",\"L\":"<<p.L<<",\"V\":"<<p.V<<",\"t\":"<<p.t<<",\"delta\":"<<p.delta
        <<",\"mu\":"<<p.mu<<",\"theta\":"<<p.theta<<",\"dt\":"<<p.dt<<",\"boundary\":\""
        <<(p.boundary?"obc":"pbc")<<"\",\"hs_scheme\":\"hs"<<p.hs<<"\",\"trial_t\":"<<p.trial_t
        <<",\"trial_delta\":"<<p.trial_delta<<",\"trial_mu\":"<<p.trial_mu
        <<",\"trial_parity\":"<<p.trial_parity<<",\"edge_splitting\":"<<p.edge_splitting
        <<",\"burn\":"<<p.burn<<",\"measurements\":"<<p.measurements<<",\"seed\":"<<p.seed
        <<",\"run_units\":\""<<p.run_units<<"\",\"hs_variable_count\":"<<indices.size()
        <<",\"burn_proposals\":"<<burnProposals
        <<",\"burn_sweep_equivalent\":"<<double(burnProposals)/indices.size()
        <<",\"measurement_stride\":"<<p.measurement_stride
        <<",\"measurement_stride_unit\":\""<<p.run_units<<"\""
        <<",\"measurement_stride_proposals\":"<<(p.run_units=="sweeps"?
            static_cast<long long>(p.measurement_stride)*static_cast<long long>(indices.size()):
            static_cast<long long>(p.measurement_stride))
        <<",\"measurement_proposals\":"<<(attempted-burnProposals)
        <<",\"measurement_sweep_equivalent\":"<<double(attempted-burnProposals)/indices.size()
        <<",\"proposal_count\":"<<attempted
        <<",\"fermion_parity_convention\":\"block_majorana_physical\""
        <<",\"trial_internal_pfaffian_sign\":"<<trialParity.internal_pfaffian_sign
        <<",\"trial_block_reordering_sign\":"<<trialParity.block_reordering_sign
        <<",\"stabilization_block\":"<<p.block<<",\"audit_interval\":"<<p.audit_interval
        <<",\"walker_mode\":\""<<p.walker_mode<<"\""
        <<",\"initialization_policy\":\""<<walker.initializationPolicy()<<"\""
        <<",\"proposal_schedule\":\""<<(p.run_units=="sweeps"?
            "complete_sweep_alternating":"legacy_forward_backward")<<"\""
        <<",\"minimum_overlap_rcond\":"<<d.minimum_overlap_rcond
        <<",\"green_fast_rebuild_relative_error_max\":"<<d.maximum_green_rebuild_error
        <<",\"ratio_reference_relative_error_max\":"<<d.maximum_ratio_reference_error
        <<",\"rebuild_count\":"<<d.rebuild_count<<",\"ratio_slow_reference_count\":"<<d.ratio_slow_reference_count
        <<",\"trust_alarm_count\":"<<d.trust_alarm_count<<",\"slow_reference_failure_count\":"<<d.slow_reference_failure_count
        <<",\"mp_same_proposal_fallback_count\":"<<d.mp_fallback_count
        <<",\"mp_same_proposal_fallback_failure_count\":"<<d.mp_fallback_failure_count
        <<",\"mp_green_double_recovery_count\":"<<d.mp_green_double_recovery_count
        <<",\"read_only_endpoint_audit_failure_count\":"<<d.read_only_endpoint_audit_failure_count
        <<",\"fast_path_seconds\":"<<d.total_fast_seconds
        <<",\"reference_seconds\":"<<d.total_reference_seconds
        <<",\"mp_fallback_seconds\":"<<d.total_mp_fallback_seconds
        <<",\"mp_fallback_average_seconds\":"<<(d.mp_fallback_count?
            d.total_mp_fallback_seconds/d.mp_fallback_count:0.0)
        <<",\"mp_canonical_input_builds\":"<<d.mp_profile.canonical_input_builds
        <<",\"mp_cache_invalidations\":"<<d.mp_profile.cache_invalidations
        <<",\"mp_operator_requests\":"<<d.mp_profile.operator_requests
        <<",\"mp_operator_cache_hits\":"<<d.mp_profile.operator_cache_hits
        <<",\"mp_operator_cache_misses\":"<<d.mp_profile.operator_cache_misses
        <<",\"mp_operator_cache_hit_rate\":"<<(d.mp_profile.operator_requests?
            double(d.mp_profile.operator_cache_hits)/d.mp_profile.operator_requests:0.0)
        <<",\"mp_sparse_apply_count\":"<<d.mp_profile.sparse_apply_count
        <<",\"mp_dense_apply_count\":"<<d.mp_profile.dense_apply_count
        <<",\"mp_canonicalization_seconds\":"<<d.mp_profile.canonicalization_seconds
        <<",\"mp_precision_160_seconds\":"<<d.mp_profile.precision_160_seconds
        <<",\"mp_precision_320_seconds\":"<<d.mp_profile.precision_320_seconds
        <<",\"mp_precision_640_seconds\":"<<d.mp_profile.precision_640_seconds
        <<",\"mp_conversion_seconds\":"<<d.mp_profile.conversion_seconds
        <<",\"mp_propagation_seconds\":"<<d.mp_profile.propagation_seconds
        <<",\"mp_thin_qr_seconds\":"<<d.mp_profile.thin_qr_seconds
        <<",\"mp_endpoint_seconds\":"<<d.mp_profile.endpoint_seconds
        <<",\"mp_local_pfaffian_seconds\":"<<d.mp_profile.local_pfaffian_seconds
        <<",\"runtime_seconds\":"<<std::chrono::duration<double>(
            std::chrono::steady_clock::now()-runStarted).count()
        <<",\"endpoint_rebuild_green_residual_max\":"<<d.maximum_endpoint_rebuild_green_residual
        <<",\"first_failure_proposal\":"<<d.first_failure_proposal<<",\"final_hs_hash\":"<<walker.configurationHash()
        <<",\"final_rng_hash\":"<<hashRng(rng)<<",\"source_commit\":\""<<jsonEscape(p.source_commit)
        <<"\",\"executable_sha256\":\""<<jsonEscape(p.executable_sha256)<<"\"}\n";
    if(!std::cout)throw std::runtime_error("completion JSON write failed");return 0;
}catch(const std::exception&e){std::cerr<<"pure_projector_driver: "<<e.what()<<'\n';return 1;}}
