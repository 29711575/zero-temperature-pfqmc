#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"

class Walker : public Spinless_tV {
public:
    int center=-1,ntrial=0,nphys=0;
    Walker(const SpinlessTvChainUtils*c,rdGenerator*r,double theta,double beta) {
        build_projector_static_contour(*this,c,r,theta,beta,center,ntrial,nphys);
    }
};

static std::vector<int> fields(const Spinless_tV&w) {
    std::vector<int> out;
    for (Operator* base:w.op_array) if (auto*o=dynamic_cast<SpinlessVOperator*>(base))
        for(int i=0;i<o->s->size();++i) out.push_back((*o->s)(i));
    return out;
}
static std::uint64_t hash_fields(const std::vector<int>&v) {
    std::uint64_t h=1469598103934665603ULL;
    for(int x:v){h^=std::uint64_t(std::int64_t(x));h*=1099511628211ULL;}
    return h;
}
static void cj(std::ostream&o,DataType z){o<<'['<<z.real()<<','<<z.imag()<<']';}
static void mj(std::ostream&o,const MatType&m){o<<'[';for(int i=0;i<m.rows();++i)for(int j=0;j<m.cols();++j){if(i||j)o<<',';cj(o,m(i,j));}o<<']';}
static double maxabs(const MatType&a,const MatType&b){return (a-b).cwiseAbs().maxCoeff();}
static double relfro(const MatType&a,const MatType&b){return (a-b).norm()/std::max(b.norm(),std::numeric_limits<double>::min());}

static void phase_row(std::ofstream&f,int sample,const char*phase,int boundary,const Spinless_tV&w,
                      PfQMC&q,const MatType&fast) {
    const auto hs=fields(w); MatType full; q.rebuildGreenFromFullContourAtBoundary(boundary,full);
    const DataType direct=q.getSignRaw();
    f<<sample<<','<<phase<<','<<boundary<<','<<hash_fields(hs)<<','
     <<q.sign.real()<<','<<q.sign.imag()<<','<<direct.real()<<','<<direct.imag()<<','<<(q.sign.real()>=0)<<','<<(direct.real()>=0)<<','
     <<maxabs(fast,full)<<','<<relfro(fast,full)<<'\n';
}

static void snapshot(const std::string&path,int sample,int boundary,const Spinless_tV&w,
                     const Walker&meta,const PfQMC&q,const MatType&fast,const MatType&full,
                     DataType direct,double L,double V,double theta,double beta,double dt,
                     int boundary_type,int hs,int seed) {
    std::ofstream f(path);f<<std::setprecision(17);const auto hv=fields(w);
    f<<"{\"sample\":"<<sample<<",\"L\":"<<int(L)<<",\"V\":"<<V<<",\"theta\":"<<theta
     <<",\"beta_trial\":"<<beta<<",\"dt\":"<<dt<<",\"boundary_type\":"<<boundary_type
     <<",\"hs_scheme\":"<<hs<<",\"seed\":"<<seed<<",\"phase\":\"after_right_sweep\""
     <<",\"exact_boundary\":"<<boundary<<",\"hs_hash\":"<<hash_fields(hv)<<",\"tracked_sign\":";cj(f,q.sign);
    f<<",\"direct_sign\":";cj(f,direct);f<<",\"fast_G\":";mj(f,fast);f<<",\"full_G\":";mj(f,full);
    f<<",\"trial_slices\":"<<meta.ntrial<<",\"physical_slices\":"<<meta.nphys<<",\"operator_count\":"<<w.op_array.size()<<",\"operators\":[";
    bool first=true;for(int j=0;j<(int)w.op_array.size();++j)if(auto*o=dynamic_cast<SpinlessVOperator*>(w.op_array[j])){
      if(!first)f<<',';first=false;f<<"{\"index\":"<<j<<",\"local_V\":"<<o->localV<<",\"bond\":"<<o->bondType<<",\"s\":[";
      for(int k=0;k<o->s->size();++k){if(k)f<<',';f<<(*o->s)(k);}f<<"]}";
    }f<<"]}\n";
}

int main(int argc,char**argv)try{
    if(argc!=12)throw std::runtime_error("usage: driver L V theta beta dt boundary hs seed burn samples outdir");
    int L=std::stoi(argv[1]),boundary_type=std::stoi(argv[6]),hs=std::stoi(argv[7]),seed=std::stoi(argv[8]),burn=std::stoi(argv[9]),samples=std::stoi(argv[10]);
    double V=std::stod(argv[2]),theta=std::stod(argv[3]),beta=std::stod(argv[4]),dt=std::stod(argv[5]);std::string out=argv[11];
    SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(theta/dt)),boundary_type,1,0,hs);rdGenerator rd(seed);Walker w(&cfg,&rd,theta,beta);PfQMC q(&w,10);
    for(int i=0;i<burn;++i){q.rightSweep();q.leftSweep();}
    std::ofstream phases(out+"/phases.csv");phases<<"sample,phase,boundary,hs_hash,tracked_sign_real,tracked_sign_imag,direct_sign_real,direct_sign_imag,tracked_pm,direct_pm,green_max_abs_error,green_fro_relative_error\n"<<std::setprecision(17);
    for(int k=0;k<samples;++k){
      phase_row(phases,k,"before_right_sweep",0,w,q,q.g);
      q.rightSweep();MatType fast=q.g,full;q.rebuildGreenFromFullContourAtBoundary(0,full);DataType direct=q.getSignRaw();
      phase_row(phases,k,"after_right_sweep",0,w,q,fast);
      snapshot(out+"/snapshot_"+std::to_string(k)+".json",k,0,w,w,q,fast,full,direct,L,V,theta,beta,dt,boundary_type,hs,seed);
      q.leftSweep();phase_row(phases,k,"after_left_sweep",0,w,q,q.g);
    }
    std::cout<<"{\"L\":"<<L<<",\"V\":"<<V<<",\"theta\":"<<theta<<",\"seed\":"<<seed<<",\"samples\":"<<samples<<"}\n";
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
