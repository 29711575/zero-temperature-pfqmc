#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"

class EnumWalker : public Spinless_tV {
public:
    int center=-1, ntrial=0, nphys=0;
    EnumWalker(const SpinlessTvChainUtils *c, rdGenerator *r, double theta, double beta) {
        build_projector_static_contour(*this,c,r,theta,beta,center,ntrial,nphys);
    }
};

struct FieldRef { SpinlessVOperator *op; int index; };

static MatType cyclic_product(const EnumWalker &w,int boundary) {
    MatType product=MatType::Identity(w.nDim,w.nDim),tmp;
    const int n=int(w.op_array.size());
    for(int off=0;off<n;++off){w.op_array[(boundary+off)%n]->left_multiply(product,tmp);product.swap(tmp);}
    return product;
}

int main(int argc,char **argv) try {
    if(argc!=12) throw std::runtime_error("usage: tiny_enumeration_driver L theta beta dt V delta mu boundary hs random_checks random_csv");
    const int L=std::stoi(argv[1]), boundary=std::stoi(argv[8]), hs=std::stoi(argv[9]), random_checks=std::stoi(argv[10]);
    const double theta=std::stod(argv[2]), beta=std::stod(argv[3]), dt=std::stod(argv[4]), V=std::stod(argv[5]), delta=std::stod(argv[6]), mu=std::stod(argv[7]);
    const int ns=2*int(std::llround(theta/dt)); rdGenerator rng(710001);
    SpinlessTvChainUtils cfg(L,dt,V,ns,boundary,delta,mu,hs); EnumWalker w(&cfg,&rng,theta,beta);
    std::vector<FieldRef> fields; std::vector<SpinlessVOperator*> interaction_ops;
    for(Operator *base:w.op_array) if(auto *op=dynamic_cast<SpinlessVOperator*>(base)){
        interaction_ops.push_back(op); for(int i=0;i<op->s->size();++i) fields.push_back({op,i});
    }
    if(fields.size()>24) throw std::runtime_error("enumeration limited to 24 HS fields");
    const unsigned long long configs=1ULL<<fields.size();
    std::vector<double> logw(configs),spi(configs),sdq(configs); std::vector<DataType> phase(configs);
    std::ofstream random(argv[11]); random<<"configuration,raw_weight_real,raw_weight_imag,raw_sign_real,raw_sign_imag,center_green_relative_error\n"<<std::setprecision(17);
    double maxlog=-std::numeric_limits<double>::infinity(),max_gerr=0,max_phase_norm_dev=0;
    for(unsigned long long mask=0;mask<configs;++mask){
        for(size_t j=0;j<fields.size();++j)(*fields[j].op->s)(fields[j].index)=((mask>>j)&1)?1:-1;
        for(auto *op:interaction_ops){op->B.setIdentity();cfg.InteractionBGenerator(op->B,*op->s,op->bondType);op->reCalcInv();}
        PfQMC q(&w,1); MatType g;q.rebuildGreenFromFullContourAtBoundary(w.center,g);
        phase[mask]=q.getSignRaw();max_phase_norm_dev=std::max(max_phase_norm_dev,std::abs(std::abs(phase[mask])-1.0));
        MatType product0=cyclic_product(w,0); Eigen::FullPivLU<MatType> lu(MatType::Identity(w.nDim,w.nDim)+product0);
        const DataType det=lu.determinant();logw[mask]=0.5*std::log(std::max(std::abs(det),std::numeric_limits<double>::min()));maxlog=std::max(maxlog,logw[mask]);
        spi[mask]=-cfg.StructureFactorCDW(g).real();sdq[mask]=-cfg.StructureFactorCDWOffset(g).real();
        if(random_checks>0 && (mask<unsigned(random_checks) || mask+unsigned(random_checks)>=configs)){
            MatType pc=cyclic_product(w,w.center);MatType direct=2.0*(MatType::Identity(w.nDim,w.nDim)+pc).inverse();
            const double err=(g-direct).norm()/std::max(direct.norm(),std::numeric_limits<double>::min());max_gerr=std::max(max_gerr,err);
            const DataType raw=phase[mask]*std::exp(logw[mask]);random<<mask<<','<<raw.real()<<','<<raw.imag()<<','<<phase[mask].real()<<','<<phase[mask].imag()<<','<<err<<'\n';
        }
    }
    DataType den=0,num_pi=0,num_dq=0; double absden=0;
    for(unsigned long long mask=0;mask<configs;++mask){const DataType weight=phase[mask]*std::exp(logw[mask]-maxlog);den+=weight;absden+=std::abs(weight);num_pi+=weight*spi[mask];num_dq+=weight*sdq[mask];}
    const DataType p=num_pi/den,d=num_dq/den;
    std::cout<<std::setprecision(17)<<"{\"method\":\"complete_HS_exact_enumeration\",\"L\":"<<L<<",\"theta\":"<<theta<<",\"beta_trial\":"<<beta<<",\"dt\":"<<dt<<",\"V\":"<<V<<",\"delta\":"<<delta<<",\"mu\":"<<mu<<",\"boundary\":"<<boundary<<",\"hs_scheme\":"<<hs<<",\"hs_fields\":"<<fields.size()<<",\"configurations\":"<<configs<<",\"S_pi\":"<<p.real()<<",\"S_pi_imag\":"<<p.imag()<<",\"S_pi_dq\":"<<d.real()<<",\"S_pi_dq_imag\":"<<d.imag()<<",\"R_cdw\":"<<(1.0-d/p).real()<<",\"average_phase_real\":"<<(den/absden).real()<<",\"average_phase_imag\":"<<(den/absden).imag()<<",\"max_raw_phase_norm_deviation\":"<<max_phase_norm_dev<<",\"max_center_green_direct_relative_error\":"<<max_gerr<<"}\n";
    return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
