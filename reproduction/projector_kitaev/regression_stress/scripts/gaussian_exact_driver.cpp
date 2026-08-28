#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"

class GaussianWalker : public Spinless_tV {
public:
    int center=-1, ntrial=0, nphys=0;
    GaussianWalker(const SpinlessTvChainUtils *c, rdGenerator *r, double theta, double beta) {
        build_projector_static_contour(*this,c,r,theta,beta,center,ntrial,nphys);
    }
};

int main(int argc,char **argv) try {
    if(argc!=9) throw std::runtime_error("usage: gaussian_exact_driver L theta beta dt delta mu boundary hs");
    const int L=std::stoi(argv[1]), boundary=std::stoi(argv[7]), hs=std::stoi(argv[8]);
    const double theta=std::stod(argv[2]), beta=std::stod(argv[3]), dt=std::stod(argv[4]), delta=std::stod(argv[5]), mu=std::stod(argv[6]);
    const int ns=2*int(std::llround(theta/dt));
    SpinlessTvChainUtils cfg(L,dt,0.0,ns,boundary,delta,mu,hs);
    rdGenerator rng(700001); GaussianWalker walker(&cfg,&rng,theta,beta); PfQMC q(&walker,10);
    MatType g; q.rebuildGreenFromFullContourAtBoundary(walker.center,g);
    const double spi=-cfg.StructureFactorCDW(g).real();
    const double sdq=-cfg.StructureFactorCDWOffset(g).real();
    std::cout<<std::setprecision(17)<<"{\"method\":\"V0_gaussian_full_contour\",\"L\":"<<L<<",\"theta\":"<<theta<<",\"beta_trial\":"<<beta<<",\"dt\":"<<dt<<",\"delta\":"<<delta<<",\"mu\":"<<mu<<",\"boundary\":"<<boundary<<",\"hs_scheme\":"<<hs<<",\"center_boundary\":"<<walker.center<<",\"S_pi\":"<<spi<<",\"S_pi_dq\":"<<sdq<<",\"R_cdw\":"<<(1-sdq/spi)<<",\"center_green_norm\":"<<g.norm()<<",\"center_green_max_imag\":"<<g.imag().cwiseAbs().maxCoeff()<<",\"raw_sign_real\":"<<q.getSignRaw().real()<<",\"raw_sign_imag\":"<<q.getSignRaw().imag()<<"}\n";
    return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
