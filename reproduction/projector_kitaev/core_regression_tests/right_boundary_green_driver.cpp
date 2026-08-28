#include "test_common.h"

#include <iostream>

namespace {
void advanceToBoundary(PfQMC &q,int requestedBoundary) {
    MatType tmp=MatType::Identity(q.nDim,q.nDim);
    MatType segment=MatType::Identity(q.nDim,q.nDim);
    int currentSegment=0;
    for (int l=0;l<requestedBoundary;++l) {
        const DataType signCurrent=q.op_array[l]->update(q.g);
        q.sign*=signCurrent;
        q.op_array[l]->left_multiply(segment,tmp); segment.swap(tmp);
        if (q.need_stabilization[(l+1)%q.op_length]) {
            if (currentSegment==0) q.udtR[currentSegment]=UDT(segment);
            else q.udtR[currentSegment]=segment*q.udtR[currentSegment-1];
            segment.setIdentity();
            if (currentSegment==q.checkpoints-1)
                q.udtR[currentSegment].onePlusInv(q.g);
            else q.g=onePlusInv(q.udtL[currentSegment+1],q.udtR[currentSegment]);
            ++currentSegment;
        } else {
            q.op_array[l]->left_propagate(q.g,tmp);
        }
    }
}
}

int main(int argc,char **argv) try {
    if (argc!=9) throw std::runtime_error(
        "usage: right_boundary_green_driver L bc V seed theta beta dt output.csv");
    const int L=std::stoi(argv[1]),bc=std::stoi(argv[2]),seed=std::stoi(argv[4]);
    const double V=std::stod(argv[3]),theta=std::stod(argv[5]);
    const double beta=std::stod(argv[6]),dt=std::stod(argv[7]);
    const int ns=2*int(std::llround(theta/dt));
    SpinlessTvChainUtils probeCfg(L,dt,V,ns,bc,1.0,0.0,0);
    rdGenerator probeRng(seed); CoreTestWalker probe(&probeCfg,&probeRng,theta,beta);
    const int boundaries=int(probe.op_array.size());
    std::ofstream out(argv[8]);
    out<<"L,bc,V,seed,boundary,boundary_type,relative_green_error,"
          "fast_structure_residual,full_structure_residual,finite,max_abs_green\n"
       <<std::setprecision(17);
    double maxError=0.0; long long nonfinite=0;
    for (int boundary=1; boundary<=boundaries; ++boundary) {
        SpinlessTvChainUtils cfg(L,dt,V,ns,bc,1.0,0.0,0);
        rdGenerator rng(seed); CoreTestWalker walker(&cfg,&rng,theta,beta);
        PfQMC q(&walker,10); MatType full;
        advanceToBoundary(q,boundary); const MatType fast=q.g;
        const int cyclicBoundary=boundary%boundaries;
        q.rebuildGreenFromFullContourAtBoundary(cyclicBoundary,full);
        const bool finite=matrixFinite(fast)&&matrixFinite(full);
        const double error=finite?relativeError(fast,full):
            std::numeric_limits<double>::infinity();
        maxError=std::max(maxError,error); nonfinite+=!finite;
        out<<L<<','<<bc<<','<<V<<','<<seed<<','<<cyclicBoundary<<','
           <<boundaryKind(walker,cyclicBoundary)<<','<<error<<','
           <<structureResidual(fast)<<','<<structureResidual(full)<<','
           <<(finite?1:0)<<','<<maxAbs(fast)<<'\n';
    }
    std::cout<<"{\"status\":\"complete\",\"boundaries\":"<<boundaries
             <<",\"nonfinite\":"<<nonfinite<<",\"max_error\":"
             <<std::setprecision(17)<<maxError<<"}\n";
    return nonfinite?3:0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
