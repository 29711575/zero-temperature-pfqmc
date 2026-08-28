#include "test_common.h"

#include <iostream>

int main(int argc,char **argv) try {
    if (argc!=11) throw std::runtime_error(
        "usage: left_green_recovery_driver L bc V seed theta beta dt sweeps threshold output.csv");
    const int L=std::stoi(argv[1]),bc=std::stoi(argv[2]),seed=std::stoi(argv[4]);
    const double V=std::stod(argv[3]),theta=std::stod(argv[5]);
    const double beta=std::stod(argv[6]),dt=std::stod(argv[7]);
    const int sweeps=std::stoi(argv[8]); const double threshold=std::stod(argv[9]);
    const int ns=2*int(std::llround(theta/dt));
    SpinlessTvChainUtils cfg(L,dt,V,ns,bc,1.0,0.0,0);
    rdGenerator rng(seed); CoreTestWalker walker(&cfg,&rng,theta,beta); PfQMC q(&walker,10);
    q.configureLeftSweepGreenRecovery(true,threshold,threshold);
    std::ofstream out(argv[10]);
    out<<"sweep,field_hash,relative_green_error,structure_residual,finite,"
          "rank_recoveries,propagation_recoveries\n"<<std::setprecision(17);
    double maxError=0.0; long long nonfinite=0;
    for (int sweep=0;sweep<sweeps;++sweep) {
        q.rightSweep(); q.leftSweep(); MatType full;
        q.rebuildGreenFromFullContourAtBoundary(0,full);
        const bool finite=matrixFinite(q.g)&&matrixFinite(full);
        const double error=finite?relativeError(q.g,full):
            std::numeric_limits<double>::infinity();
        maxError=std::max(maxError,error); nonfinite+=!finite;
        out<<sweep<<','<<hsHash(q.op_array)<<','<<error<<','
           <<structureResidual(q.g)<<','<<(finite?1:0)<<','
           <<q.left_recovery_rank_update_count<<','
           <<q.left_recovery_propagation_count<<'\n';
    }
    std::cout<<"{\"status\":\"complete\",\"nonfinite\":"<<nonfinite
             <<",\"max_error\":"<<std::setprecision(17)<<maxError
             <<",\"rank_recoveries\":"<<q.left_recovery_rank_update_count
             <<",\"propagation_recoveries\":"
             <<q.left_recovery_propagation_count<<"}\n";
    return nonfinite?3:0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
