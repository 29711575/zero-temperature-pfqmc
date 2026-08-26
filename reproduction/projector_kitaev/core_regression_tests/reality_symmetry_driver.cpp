#define main pfqmc_embedded_mp_oracle_main
#define PFQMC_MP_ORACLE_SMALL_ONLY
#include "mp_ratio_oracle_driver.cpp"
#undef PFQMC_MP_ORACLE_SMALL_ONLY
#undef main

int main(int argc,char **argv) try {
    if(argc!=6)throw std::runtime_error("usage: reality_symmetry_driver L V seed samples output.csv");
    const int L=std::stoi(argv[1]),seed=std::stoi(argv[3]),samples=std::stoi(argv[4]);const double V=std::stod(argv[2]);
    constexpr double dt=.1,beta=2.,theta=2.;
    SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(theta/dt)),0,1.0,0.0,0);
    rdGenerator rd(seed);Walker w(&cfg,&rd,theta,beta);PfQMC q(&w,10);
    std::vector<int> hsops;for(int i=0;i<int(q.op_array.size());++i)if(dynamic_cast<SpinlessVOperator*>(q.op_array[i]))hsops.push_back(i);
    std::ofstream out(argv[5]);out<<"L,V,seed,sample,boundary,aux,field_hash,precision_low,precision_high,single_operator_symmetry_residual,trial_operator_symmetry_residual,kinetic_half_symmetry_residual,hs_operator_symmetry_residual,cyclic_product_symmetry_residual,ratio_low_real,ratio_low_imag,ratio_high_real,ratio_high_imag,ratio_reality_residual,ratio_precision_relative_difference,r_squared_Q_relative,solve_residual,condition_estimate,finite\n";
    for(int sample=0;sample<samples;++sample){
        if(sample){q.rightSweep();q.leftSweep();}
        int boundary=hsops[(std::size_t(sample)*hsops.size())/samples];auto *op=dynamic_cast<SpinlessVOperator*>(q.op_array[boundary]);int aux=sample%op->s->size();
        auto lo=evaluateOracle<80>(cfg,q.op_array,w.ntrial,boundary,aux);auto hi=evaluateOracle<160>(cfg,q.op_array,w.ntrial,boundary,aux);
        using X=MP<160>;X::Real lor(MP<80>::str(lo.direct.r)),loi(MP<80>::str(lo.direct.i));
        X::C low(lor,loi);
        const X::Real directScale=std::max(X::Real(X::absc(hi.direct)),X::Real("1e-100"));
        const X::Real pd=X::absc(X::sub(low,hi.direct))/directScale;
        const X::Real reality=X::absc(X::C(0,hi.direct.i))/directScale;
        auto single=std::max(hi.symKinetic,hi.symHS);bool finite=std::isfinite(hi.direct.r.convert_to<double>())&&std::isfinite(hi.direct.i.convert_to<double>());
        out<<L<<','<<V<<','<<seed<<','<<sample<<','<<boundary<<','<<aux<<','<<fieldHash(q.op_array)<<",80,160,"
           <<X::str(single)<<','<<X::str(hi.symKinetic)<<','<<X::str(hi.symKinetic)<<','<<X::str(hi.symHS)<<','<<X::str(hi.symProduct)<<','
           <<MP<80>::str(lo.direct.r)<<','<<MP<80>::str(lo.direct.i)<<','<<X::str(hi.direct.r)<<','<<X::str(hi.direct.i)<<','<<X::str(reality)<<','<<X::str(pd)<<','<<X::str(hi.r2q)<<','<<X::str(hi.solveResidual)<<','<<X::str(hi.condition)<<','<<(finite?1:0)<<'\n';
    }
    std::cout<<"{\"status\":\"complete\",\"samples\":"<<samples<<"}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
