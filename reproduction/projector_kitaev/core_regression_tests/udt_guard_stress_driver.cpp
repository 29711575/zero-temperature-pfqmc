#define main udt_orthogonality_embedded_main
#include "udt_orthogonality_driver.cpp"
#undef main

#include <iostream>

namespace {
void runCase(int n,int range,int seed,bool expectFail,std::ofstream &out) {
    std::mt19937 rng(seed); UDT right=makeCanonical(n,range,rng,false);
    UDT left=makeCanonical(n,range,rng,true,&right);
    resetScaleSafeQRGuardDiagnostics(); std::string outcome="success",message;
    double orthogonality=std::numeric_limits<double>::quiet_NaN();
    try {
        UDT product=left*right; orthogonality=metric(product.U).adj;
        if (!std::isfinite(orthogonality) || orthogonality>1e-6)
            throw std::runtime_error("finite nonunitary U escaped guard");
    } catch(const std::exception &e) { outcome="fail_closed"; message=e.what(); }
    const auto diagnostics=scaleSafeQRGuardDiagnostics();
    const bool passed=expectFail ? outcome=="fail_closed" : outcome=="success";
    out<<n<<','<<range<<','<<seed<<','<<(expectFail?1:0)<<','<<outcome<<",\""
       <<message<<"\","<<std::setprecision(17)<<orthogonality<<','
       <<diagnostics.trigger_count<<','<<diagnostics.max_lost_bits<<','
       <<diagnostics.min_guard_margin<<','<<(passed?1:0)<<'\n';
    if (!passed) throw std::runtime_error("unexpected guard outcome");
}
}

int main(int argc,char **argv) try {
    if (argc!=2) throw std::runtime_error("usage: udt_guard_stress_driver output.csv");
    std::ofstream out(argv[1]);
    out<<"size,exponent_range,seed,expect_fail_closed,outcome,message,"
          "output_orthogonality,trigger_count,max_lost_bits,min_guard_margin,pass\n";
    runCase(12,20,900012,false,out);
    runCase(24,40,900024,true,out);
    runCase(12,500,900012,true,out);
    runCase(52,1500,903520,true,out);
    runCase(52,2000,903520,true,out);
    std::cout<<"{\"status\":\"complete\",\"cases\":5}\n";
    return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
