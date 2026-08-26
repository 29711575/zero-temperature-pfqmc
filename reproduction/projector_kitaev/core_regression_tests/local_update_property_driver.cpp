#include "test_common.h"

#include <iostream>
#include <random>

namespace {

MatType cyclicProduct(const std::vector<Operator*> &ops,int boundary) {
    const int n=int(ops.size()); MatType p=MatType::Identity(ops[0]->getAuxField()?0:1,1);
    MatType probe; ops[0]->getGreensMat(probe); p=MatType::Identity(probe.rows(),probe.rows());
    MatType tmp;
    for(int off=0;off<n;++off){ops[(boundary+off)%n]->left_multiply(p,tmp);p.swap(tmp);}
    return p;
}

void rank2(MatType &g,int a,int b,DataType alpha) {
    cVecType x1=-g.col(a),x2=-g.col(b); x1(a)+=2.0; x2(b)+=2.0;
    g += alpha*(x1*x2.transpose()-x2*x1.transpose());
}

DataType sequentialRatio(const SpinlessTvChainUtils &cfg,const SpinlessVOperator &op,
                         int aux,const MatType &green) {
    int a,b,c,d; cfg.aux2MajoranaIdx(aux,0,op.bondType,a,b); cfg.aux2MajoranaIdx(aux,1,op.bondType,c,d);
    const int sigma=(*op.s)(aux); const DataType I(0,1); MatType g=green;
    if(op.hsScheme==0){
        DataType d1=1.0-I*op.thlV*double(sigma)*g(a,b); rank2(g,a,b,I*double(sigma)*op.thlV/d1);
        DataType d2=1.0-I*op.thlV*double(sigma)*g(c,d); return (d1*op.chlV)*(d2*op.chlV);
    }
    DataType d1=1.0-I*op.thlV*double(sigma)*g(a,c); rank2(g,a,c,I*double(sigma)*op.thlV/d1);
    DataType d2=1.0+I*op.thlV*double(sigma)*g(b,d); return (d1*op.chlV)*(d2*op.chlV);
}

double scalarRel(DataType a,DataType b) {
    return std::abs(a-b)/std::max(std::abs(b),std::numeric_limits<double>::min());
}

}

int main(int argc,char **argv) try {
    if(argc!=10)throw std::runtime_error(
        "usage: local_update_property_driver L boundary hs V seed flips theta beta output.csv");
    const int L=std::stoi(argv[1]),bc=std::stoi(argv[2]),hs=std::stoi(argv[3]),seed=std::stoi(argv[5]);
    const double V=std::stod(argv[4]),theta=std::stod(argv[7]),beta=std::stod(argv[8]);
    const int flips=std::stoi(argv[6]); constexpr double dt=.1;
    SpinlessTvChainUtils cfg(L,dt,V,2*int(std::llround(theta/dt)),bc,1.0,0.0,hs);
    rdGenerator rng(seed); CoreTestWalker w(&cfg,&rng,theta,beta); PfQMC q(&w,10);
    std::vector<int> boundaries; for(int i=0;i<int(q.op_array.size());++i)
        if(dynamic_cast<SpinlessVOperator*>(q.op_array[i]))boundaries.push_back(i);
    if(boundaries.empty())throw std::runtime_error("no HS operators");
    std::mt19937 pick(seed+99173); std::uniform_int_distribution<int> bop(0,int(boundaries.size())-1);
    std::ofstream out(argv[9]); if(!out)throw std::runtime_error("cannot open output");
    out<<"L,boundary_condition,hs_scheme,V,seed,flip,boundary,aux,field_hash_before,core_rcond,direct_real,direct_imag,sequential_real,sequential_imag,full_rebuild_ratio_real,full_rebuild_ratio_imag,direct_sequential_relative,direct_full_relative,direct_squared_Q_relative,rank_update_green_relative,imaginary_residual,finite,conditioning_class\n"<<std::setprecision(17);
    double maxRatio=0,maxGreen=0,maxImag=0; long long bad=0;
    for(int k=0;k<flips;++k){
        int boundary=boundaries[bop(pick)]; auto *op=dynamic_cast<SpinlessVOperator*>(q.op_array[boundary]);
        std::uniform_int_distribution<int> apick(0,op->s->size()-1); int aux=apick(pick);
        MatType full; q.rebuildGreenFromFullContourAtBoundary(boundary,full); q.g=full;
        MatType oldProduct=cyclicProduct(q.op_array,boundary);
        MatType oldCore=MatType::Identity(q.nDim,q.nDim)+oldProduct;
        DataType detOld=oldCore.fullPivLu().determinant();
        const double rcond=q.fullContourCoreConditionAtBoundary(boundary);
        double unused=0; if(!op->prepareSingleFlip(q.g,aux,&unused))throw std::runtime_error("prepare failed");
        const DataType direct=op->preparedRatio(),sequential=sequentialRatio(cfg,*op,aux,q.g);
        MatType independentFull; q.rebuildGreenFromFullContourAtBoundary(boundary,independentFull);
        const DataType fullRebuildRatio=op->recomputePreparedRatio(independentFull);
        const std::uint64_t hashBefore=hsHash(q.op_array);
        op->finishSingleFlip(q.g,true,true); MatType rankUpdated=q.g;
        MatType acceptedFull; q.rebuildGreenFromFullContourAtBoundary(boundary,acceptedFull);
        MatType newProduct=cyclicProduct(q.op_array,boundary);
        DataType detNew=(MatType::Identity(q.nDim,q.nDim)+newProduct).fullPivLu().determinant();
        DataType Q=detNew/detOld;
        const double ds=scalarRel(direct,sequential),df=scalarRel(direct,fullRebuildRatio);
        const double rq=scalarRel(direct*direct,Q),ge=relativeError(rankUpdated,acceptedFull);
        const double imag=std::max({std::abs(direct.imag()),std::abs(sequential.imag()),std::abs(fullRebuildRatio.imag())});
        const bool finite=std::isfinite(ds)&&std::isfinite(df)&&std::isfinite(rq)&&std::isfinite(ge)&&matrixFinite(rankUpdated)&&matrixFinite(acceptedFull);
        maxRatio=std::max({maxRatio,ds,df,rq});maxGreen=std::max(maxGreen,ge);maxImag=std::max(maxImag,imag);
        // The explicit dense determinant is retained as a conditioning
        // diagnostic only; the stabilized full rebuild is the trusted gate.
        // Q uses an unstabilized dense determinant and remains a conditioning
        // diagnostic. The integration gate is the independent stabilized
        // ratio plus the accepted-configuration stabilized Green.
        const bool anomaly=!finite||ds>1e-8||df>1e-8||ge>1e-6;
        if(anomaly&&bad++==0){std::ostringstream h;h<<"L="<<L<<" bc="<<bc<<" hs="<<hs<<" V="<<V<<" seed="<<seed<<" flip="<<k<<" boundary="<<boundary<<" aux="<<aux;dumpConfiguration(q.op_array,"local_update_min_repro.txt",h.str());}
        out<<L<<','<<bc<<','<<hs<<','<<V<<','<<seed<<','<<k<<','<<boundary<<','<<aux<<','<<hashBefore<<','<<rcond<<','
           <<direct.real()<<','<<direct.imag()<<','<<sequential.real()<<','<<sequential.imag()<<','<<fullRebuildRatio.real()<<','<<fullRebuildRatio.imag()<<','
           <<ds<<','<<df<<','<<rq<<','<<ge<<','<<imag<<','<<(finite?1:0)<<','<<(rcond<1e-12?"poor":"clean")<<'\n';
        q.g=acceptedFull;
    }
    std::cout<<"{\"status\":\"complete\",\"flips\":"<<flips<<",\"anomalies\":"<<bad
             <<",\"max_ratio_relative\":"<<std::setprecision(17)<<maxRatio
             <<",\"max_green_relative\":"<<maxGreen<<",\"max_imaginary\":"<<maxImag<<"}\n";
    return bad?3:0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
