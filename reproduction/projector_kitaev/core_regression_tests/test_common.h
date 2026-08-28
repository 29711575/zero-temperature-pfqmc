#ifndef PFQMC_CORE_REGRESSION_TEST_COMMON_H
#define PFQMC_CORE_REGRESSION_TEST_COMMON_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"

class CoreTestWalker : public Spinless_tV {
public:
    int center=-1, ntrial=0, nphys=0;
    CoreTestWalker(const SpinlessTvChainUtils *c, rdGenerator *r,
                   double theta, double beta) {
        build_projector_static_contour(*this,c,r,theta,beta,center,ntrial,nphys);
    }
};

inline bool matrixFinite(const MatType &a) {
    for (Eigen::Index j=0;j<a.cols();++j)
        for (Eigen::Index i=0;i<a.rows();++i)
            if (!std::isfinite(a(i,j).real()) || !std::isfinite(a(i,j).imag())) return false;
    return true;
}

inline double relativeError(const MatType &a,const MatType &b) {
    return (a-b).norm()/std::max(b.norm(),std::numeric_limits<double>::min());
}

inline double structureResidual(const MatType &g) {
    MatType r=g+g.transpose()-2.0*MatType::Identity(g.rows(),g.cols());
    return r.norm()/std::max(g.norm(),std::numeric_limits<double>::min());
}

inline double maxAbs(const MatType &a) {
    double x=0;
    for (Eigen::Index j=0;j<a.cols();++j)
        for (Eigen::Index i=0;i<a.rows();++i) x=std::max(x,std::abs(a(i,j)));
    return x;
}

inline std::uint64_t hsHash(const std::vector<Operator*> &ops) {
    std::uint64_t h=1469598103934665603ULL;
    for (Operator *op:ops) if (iVecType *s=op->getAuxField())
        for (int i=0;i<s->size();++i) { h^=std::uint64_t((*s)(i)+2); h*=1099511628211ULL; }
    return h;
}

inline std::string boundaryKind(const CoreTestWalker &w,int boundary) {
    if (boundary==0) return "cyclic_boundary_0";
    if (boundary==w.center) return "center";
    if (boundary==w.ntrial) return "trial_physical_interface";
    if (boundary<w.ntrial) return "trial";
    return "physical";
}

inline std::string operatorKind(const CoreTestWalker &w,int boundary) {
    if (boundary<0 || boundary>=int(w.op_array.size())) return "cyclic";
    if (boundary<w.ntrial) return "trial_dense";
    if (dynamic_cast<SpinlessVOperator*>(w.op_array[boundary])) return "HS";
    return "physical_kinetic_half";
}

inline void dumpConfiguration(const std::vector<Operator*> &ops,
                              const std::string &path,
                              const std::string &header) {
    std::ofstream out(path.c_str());
    out<<header<<"\noperator,aux,sigma\n";
    for (std::size_t opi=0;opi<ops.size();++opi) if (iVecType *s=ops[opi]->getAuxField())
        for (int i=0;i<s->size();++i) out<<opi<<','<<i<<','<<(*s)(i)<<'\n';
}

#endif
