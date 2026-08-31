#ifndef PURE_PROJECTOR_MP_H
#define PURE_PROJECTOR_MP_H

#include "gaussian_trial_state.h"
#include "pure_projector_weight.h"

#include <boost/serialization/nvp.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

enum class PureMpProposalStatus {
    trusted,
    invalid_proposal,
    endpoint_untrusted,
    zero_ratio,
    nonfinite_ratio,
    complex_ratio,
    pfaffian_untrusted,
    precision_disagreement,
    unavailable
};

inline const char *pureMpProposalStatusName(PureMpProposalStatus status) {
    switch (status) {
    case PureMpProposalStatus::trusted: return "trusted";
    case PureMpProposalStatus::invalid_proposal: return "invalid_proposal";
    case PureMpProposalStatus::endpoint_untrusted: return "endpoint_untrusted";
    case PureMpProposalStatus::zero_ratio: return "zero_ratio";
    case PureMpProposalStatus::nonfinite_ratio: return "nonfinite_ratio";
    case PureMpProposalStatus::complex_ratio: return "complex_ratio";
    case PureMpProposalStatus::pfaffian_untrusted: return "pfaffian_untrusted";
    case PureMpProposalStatus::precision_disagreement: return "precision_disagreement";
    case PureMpProposalStatus::unavailable: return "unavailable";
    }
    return "unknown";
}

struct PureMpPerformanceProfile {
    std::uint64_t canonical_input_builds = 0;
    std::uint64_t cache_invalidations = 0;
    std::uint64_t operator_requests = 0;
    std::uint64_t operator_cache_hits = 0;
    std::uint64_t operator_cache_misses = 0;
    std::uint64_t sparse_apply_count = 0;
    std::uint64_t dense_apply_count = 0;
    double canonicalization_seconds = 0.0;
    double precision_160_seconds = 0.0;
    double precision_320_seconds = 0.0;
    double precision_640_seconds = 0.0;
    double conversion_seconds = 0.0;
    double propagation_seconds = 0.0;
    double thin_qr_seconds = 0.0;
    double endpoint_seconds = 0.0;
    double local_pfaffian_seconds = 0.0;

    PureMpPerformanceProfile &operator+=(const PureMpPerformanceProfile &other) {
        operator_requests+=other.operator_requests;
        canonical_input_builds+=other.canonical_input_builds;
        cache_invalidations+=other.cache_invalidations;
        operator_cache_hits+=other.operator_cache_hits;
        operator_cache_misses+=other.operator_cache_misses;
        sparse_apply_count+=other.sparse_apply_count;
        dense_apply_count+=other.dense_apply_count;
        canonicalization_seconds+=other.canonicalization_seconds;
        precision_160_seconds+=other.precision_160_seconds;
        precision_320_seconds+=other.precision_320_seconds;
        precision_640_seconds+=other.precision_640_seconds;
        conversion_seconds+=other.conversion_seconds;
        propagation_seconds+=other.propagation_seconds;
        thin_qr_seconds+=other.thin_qr_seconds;
        endpoint_seconds+=other.endpoint_seconds;
        local_pfaffian_seconds+=other.local_pfaffian_seconds;
        return *this;
    }
};

struct PureMpProposalResult {
    DataType ratio = DataType(0.0,0.0);
    int ratio_z2 = 0;
    double log_abs_ratio = -std::numeric_limits<double>::infinity();
    PureMpProposalStatus status = PureMpProposalStatus::unavailable;
    int precision_digits = 0;
    bool converged = false;
    double reality_error = std::numeric_limits<double>::infinity();
    double pre_endpoint_rcond = 0.0;
    double post_endpoint_rcond = 0.0;
    double pre_endpoint_residual = std::numeric_limits<double>::infinity();
    double post_endpoint_residual = std::numeric_limits<double>::infinity();
    MatType trusted_pre_green;
    MatType trusted_post_green;
    std::string canonical_contour_order =
        "ket_action_then_bra_strict_adjoint_reverse";
    std::string message;
    PureMpPerformanceProfile profile;
    bool ok() const { return status==PureMpProposalStatus::trusted && converged; }
};

struct PureMpProposalOptions {
    double minimum_endpoint_rcond = 1e-12;
    double residual_tolerance = 1e-10;
    double reality_tolerance = 1e-10;
    double zero_tolerance = 1e-14;
    double agreement_tolerance = 1e-12;
    bool real_z2 = true;
    bool enable_operator_cache = true;
};

namespace pure_projector_mp_detail {

struct CanonicalDelta {
    int row=0,column=0;
    DataType value=DataType(0.0,0.0);
};

struct CanonicalOperator {
    MatType matrix;
    std::vector<CanonicalDelta> delta_from_identity;
    bool sparse_delta=false;
};

struct CanonicalInput {
    std::vector<CanonicalOperator> operators;
    std::vector<int> slice_operator_ids;
    int candidate_operator_id=-1;
};

inline std::uint64_t canonicalMatrixHash(const MatType &matrix) {
    std::uint64_t hash=1469598103934665603ULL;
    auto mix=[&](std::uint64_t value){hash^=value;hash*=1099511628211ULL;};
    mix(std::uint64_t(matrix.rows()));mix(std::uint64_t(matrix.cols()));
    for(int c=0;c<matrix.cols();++c)for(int r=0;r<matrix.rows();++r){
        double re=matrix(r,c).real(),im=matrix(r,c).imag();
        if(re==0.0)re=0.0;if(im==0.0)im=0.0;
        std::uint64_t a=0,b=0;std::memcpy(&a,&re,sizeof(double));
        std::memcpy(&b,&im,sizeof(double));mix(a);mix(b);
    }
    return hash;
}

inline bool canonicalMatrixEqual(const MatType&a,const MatType&b) {
    if(a.rows()!=b.rows()||a.cols()!=b.cols())return false;
    for(int c=0;c<a.cols();++c)for(int r=0;r<a.rows();++r)
        if(a(r,c)!=b(r,c))return false;
    return true;
}

inline CanonicalOperator makeCanonicalOperator(const MatType &matrix) {
    CanonicalOperator result;result.matrix=matrix;
    if(matrix.rows()!=matrix.cols())return result;
    for(int r=0;r<matrix.rows();++r)for(int c=0;c<matrix.cols();++c){
        const DataType identity=r==c?DataType(1.0,0.0):DataType(0.0,0.0);
        const DataType delta=matrix(r,c)-identity;
        if(delta!=DataType(0.0,0.0))result.delta_from_identity.push_back({r,c,delta});
    }
    result.sparse_delta=result.delta_from_identity.size()<
        std::size_t(matrix.rows()*matrix.cols()/2);
    return result;
}

inline CanonicalInput canonicalizeInput(
    const std::vector<PureProjectorSlice>&slices,const MatType&candidate) {
    CanonicalInput result;std::unordered_map<std::uint64_t,std::vector<int>> buckets;
    auto intern=[&](const MatType &matrix){const std::uint64_t hash=canonicalMatrixHash(matrix);
        auto &ids=buckets[hash];for(int id:ids)if(canonicalMatrixEqual(
            result.operators[std::size_t(id)].matrix,matrix))return id;
        const int id=int(result.operators.size());
        result.operators.push_back(makeCanonicalOperator(matrix));ids.push_back(id);return id;};
    result.slice_operator_ids.reserve(slices.size());
    for(const auto &slice:slices)result.slice_operator_ids.push_back(intern(slice.matrix));
    result.candidate_operator_id=intern(candidate);return result;
}

template<unsigned Digits> struct Number {
    using Real=boost::multiprecision::number<
        boost::multiprecision::cpp_dec_float<Digits>>;
    struct Complex {
        Real r=0,i=0;
        Complex()=default;
        Complex(const Real &real):r(real),i(0){}
        Complex(const Real &real,const Real &imag):r(real),i(imag){}
    };
    struct Matrix {
        int rows=0,cols=0;
        std::vector<Complex> values;
        Matrix()=default;
        Matrix(int r,int c):rows(r),cols(c),values(std::size_t(r)*c){}
        Complex &operator()(int r,int c){return values[std::size_t(r)*cols+c];}
        const Complex &operator()(int r,int c)const{return values[std::size_t(r)*cols+c];}
        static Matrix identity(int n){Matrix x(n,n);for(int i=0;i<n;++i)x(i,i)=Complex(1);return x;}
    };
    struct SolveResult {bool ok=false;Matrix solution;Real residual=0,rcond=0;};

    static Complex add(Complex a,const Complex&b){a.r+=b.r;a.i+=b.i;return a;}
    static Complex sub(Complex a,const Complex&b){a.r-=b.r;a.i-=b.i;return a;}
    static Complex neg(Complex a){a.r=-a.r;a.i=-a.i;return a;}
    static Complex conjugate(Complex a){a.i=-a.i;return a;}
    static Complex multiply(const Complex&a,const Complex&b){
        return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};}
    static Real abs2(const Complex&a){return a.r*a.r+a.i*a.i;}
    static Real abs(const Complex&a){using boost::multiprecision::sqrt;return sqrt(abs2(a));}
    static Complex divide(const Complex&a,const Complex&b){const Real d=abs2(b);
        return {(a.r*b.r+a.i*b.i)/d,(a.i*b.r-a.r*b.i)/d};}
    static Complex scale(Complex a,const Real&s){a.r*=s;a.i*=s;return a;}
    static Complex unit(const Complex&a){return divide(a,Complex(abs(a)));}
    static Complex squareRoot(const Complex&z){
        using boost::multiprecision::sqrt;
        const Real magnitude=abs(z);
        Real positiveReal=(magnitude+z.r)/2;
        Real positiveImag=(magnitude-z.r)/2;
        if(positiveReal<0)positiveReal=0;if(positiveImag<0)positiveImag=0;
        Real realPart=sqrt(positiveReal);
        Real imagPart=sqrt(positiveImag);
        if(z.i<0)imagPart=-imagPart;
        return {realPart,imagPart};
    }
    static Matrix fromDouble(const MatType&a){Matrix x(a.rows(),a.cols());
        for(int r=0;r<a.rows();++r)for(int c=0;c<a.cols();++c)
            x(r,c)=Complex(Real(a(r,c).real()),Real(a(r,c).imag()));return x;}
    static MatType toDouble(const Matrix&a){MatType x(a.rows,a.cols);
        for(int r=0;r<a.rows;++r)for(int c=0;c<a.cols;++c)
            x(r,c)=DataType(a(r,c).r.template convert_to<double>(),
                            a(r,c).i.template convert_to<double>());return x;}
    static Matrix adjoint(const Matrix&a){Matrix x(a.cols,a.rows);
        for(int r=0;r<a.rows;++r)for(int c=0;c<a.cols;++c)x(c,r)=conjugate(a(r,c));return x;}
    static Matrix transpose(const Matrix&a){Matrix x(a.cols,a.rows);
        for(int r=0;r<a.rows;++r)for(int c=0;c<a.cols;++c)x(c,r)=a(r,c);return x;}
    static Matrix multiply(const Matrix&a,const Matrix&b){Matrix x(a.rows,b.cols);
        for(int r=0;r<a.rows;++r)for(int k=0;k<a.cols;++k){const Complex z=a(r,k);
            for(int c=0;c<b.cols;++c)x(r,c)=add(x(r,c),multiply(z,b(k,c)));}return x;}
    static Matrix add(const Matrix&a,const Matrix&b){Matrix x(a.rows,a.cols);
        for(std::size_t i=0;i<x.values.size();++i)x.values[i]=add(a.values[i],b.values[i]);return x;}
    static Matrix sub(const Matrix&a,const Matrix&b){Matrix x(a.rows,a.cols);
        for(std::size_t i=0;i<x.values.size();++i)x.values[i]=sub(a.values[i],b.values[i]);return x;}
    static Matrix scale(Matrix a,const Real&s){for(auto&z:a.values)z=scale(z,s);return a;}
    static Real norm(const Matrix&a){Real total=0;for(const auto&z:a.values)total+=abs2(z);
        using boost::multiprecision::sqrt;return sqrt(total);}

    // Rectangular modified Gram-Schmidt with reorthogonalization.  It changes
    // only the column gauge of the propagated pure subspace; the same gauge is
    // shared by pre/post local overlaps and cancels from their determinant
    // ratio.  No Green function is constructed at these intermediate points.
    static bool orthonormalize(Matrix &a){
        if(a.rows<a.cols||a.cols<=0)return false;
        const Real floor=boost::multiprecision::pow(Real(10),-int(Digits*2/3));
        for(int column=0;column<a.cols;++column){
            for(int pass=0;pass<2;++pass)for(int prior=0;prior<column;++prior){
                Complex projection;
                for(int row=0;row<a.rows;++row)projection=add(projection,
                    multiply(conjugate(a(row,prior)),a(row,column)));
                for(int row=0;row<a.rows;++row)a(row,column)=sub(a(row,column),
                    multiply(a(row,prior),projection));
            }
            Real length2=0;for(int row=0;row<a.rows;++row)length2+=abs2(a(row,column));
            using boost::multiprecision::sqrt;const Real length=sqrt(length2);
            if(!(length>floor))return false;
            for(int row=0;row<a.rows;++row)a(row,column)=scale(a(row,column),Real(1)/length);
        }
        return true;
    }

    static SolveResult solve(const Matrix&input,const Matrix&rhs){
        SolveResult result;if(input.rows!=input.cols||rhs.rows!=input.rows)return result;
        const int n=input.rows;Matrix a=input,b=rhs;
        const Real pivotFloor=boost::multiprecision::pow(Real(10),-int(Digits*3/4));
        for(int k=0;k<n;++k){int pivot=k;Real best=abs(a(k,k));
            for(int r=k+1;r<n;++r)if(abs(a(r,k))>best){best=abs(a(r,k));pivot=r;}
            if(!(best>pivotFloor))return result;
            if(pivot!=k)for(int c=0;c<n;++c)std::swap(a(k,c),a(pivot,c));
            if(pivot!=k)for(int c=0;c<b.cols;++c)std::swap(b(k,c),b(pivot,c));
            const Complex diagonal=a(k,k);
            for(int c=k;c<n;++c)a(k,c)=divide(a(k,c),diagonal);
            for(int c=0;c<b.cols;++c)b(k,c)=divide(b(k,c),diagonal);
            for(int r=0;r<n;++r)if(r!=k){const Complex factor=a(r,k);
                if(abs2(factor)==0)continue;
                for(int c=k;c<n;++c)a(r,c)=sub(a(r,c),multiply(factor,a(k,c)));
                for(int c=0;c<b.cols;++c)b(r,c)=sub(b(r,c),multiply(factor,b(k,c)));}}
        const Matrix residualMatrix=sub(multiply(input,b),rhs);
        result.residual=norm(residualMatrix)/std::max(Real(1),norm(rhs));
        if(rhs.rows==rhs.cols){const Real condition=norm(input)*norm(b);
            result.rcond=condition>0?Real(1)/condition:Real(0);}
        result.solution=std::move(b);result.ok=true;return result;
    }

    static Complex determinant(Matrix a,bool&ok){ok=false;if(a.rows!=a.cols)return {};
        Complex value(1);int parity=1;const int n=a.rows;
        const Real pivotFloor=boost::multiprecision::pow(Real(10),-int(Digits*3/4));
        for(int k=0;k<n;++k){int pivot=k;Real best=abs(a(k,k));
            for(int r=k+1;r<n;++r)if(abs(a(r,k))>best){best=abs(a(r,k));pivot=r;}
            if(!(best>pivotFloor))return {};
            if(pivot!=k){for(int c=k;c<n;++c)std::swap(a(k,c),a(pivot,c));parity=-parity;}
            const Complex diagonal=a(k,k);value=multiply(value,diagonal);
            for(int r=k+1;r<n;++r){const Complex factor=divide(a(r,k),diagonal);
                for(int c=k+1;c<n;++c)a(r,c)=sub(a(r,c),multiply(factor,a(k,c)));}}
        if(parity<0)value=neg(value);ok=true;return value;
    }

    static Complex pfaffian(Matrix a,bool&ok){ok=false;if(a.rows!=a.cols||a.rows%2)return {};
        Complex value(1);const int n=a.rows;
        const Real pivotFloor=boost::multiprecision::pow(Real(10),-int(Digits*3/4));
        for(int k=0;k<n;k+=2){int pivot=k+1;Real best=abs(a(k,pivot));
            for(int c=k+2;c<n;++c)if(abs(a(k,c))>best){best=abs(a(k,c));pivot=c;}
            if(!(best>pivotFloor))return {};
            if(pivot!=k+1){for(int c=0;c<n;++c)std::swap(a(k+1,c),a(pivot,c));
                for(int r=0;r<n;++r)std::swap(a(r,k+1),a(r,pivot));value=neg(value);}
            const Complex p=a(k,k+1);value=multiply(value,p);
            for(int i=k+2;i<n;++i)for(int j=i+1;j<n;++j){
                const Complex correction=divide(sub(multiply(a(k,i),a(k+1,j)),
                    multiply(a(k,j),a(k+1,i))),p);
                a(i,j)=sub(a(i,j),correction);a(j,i)=neg(a(i,j));}}
        ok=true;return value;
    }
};

template<unsigned Digits> PureMpProposalResult evaluateAtPrecisionLegacy(
    const GaussianTrialState &trial,const std::vector<PureProjectorSlice>&slices,
    int index,const MatType&newFactor,DataType newEta,
    const PureMpProposalOptions&options) {
    using X=Number<Digits>;using C=typename X::Complex;using M=typename X::Matrix;
    PureMpProposalResult result;result.precision_digits=Digits;
    try {
        if(index<0||index>=int(slices.size())||newFactor.rows()!=trial.Phi.rows()||
           newFactor.cols()!=trial.Phi.rows()) {result.status=PureMpProposalStatus::invalid_proposal;
            result.message="invalid frozen proposal";return result;}
        M right=X::fromDouble(trial.Phi);int rightSteps=0;
        for(int i=0;i<index;++i){right=X::multiply(X::fromDouble(slices[i].matrix),right);
            if(++rightSteps%8==0&&!X::orthonormalize(right)){
                result.status=PureMpProposalStatus::endpoint_untrusted;
                result.message="MP right thin-subspace QR failed";return result;}}
        if(!X::orthonormalize(right)){result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP right target subspace QR failed";return result;}
        M left=X::fromDouble(trial.Phi);
        int leftSteps=0;for(int i=int(slices.size())-1;i>index;--i){
            left=X::multiply(X::adjoint(X::fromDouble(slices[i].matrix)),left);
            if(++leftSteps%8==0&&!X::orthonormalize(left)){
                result.status=PureMpProposalStatus::endpoint_untrusted;
                result.message="MP left thin-subspace QR failed";return result;}}
        if(!X::orthonormalize(left)){result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP left target subspace QR failed";return result;}
        const M leftAdj=X::adjoint(left);
        const M oldFactor=X::fromDouble(slices[index].matrix),candidate=X::fromDouble(newFactor);
        const M oldRight=X::multiply(oldFactor,right),newRight=X::multiply(candidate,right);
        const M oldOverlap=X::multiply(leftAdj,oldRight),newOverlap=X::multiply(leftAdj,newRight);
        auto oldSolve=X::solve(oldOverlap,leftAdj),newSolve=X::solve(newOverlap,leftAdj);
        auto oldCond=X::solve(oldOverlap,M::identity(oldOverlap.rows));
        auto newCond=X::solve(newOverlap,M::identity(newOverlap.rows));
        if(!oldSolve.ok||!newSolve.ok||!oldCond.ok||!newCond.ok){
            result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP endpoint overlap solve failed";return result;}
        result.pre_endpoint_rcond=oldCond.rcond.template convert_to<double>();
        result.post_endpoint_rcond=newCond.rcond.template convert_to<double>();
        result.pre_endpoint_residual=oldSolve.residual.template convert_to<double>();
        result.post_endpoint_residual=newSolve.residual.template convert_to<double>();
        if(result.pre_endpoint_rcond<options.minimum_endpoint_rcond||
           result.post_endpoint_rcond<options.minimum_endpoint_rcond||
           result.pre_endpoint_residual>options.residual_tolerance||
           result.post_endpoint_residual>options.residual_tolerance){
            result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP endpoint solve failed trust thresholds";return result;}
        bool oldDetOk=false,newDetOk=false;
        const C oldDet=X::determinant(oldOverlap,oldDetOk);
        const C newDet=X::determinant(newOverlap,newDetOk);
        if(!oldDetOk||!newDetOk){result.status=PureMpProposalStatus::zero_ratio;
            result.message="MP endpoint determinant is zero";return result;}
        const C etaRatio(C(typename X::Real(newEta.real()),typename X::Real(newEta.imag())));
        const C oldEta(typename X::Real(slices[index].eta.real()),
                       typename X::Real(slices[index].eta.imag()));
        const C normalizedEta=X::divide(etaRatio,oldEta);
        C ratioSquared=X::multiply(X::multiply(normalizedEta,normalizedEta),X::divide(newDet,oldDet));
        C ratio=X::squareRoot(ratioSquared);

        auto oldInverse=X::solve(oldFactor,M::identity(oldFactor.rows));
        if(!oldInverse.ok){result.status=PureMpProposalStatus::pfaffian_untrusted;
            result.message="MP local old factor solve failed";return result;}
        const M transform=X::multiply(candidate,oldInverse.solution);
        auto cayley=X::solve(X::add(M::identity(transform.rows),transform),
                             X::scale(M::identity(transform.rows),typename X::Real(2)));
        if(!cayley.ok){result.status=PureMpProposalStatus::pfaffian_untrusted;
            result.message="MP local Cayley solve failed";return result;}
        const M transformGreen=X::sub(cayley.solution,M::identity(transform.rows));
        const M oldGreen=X::sub(M::identity(oldRight.rows),
            X::scale(X::multiply(oldRight,oldSolve.solution),typename X::Real(2)));
        M block(2*oldGreen.rows,2*oldGreen.rows);
        for(int r=0;r<oldGreen.rows;++r)for(int c=0;c<oldGreen.cols;++c){
            block(r,c)=transformGreen(r,c);block(oldGreen.rows+r,oldGreen.cols+c)=oldGreen(r,c);}
        for(int i=0;i<oldGreen.rows;++i){block(i,oldGreen.rows+i)=C(-1);
            block(oldGreen.rows+i,i)=C(1);}
        bool pfOk=false;const C pf=X::pfaffian(block,pfOk);
        if(!pfOk){result.status=PureMpProposalStatus::pfaffian_untrusted;
            result.message="MP low-rank Pfaffian pivot failed";return result;}
        C target=X::unit(X::multiply(normalizedEta,pf));
        if((trial.Phi.cols()%2)!=0)target=X::neg(target);
        const C aligned=X::multiply(X::conjugate(target),ratio);
        if(aligned.r<0)ratio=X::neg(ratio);
        const typename X::Real ratioMagnitude=X::abs(ratio);
        if(!(ratioMagnitude>typename X::Real(options.zero_tolerance))){
            result.status=PureMpProposalStatus::zero_ratio;result.message="MP ratio is zero";return result;}
        result.ratio=DataType(ratio.r.template convert_to<double>(),ratio.i.template convert_to<double>());
        result.ratio_z2=ratio.r<0?-1:1;
        using boost::multiprecision::log;
        result.log_abs_ratio=log(ratioMagnitude).template convert_to<double>();
        result.reality_error=(X::abs(C(ratio.i))/std::max(X::abs(C(ratio.r)),
            typename X::Real("1e-1000"))).template convert_to<double>();
        if(options.real_z2&&result.reality_error>options.reality_tolerance){
            result.status=PureMpProposalStatus::complex_ratio;
            result.message="MP ratio is not resolved on the real axis";return result;}
        const M newGreen=X::sub(M::identity(newRight.rows),
            X::scale(X::multiply(newRight,newSolve.solution),typename X::Real(2)));
        result.trusted_pre_green=X::toDouble(oldGreen);
        result.trusted_post_green=X::toDouble(newGreen);
        result.status=PureMpProposalStatus::precision_disagreement;
        result.message="single precision requires consecutive-precision agreement";
        return result;
    } catch(const std::exception&e) {
        result.status=PureMpProposalStatus::unavailable;result.message=e.what();return result;
    }
}

template<unsigned Digits> PureMpProposalResult evaluateAtPrecisionCached(
    const GaussianTrialState &trial,const std::vector<PureProjectorSlice>&slices,
    int index,const MatType&newFactor,DataType newEta,
    const PureMpProposalOptions&options,const CanonicalInput&canonical) {
    using X=Number<Digits>;using C=typename X::Complex;using M=typename X::Matrix;
    struct Entry {int row=0,column=0;C value;};
    struct Operator {M dense;std::vector<Entry> delta;bool sparse=false;};
    PureMpProposalResult result;result.precision_digits=Digits;
    try {
        if(index<0||index>=int(slices.size())||newFactor.rows()!=trial.Phi.rows()||
           newFactor.cols()!=trial.Phi.rows()||canonical.slice_operator_ids.size()!=slices.size()||
           canonical.candidate_operator_id<0) {
            result.status=PureMpProposalStatus::invalid_proposal;
            result.message="invalid cached frozen proposal";return result;}

        const auto conversionStarted=std::chrono::steady_clock::now();
        std::vector<Operator> operators;operators.reserve(canonical.operators.size());
        for(const auto &source:canonical.operators){Operator converted;
            converted.dense=X::fromDouble(source.matrix);converted.sparse=source.sparse_delta;
            converted.delta.reserve(source.delta_from_identity.size());
            for(const auto &entry:source.delta_from_identity)
                converted.delta.push_back({entry.row,entry.column,
                    C(typename X::Real(entry.value.real()),typename X::Real(entry.value.imag()))});
            operators.push_back(std::move(converted));}
        result.profile.conversion_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-conversionStarted).count();
        result.profile.operator_requests=slices.size()+1;
        result.profile.operator_cache_misses=operators.size();
        result.profile.operator_cache_hits=result.profile.operator_requests>
            result.profile.operator_cache_misses?result.profile.operator_requests-
            result.profile.operator_cache_misses:0;

        auto apply=[&](const Operator&op,const M&input,bool useAdjoint){
            const auto started=std::chrono::steady_clock::now();M output;
            if(op.sparse){++result.profile.sparse_apply_count;output=input;
                if(!useAdjoint){for(const auto&e:op.delta)for(int c=0;c<input.cols;++c)
                    output(e.row,c)=X::add(output(e.row,c),X::multiply(e.value,input(e.column,c)));}
                else{for(const auto&e:op.delta)for(int c=0;c<input.cols;++c)
                    output(e.column,c)=X::add(output(e.column,c),
                        X::multiply(X::conjugate(e.value),input(e.row,c)));}}
            else{++result.profile.dense_apply_count;output=X::multiply(
                useAdjoint?X::adjoint(op.dense):op.dense,input);}
            result.profile.propagation_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count();return output;};
        auto orthonormalize=[&](M&subspace){const auto started=std::chrono::steady_clock::now();
            const bool ok=X::orthonormalize(subspace);
            result.profile.thin_qr_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count();return ok;};

        M right=X::fromDouble(trial.Phi);int rightSteps=0;
        for(int i=0;i<index;++i){right=apply(operators[std::size_t(
                canonical.slice_operator_ids[std::size_t(i)])],right,false);
            if(++rightSteps%8==0&&!orthonormalize(right)){
                result.status=PureMpProposalStatus::endpoint_untrusted;
                result.message="MP right thin-subspace QR failed";return result;}}
        if(!orthonormalize(right)){result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP right target subspace QR failed";return result;}
        M left=X::fromDouble(trial.Phi);int leftSteps=0;
        for(int i=int(slices.size())-1;i>index;--i){left=apply(operators[std::size_t(
                canonical.slice_operator_ids[std::size_t(i)])],left,true);
            if(++leftSteps%8==0&&!orthonormalize(left)){
                result.status=PureMpProposalStatus::endpoint_untrusted;
                result.message="MP left thin-subspace QR failed";return result;}}
        if(!orthonormalize(left)){result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP left target subspace QR failed";return result;}

        const auto endpointStarted=std::chrono::steady_clock::now();
        const M leftAdj=X::adjoint(left);
        const Operator&oldOperator=operators[std::size_t(
            canonical.slice_operator_ids[std::size_t(index)])];
        const Operator&candidateOperator=operators[std::size_t(canonical.candidate_operator_id)];
        const M oldRight=apply(oldOperator,right,false),newRight=apply(candidateOperator,right,false);
        const M oldOverlap=X::multiply(leftAdj,oldRight),newOverlap=X::multiply(leftAdj,newRight);
        auto oldSolve=X::solve(oldOverlap,leftAdj),newSolve=X::solve(newOverlap,leftAdj);
        auto oldCond=X::solve(oldOverlap,M::identity(oldOverlap.rows));
        auto newCond=X::solve(newOverlap,M::identity(newOverlap.rows));
        if(!oldSolve.ok||!newSolve.ok||!oldCond.ok||!newCond.ok){
            result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP endpoint overlap solve failed";return result;}
        result.pre_endpoint_rcond=oldCond.rcond.template convert_to<double>();
        result.post_endpoint_rcond=newCond.rcond.template convert_to<double>();
        result.pre_endpoint_residual=oldSolve.residual.template convert_to<double>();
        result.post_endpoint_residual=newSolve.residual.template convert_to<double>();
        if(result.pre_endpoint_rcond<options.minimum_endpoint_rcond||
           result.post_endpoint_rcond<options.minimum_endpoint_rcond||
           result.pre_endpoint_residual>options.residual_tolerance||
           result.post_endpoint_residual>options.residual_tolerance){
            result.status=PureMpProposalStatus::endpoint_untrusted;
            result.message="MP endpoint solve failed trust thresholds";return result;}
        bool oldDetOk=false,newDetOk=false;
        const C oldDet=X::determinant(oldOverlap,oldDetOk);
        const C newDet=X::determinant(newOverlap,newDetOk);
        if(!oldDetOk||!newDetOk){result.status=PureMpProposalStatus::zero_ratio;
            result.message="MP endpoint determinant is zero";return result;}
        const C etaRatio(C(typename X::Real(newEta.real()),typename X::Real(newEta.imag())));
        const C oldEta(typename X::Real(slices[index].eta.real()),
                       typename X::Real(slices[index].eta.imag()));
        const C normalizedEta=X::divide(etaRatio,oldEta);
        C ratioSquared=X::multiply(X::multiply(normalizedEta,normalizedEta),X::divide(newDet,oldDet));
        C ratio=X::squareRoot(ratioSquared);
        result.profile.endpoint_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-endpointStarted).count();

        const auto localStarted=std::chrono::steady_clock::now();
        auto oldInverse=X::solve(oldOperator.dense,M::identity(oldOperator.dense.rows));
        if(!oldInverse.ok){result.status=PureMpProposalStatus::pfaffian_untrusted;
            result.message="MP local old factor solve failed";return result;}
        const M transform=X::multiply(candidateOperator.dense,oldInverse.solution);
        auto cayley=X::solve(X::add(M::identity(transform.rows),transform),
                             X::scale(M::identity(transform.rows),typename X::Real(2)));
        if(!cayley.ok){result.status=PureMpProposalStatus::pfaffian_untrusted;
            result.message="MP local Cayley solve failed";return result;}
        const M transformGreen=X::sub(cayley.solution,M::identity(transform.rows));
        const M oldGreen=X::sub(M::identity(oldRight.rows),
            X::scale(X::multiply(oldRight,oldSolve.solution),typename X::Real(2)));
        M block(2*oldGreen.rows,2*oldGreen.rows);
        for(int r=0;r<oldGreen.rows;++r)for(int c=0;c<oldGreen.cols;++c){
            block(r,c)=transformGreen(r,c);block(oldGreen.rows+r,oldGreen.cols+c)=oldGreen(r,c);}
        for(int i=0;i<oldGreen.rows;++i){block(i,oldGreen.rows+i)=C(-1);
            block(oldGreen.rows+i,i)=C(1);}
        bool pfOk=false;const C pf=X::pfaffian(block,pfOk);
        if(!pfOk){result.status=PureMpProposalStatus::pfaffian_untrusted;
            result.message="MP low-rank Pfaffian pivot failed";return result;}
        C target=X::unit(X::multiply(normalizedEta,pf));
        if((trial.Phi.cols()%2)!=0)target=X::neg(target);
        const C aligned=X::multiply(X::conjugate(target),ratio);
        if(aligned.r<0)ratio=X::neg(ratio);
        const typename X::Real ratioMagnitude=X::abs(ratio);
        if(!(ratioMagnitude>typename X::Real(options.zero_tolerance))){
            result.status=PureMpProposalStatus::zero_ratio;result.message="MP ratio is zero";return result;}
        result.ratio=DataType(ratio.r.template convert_to<double>(),ratio.i.template convert_to<double>());
        result.ratio_z2=ratio.r<0?-1:1;
        using boost::multiprecision::log;
        result.log_abs_ratio=log(ratioMagnitude).template convert_to<double>();
        result.reality_error=(X::abs(C(ratio.i))/std::max(X::abs(C(ratio.r)),
            typename X::Real("1e-1000"))).template convert_to<double>();
        if(options.real_z2&&result.reality_error>options.reality_tolerance){
            result.status=PureMpProposalStatus::complex_ratio;
            result.message="MP ratio is not resolved on the real axis";return result;}
        const M newGreen=X::sub(M::identity(newRight.rows),
            X::scale(X::multiply(newRight,newSolve.solution),typename X::Real(2)));
        result.trusted_pre_green=X::toDouble(oldGreen);
        result.trusted_post_green=X::toDouble(newGreen);
        result.profile.local_pfaffian_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-localStarted).count();
        result.status=PureMpProposalStatus::precision_disagreement;
        result.message="single precision requires consecutive-precision agreement";
        return result;
    } catch(const std::exception&e) {
        result.status=PureMpProposalStatus::unavailable;result.message=e.what();return result;
    }
}

inline bool resolved(const PureMpProposalResult&r) {
    return r.status==PureMpProposalStatus::precision_disagreement&&
        std::isfinite(r.ratio.real())&&std::isfinite(r.ratio.imag())&&r.ratio_z2!=0;
}

inline bool agree(const PureMpProposalResult&a,const PureMpProposalResult&b,
                  const PureMpProposalOptions&o) {
    if(!resolved(a)||!resolved(b)||a.ratio_z2!=b.ratio_z2)return false;
    const double relative=std::abs(a.ratio-b.ratio)/std::max(o.zero_tolerance,std::abs(b.ratio));
    return std::isfinite(relative)&&relative<=o.agreement_tolerance&&
        (!o.real_z2||(a.reality_error<=o.reality_tolerance&&b.reality_error<=o.reality_tolerance));
}

} // namespace pure_projector_mp_detail

inline PureMpProposalResult pureProjectorMpSameProposal(
    const GaussianTrialState &trial,const std::vector<PureProjectorSlice>&slices,
    int index,const MatType&newFactor,DataType newEta,
    const PureMpProposalOptions&options=PureMpProposalOptions()) {
    using namespace pure_projector_mp_detail;
    // Phase 3E only optimizes the production RealZ2 oracle.  Keep the
    // GenericComplex implementation bit-for-bit on the legacy path.
    const bool useCache=options.enable_operator_cache&&options.real_z2;
    CanonicalInput canonical;double canonicalizationSeconds=0.0;
    if(useCache){const auto started=std::chrono::steady_clock::now();
        canonical=canonicalizeInput(slices,newFactor);
        canonicalizationSeconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();}
    auto run160=[&](){const auto started=std::chrono::steady_clock::now();
        PureMpProposalResult value=useCache?
            evaluateAtPrecisionCached<160>(trial,slices,index,newFactor,newEta,options,canonical):
            evaluateAtPrecisionLegacy<160>(trial,slices,index,newFactor,newEta,options);
        value.profile.precision_160_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();return value;};
    auto run320=[&](){const auto started=std::chrono::steady_clock::now();
        PureMpProposalResult value=useCache?
            evaluateAtPrecisionCached<320>(trial,slices,index,newFactor,newEta,options,canonical):
            evaluateAtPrecisionLegacy<320>(trial,slices,index,newFactor,newEta,options);
        value.profile.precision_320_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();return value;};
    PureMpProposalResult r160=run160(),r320=run320();
    PureMpPerformanceProfile combined=r160.profile;combined+=r320.profile;
    if(useCache){combined.canonical_input_builds=1;
        combined.cache_invalidations=1;}
    combined.canonicalization_seconds=canonicalizationSeconds;
    if(agree(r160,r320,options)){PureMpProposalResult result=r320;result.profile=combined;
        result.status=PureMpProposalStatus::trusted;result.converged=true;
        result.message="MP160/MP320 consecutive-precision agreement";return result;}
    const auto started640=std::chrono::steady_clock::now();
    PureMpProposalResult r640=useCache?
        evaluateAtPrecisionCached<640>(trial,slices,index,newFactor,newEta,options,canonical):
        evaluateAtPrecisionLegacy<640>(trial,slices,index,newFactor,newEta,options);
    r640.profile.precision_640_seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started640).count();combined+=r640.profile;
    combined.canonicalization_seconds=canonicalizationSeconds;
    if(agree(r320,r640,options)){PureMpProposalResult result=r640;result.profile=combined;
        result.status=PureMpProposalStatus::trusted;result.converged=true;
        result.message="MP320/MP640 consecutive-precision agreement";return result;}
    PureMpProposalResult result=r640;result.profile=combined;
    result.status=PureMpProposalStatus::precision_disagreement;
    result.converged=false;result.message="MP same-proposal oracle failed consecutive-precision agreement";
    return result;
}

#endif
