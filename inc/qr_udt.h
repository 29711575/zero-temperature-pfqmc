#ifndef QR_UDT_H
#define QR_UDT_H

#include "types.h"

#ifdef PFQMC_SCALE_SAFE_UDT

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Enabled by default for scale-safe UDT builds. A no-guard build is retained
// only for matched validation against the already validated production path.
#ifndef PFQMC_UDT_RANK_LOSS_GUARD_BITS
#define PFQMC_UDT_RANK_LOSS_GUARD_BITS 45
#endif

static constexpr double scaleSafeUDTRankLossGuardBits =
    double(PFQMC_UDT_RANK_LOSS_GUARD_BITS);
// The unchanged production MGS path reaches 1.8627e-7 on the audited
// n=12, +/-20 control while retaining 1e-16 reconstruction/solve residuals.
// The first known unsafe case is O(1).  Keep a conservative separating gate
// without changing the factorization through reorthogonalization.
static constexpr double scaleSafeUDTPartialQOrthogonalityLimit = 1.0e-6;

#ifndef PFQMC_UDT_ORTHOGONALITY_PRECHECK_BITS
#define PFQMC_UDT_ORTHOGONALITY_PRECHECK_BITS 32
#endif
// Real QMC reached 9.081953 bits; the first audited orthogonality alarm was
// 37.549 bits.  32 bits leaves >22 bits above real QMC and >5 bits before that
// synthetic alarm while retaining the formal 45-bit fail threshold unchanged.
static constexpr double scaleSafeUDTOrthogonalityPrecheckBits =
    double(PFQMC_UDT_ORTHOGONALITY_PRECHECK_BITS);

class ScaleSafeQRGuardFailure : public std::runtime_error
{
public:
    explicit ScaleSafeQRGuardFailure(const std::string &message)
        : std::runtime_error(message) {}
};

struct ScaleSafeQRGuardDiagnostics
{
    std::uint64_t trigger_count = 0;
    double max_lost_bits = 0.0;
    double min_guard_margin = std::numeric_limits<double>::infinity();
    int last_matrix_size = 0;
    int last_pivot = -1;
    double last_lost_bits = 0.0;
    double last_exponent_span = 0.0;
    double last_orthogonality_residual = 0.0;
};

inline ScaleSafeQRGuardDiagnostics &scaleSafeQRGuardDiagnostics()
{
    static ScaleSafeQRGuardDiagnostics diagnostics;
    return diagnostics;
}

inline void resetScaleSafeQRGuardDiagnostics()
{
    scaleSafeQRGuardDiagnostics() = ScaleSafeQRGuardDiagnostics();
}

[[noreturn]] inline void scaleSafeQRGuardFail(const char *reason, int n, int pivot,
    double lost_bits, double exponent_span, double orthogonality_residual)
{
    ScaleSafeQRGuardDiagnostics &diagnostics = scaleSafeQRGuardDiagnostics();
    ++diagnostics.trigger_count;
    diagnostics.last_matrix_size = n;
    diagnostics.last_pivot = pivot;
    diagnostics.last_lost_bits = lost_bits;
    diagnostics.last_exponent_span = exponent_span;
    diagnostics.last_orthogonality_residual = orthogonality_residual;
    std::ostringstream message;
    message << std::setprecision(17) << "scale-safe exponent QR guard: " << reason
            << "; n=" << n << "; pivot=" << pivot << "; lost_bits=" << lost_bits
            << "; exponent_span=" << exponent_span
            << "; orthogonality_residual=" << orthogonality_residual;
    throw ScaleSafeQRGuardFailure(message.str());
}

inline bool scaleSafeFinite(const MatType &a)
{
    for (Eigen::Index j = 0; j < a.cols(); ++j)
        for (Eigen::Index i = 0; i < a.rows(); ++i)
            if (!std::isfinite(a(i,j).real()) || !std::isfinite(a(i,j).imag())) return false;
    return true;
}

inline double scaleSafePow2(int exponent)
{
    return std::ldexp(1.0, exponent);
}

inline DataType scaleSafeLdexp(const DataType &z, int exponent)
{
    return DataType(std::ldexp(z.real(), exponent), std::ldexp(z.imag(), exponent));
}

// Solve A X = B after exact power-of-two row/column equilibration.  The
// equilibration changes only the numerical representation, not the equation.
inline MatType scaleSafeCheckedSolve(const MatType &A, const MatType &B,
                                     const char *site)
{
    if (A.rows() != A.cols() || A.rows() != B.rows() || !scaleSafeFinite(A) || !scaleSafeFinite(B))
        throw std::runtime_error(std::string("scale-safe solve invalid/nonfinite input at ") + site);
    const int n = static_cast<int>(A.rows());
    MatType As = A, Bs = B;
    std::vector<int> row_exp(n, 0), col_exp(n, 0);
    for (int i = 0; i < n; ++i) {
        double m = 0.0;
        for (int j = 0; j < n; ++j) m = std::max(m, std::abs(As(i,j)));
        if (!(m > 0.0) || !std::isfinite(m))
            throw std::runtime_error(std::string("scale-safe solve zero/nonfinite row at ") + site);
        int e = 0; std::frexp(m, &e); row_exp[i] = -e;
        const double s = scaleSafePow2(row_exp[i]); As.row(i) *= s; Bs.row(i) *= s;
    }
    for (int j = 0; j < n; ++j) {
        double m = 0.0;
        for (int i = 0; i < n; ++i) m = std::max(m, std::abs(As(i,j)));
        if (!(m > 0.0) || !std::isfinite(m))
            throw std::runtime_error(std::string("scale-safe solve zero/nonfinite column at ") + site);
        int e = 0; std::frexp(m, &e); col_exp[j] = -e;
        As.col(j) *= scaleSafePow2(col_exp[j]);
    }
    MatType lu = As;
    std::vector<lapack_int> piv(n);
    lapack_int info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, lu.data(), n, piv.data());
    if (info != 0)
        throw std::runtime_error(std::string("scale-safe zgetrf info=") + std::to_string(info) + " at " + site);
    const double anorm = LAPACKE_zlange(LAPACK_COL_MAJOR, '1', n, n, As.data(), n);
    double rcond = 0.0;
    info = LAPACKE_zgecon(LAPACK_COL_MAJOR, '1', n, lu.data(), n, anorm, &rcond);
    if (info != 0 || !std::isfinite(rcond) || rcond <= 0.0)
        throw std::runtime_error(std::string("scale-safe zgecon invalid at ") + site);
    MatType X = Bs;
    info = LAPACKE_zgetrs(LAPACK_COL_MAJOR, 'N', n, static_cast<int>(X.cols()),
                          lu.data(), n, piv.data(), X.data(), n);
    if (info != 0)
        throw std::runtime_error(std::string("scale-safe zgetrs info=") + std::to_string(info) + " at " + site);
    for (int i = 0; i < n; ++i) X.row(i) *= scaleSafePow2(col_exp[i]);
    if (!scaleSafeFinite(X))
        throw std::runtime_error(std::string("scale-safe solve produced nonfinite output at ") + site);
    return X;
}

inline MatType scaleSafeCheckedRightSolve(const MatType &B, const MatType &A,
                                          const char *site)
{
    return scaleSafeCheckedSolve(A.adjoint(), B.adjoint(), site).adjoint();
}

class UDT
{
public:
    int nDim = 0;
    MatType U;
    // D stores a positive mantissa in [0.5,1); Dexp stores its base-2 exponent.
    dVecType D;
    iVecType Dexp;
    MatType T;

    UDT() = default;
    explicit UDT(int n) : nDim(n), U(MatType::Identity(n,n)), D(dVecType::Constant(n,0.5)),
                          Dexp(iVecType::Ones(n)), T(MatType::Identity(n,n))
    {}

    UDT(const MatType &_U, const dVecType &_D, const MatType &_T)
        : nDim(static_cast<int>(_U.rows())), U(_U), D(nDim), Dexp(nDim), T(_T)
    {
        for (int i=0;i<nDim;++i) setActualD(i,_D(i));
    }

    explicit UDT(MatType &A) { factorDense(A); }

    void setActualD(int i, double value)
    {
        if (!(value > 0.0) || !std::isfinite(value))
            throw std::runtime_error("scale-safe UDT invalid QR diagonal");
        int e=0; D(i)=std::frexp(value,&e); Dexp(i)=e;
    }

    bool dGreaterThanOne(int i) const
    {
        return Dexp(i)>1 || (Dexp(i)==1 && D(i)>0.5);
    }

    double dSmallPart(int i) const
    {
        return dGreaterThanOne(i) ? 1.0 : std::ldexp(D(i),Dexp(i));
    }

    double dLargeInverse(int i) const
    {
        return dGreaterThanOne(i) ? std::ldexp(1.0/D(i),-Dexp(i)) : 1.0;
    }

    bool materializable(int limit=500) const
    {
        for (int i=0;i<nDim;++i) if (Dexp(i)<-limit || Dexp(i)>limit) return false;
        return true;
    }

    double actualD(int i) const { return std::ldexp(D(i),Dexp(i)); }

    void factorDense(MatType &A)
    {
        if (A.rows()!=A.cols() || !scaleSafeFinite(A))
            throw std::runtime_error("scale-safe dense QR invalid input");
        nDim=static_cast<int>(A.rows()); T=MatType::Zero(nDim,nDim); D.resize(nDim); Dexp.resize(nDim);
        std::vector<lapack_int> jpvt(nDim,0); std::vector<DataType> tau(nDim);
        lapack_int info=LAPACKE_zgeqp3(LAPACK_COL_MAJOR,nDim,nDim,A.data(),nDim,jpvt.data(),tau.data());
        if(info!=0) throw std::runtime_error("scale-safe zgeqp3 info="+std::to_string(info));
        for(int i=0;i<nDim;++i) {
            const double diag=std::abs(A(i,i).real());
            setActualD(i,diag);
            const double inv=1.0/diag;
            for(int j=i;j<nDim;++j) A(i,j)*=inv;
        }
        for(int i=0;i<nDim;++i) T(i,jpvt[i]-1)=1.0;
        cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
                    nDim,nDim,&one,A.data(),nDim,T.data(),nDim);
        info=LAPACKE_zungqr(LAPACK_COL_MAJOR,nDim,nDim,nDim,A.data(),nDim,tau.data());
        if(info!=0) throw std::runtime_error("scale-safe zungqr info="+std::to_string(info));
        U=A;
    }

    // CPQR for columns represented as work.col(j)*2^exponent(j).
    static UDT factorPreparedScaled(MatType work, iVecType exponent)
    {
        const int n=static_cast<int>(work.rows());
        if(work.cols()!=n || exponent.size()!=n || !scaleSafeFinite(work))
            throw std::runtime_error("scale-safe exponent QR invalid prepared columns");
        MatType q=MatType::Zero(n,n), t=MatType::Zero(n,n);
        std::vector<int> permutation(n);
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
        ScaleSafeQRGuardDiagnostics &guard=scaleSafeQRGuardDiagnostics();
        std::vector<double> original_log2_norm(n,0.0);
        double exponent_span=std::numeric_limits<double>::infinity();
        int active_pivot=-1;
        double active_lost_bits=std::numeric_limits<double>::infinity();
        double active_orthogonality=std::numeric_limits<double>::infinity();
#endif
        for(int j=0;j<n;++j) permutation[j]=j;
        const auto renormalize=[&](int j,double known_norm) {
            const double norm=known_norm>0.0?known_norm:work.col(j).stableNorm();
            if(!(norm>0.0) || !std::isfinite(norm)) {
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
                scaleSafeQRGuardFail("zero/nonfinite column during power-of-two normalization",
                    n,active_pivot>=0?active_pivot:j,active_lost_bits,
                    exponent_span,active_orthogonality);
#else
                throw std::runtime_error("scale-safe exponent QR rank loss/nonfinite column");
#endif
            }
            int e=0; std::frexp(norm,&e);
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
            const long long updated_exponent=static_cast<long long>(exponent(j))+e;
            if(updated_exponent<std::numeric_limits<int>::min() ||
               updated_exponent>std::numeric_limits<int>::max())
                scaleSafeQRGuardFail("column exponent bookkeeping overflow",n,
                    active_pivot>=0?active_pivot:j,active_lost_bits,
                    exponent_span,active_orthogonality);
#endif
            work.col(j)*=scaleSafePow2(-e); exponent(j)+=e;
        };
        for(int j=0;j<n;++j) {
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
            const double original_norm=work.col(j).stableNorm();
            if(!(original_norm>0.0) || !std::isfinite(original_norm))
                scaleSafeQRGuardFail("zero/nonfinite original column",n,j,
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity());
            original_log2_norm[j]=double(exponent(j))+std::log2(original_norm);
            renormalize(j,original_norm);
#else
            renormalize(j,-1.0);
#endif
        }
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
        exponent_span=*std::max_element(original_log2_norm.begin(),original_log2_norm.end())-
                      *std::min_element(original_log2_norm.begin(),original_log2_norm.end());
#endif
        UDT out; out.nDim=n; out.D.resize(n); out.Dexp.resize(n);
        for(int k=0;k<n;++k) {
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
            active_pivot=k;
            active_lost_bits=std::numeric_limits<double>::infinity();
            active_orthogonality=0.0;
#endif
            int pivot=k;
            for(int j=k+1;j<n;++j) {
                if(exponent(j)>exponent(pivot) ||
                   (exponent(j)==exponent(pivot) && work.col(j).squaredNorm()>work.col(pivot).squaredNorm())) pivot=j;
            }
            if(pivot!=k) {
                work.col(k).swap(work.col(pivot)); std::swap(exponent(k),exponent(pivot));
                std::swap(permutation[k],permutation[pivot]);
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
                std::swap(original_log2_norm[k],original_log2_norm[pivot]);
#endif
            }
            const double norm=work.col(k).stableNorm();
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
            if(!(norm>0.0) || !std::isfinite(norm))
                scaleSafeQRGuardFail("zero/nonfinite pivot residual",n,k,
                    std::numeric_limits<double>::infinity(),exponent_span,
                    std::numeric_limits<double>::infinity());
            const double residual_log2_norm=double(exponent(k))+std::log2(norm);
            const double lost_bits=original_log2_norm[k]-residual_log2_norm;
            active_lost_bits=lost_bits;
            guard.max_lost_bits=std::max(guard.max_lost_bits,lost_bits);
            guard.min_guard_margin=std::min(guard.min_guard_margin,
                scaleSafeUDTRankLossGuardBits-lost_bits);
            if(!std::isfinite(lost_bits) || lost_bits>scaleSafeUDTRankLossGuardBits)
                scaleSafeQRGuardFail("rank loss before normalization",n,k,lost_bits,
                    exponent_span,0.0);
#endif
            int diag_e=0; const double diag_m=std::frexp(norm,&diag_e);
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
            const long long combined_exponent=static_cast<long long>(exponent(k))+diag_e;
            if(!(diag_m>=0.5 && diag_m<1.0) || !std::isfinite(diag_m) ||
               combined_exponent<std::numeric_limits<int>::min() ||
               combined_exponent>std::numeric_limits<int>::max())
                scaleSafeQRGuardFail("invalid D mantissa/exponent bookkeeping",n,k,
                    lost_bits,exponent_span,0.0);
#endif
            out.D(k)=diag_m; out.Dexp(k)=exponent(k)+diag_e;
            q.col(k)=work.col(k)/norm; t(k,permutation[k])=1.0;
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
            double orthogonality_residual=0.0;
            for(int i=0;i<n;++i)
                if(!std::isfinite(q(i,k).real()) || !std::isfinite(q(i,k).imag()))
                    scaleSafeQRGuardFail("nonfinite Q column after normalization",n,k,
                        lost_bits,exponent_span,std::numeric_limits<double>::infinity());
            if(lost_bits>=scaleSafeUDTOrthogonalityPrecheckBits) {
                const MatType partial_q=q.leftCols(k+1);
                const MatType partial_gram=partial_q.adjoint()*partial_q-
                    MatType::Identity(k+1,k+1);
                orthogonality_residual=partial_gram.norm()/std::sqrt(double(k+1));
                active_orthogonality=orthogonality_residual;
                if(!std::isfinite(orthogonality_residual) ||
                   orthogonality_residual>scaleSafeUDTPartialQOrthogonalityLimit)
                    scaleSafeQRGuardFail("nonunitary staged partial Q after normalization",n,k,
                        lost_bits,exponent_span,orthogonality_residual);
            }
#endif
            for(int j=k+1;j<n;++j) {
                DataType c1=q.col(k).dot(work.col(j)); work.col(j)-=q.col(k)*c1;
                DataType c2=q.col(k).dot(work.col(j)); work.col(j)-=q.col(k)*c2;
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
                const long long shift_wide=static_cast<long long>(exponent(j))-out.Dexp(k);
                if(shift_wide<std::numeric_limits<int>::min() ||
                   shift_wide>std::numeric_limits<int>::max())
                    scaleSafeQRGuardFail("T exponent bookkeeping overflow",n,k,lost_bits,
                        exponent_span,orthogonality_residual);
                const int shift=static_cast<int>(shift_wide);
#else
                const int shift=exponent(j)-out.Dexp(k);
#endif
                t(k,permutation[j])=scaleSafeLdexp((c1+c2)/diag_m,shift);
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
                const DataType coefficient=t(k,permutation[j]);
                if(!std::isfinite(coefficient.real()) || !std::isfinite(coefficient.imag()))
                    scaleSafeQRGuardFail("nonfinite T bookkeeping",n,k,lost_bits,
                        exponent_span,orthogonality_residual);
#endif
                renormalize(j,-1.0);
            }
        }
#ifndef PFQMC_UDT_DISABLE_RANK_LOSS_GUARD
        const MatType final_gram=q.adjoint()*q-MatType::Identity(n,n);
        const double final_orthogonality_residual=final_gram.norm()/std::sqrt(double(n));
        if(!std::isfinite(final_orthogonality_residual) ||
           final_orthogonality_residual>scaleSafeUDTPartialQOrthogonalityLimit)
            scaleSafeQRGuardFail("nonunitary final Q before factorization return",n,n-1,
                guard.max_lost_bits,exponent_span,final_orthogonality_residual);
#endif
        out.U=q; out.T=t;
        if(!scaleSafeFinite(out.U)||!scaleSafeFinite(out.T))
            throw std::runtime_error("scale-safe exponent QR produced nonfinite factors");
        return out;
    }

    // CPQR for X*diag(mantissa*2^exponent), without materializing the powers.
    static UDT factorColumnScaled(const MatType &X, const UDT &right)
    {
        const int n=right.nDim;
        if(X.rows()!=n || X.cols()!=n || !scaleSafeFinite(X))
            throw std::runtime_error("scale-safe exponent QR invalid base matrix");
        if(right.materializable()) {
            MatType dense=X;
            for(int j=0;j<n;++j) dense.col(j)*=right.actualD(j);
            UDT out(dense); return out;
        }
        MatType work=X;
        for(int j=0;j<n;++j) work.col(j)*=right.D(j);
        return factorPreparedScaled(work,right.Dexp);
    }

    // CPQR for diag(left.D)*X*diag(right.D), with both exponent vectors kept
    // symbolic.  Each input column is normalized by an exact power of two
    // before any floating-point QR operation.
    static UDT factorBiScaled(const UDT &left, const MatType &X, const UDT &right)
    {
        const int n=left.nDim;
        if(right.nDim!=n || X.rows()!=n || X.cols()!=n || !scaleSafeFinite(X))
            throw std::runtime_error("scale-safe two-sided QR invalid input");
        MatType work=MatType::Zero(n,n);
        iVecType exponent(n);
        for(int j=0;j<n;++j) {
            bool found=false;
            int column_exponent=std::numeric_limits<int>::min();
            for(int i=0;i<n;++i) {
                const DataType z=X(i,j)*left.D(i)*right.D(j);
                const double magnitude=std::abs(z);
                if(magnitude==0.0) continue;
                if(!std::isfinite(magnitude))
                    throw std::runtime_error("scale-safe two-sided QR nonfinite mantissa");
                int local_exponent=0;
                std::frexp(magnitude,&local_exponent);
                column_exponent=std::max(column_exponent,
                    left.Dexp(i)+right.Dexp(j)+local_exponent);
                found=true;
            }
            if(!found) throw std::runtime_error("scale-safe two-sided QR zero column");
            exponent(j)=column_exponent;
            for(int i=0;i<n;++i) {
                const DataType z=X(i,j)*left.D(i)*right.D(j);
                work(i,j)=scaleSafeLdexp(z,left.Dexp(i)+right.Dexp(j)-column_exponent);
            }
        }
        return factorPreparedScaled(work,exponent);
    }

    inline void onePlusInv(MatType &g) const
    {
        const MatType I=MatType::Identity(nDim,nDim);
        MatType Xinv=scaleSafeCheckedSolve(T,I,"UDT::onePlusInv/T");
        dVecType dp(nDim),dm(nDim);
        for(int i=0;i<nDim;++i) { dp(i)=dLargeInverse(i); dm(i)=dSmallPart(i); }
        MatType left=Xinv*dp.asDiagonal();
        MatType core=left+U*dm.asDiagonal();
        g=2.0*scaleSafeCheckedRightSolve(left,core,"UDT::onePlusInv/core");
    }
};

inline UDT operator*(const MatType &B, const UDT &right)
{
    MatType base=B*right.U;
    UDT out=UDT::factorColumnScaled(base,right);
    out.T=out.T*right.T;
    if(!scaleSafeFinite(out.T)) throw std::runtime_error("scale-safe MatType*UDT nonfinite T");
    return out;
}

inline UDT operator*(const UDT &left, const UDT &right)
{
    MatType bridge=left.T*right.U;
    UDT out=UDT::factorBiScaled(left,bridge,right);
    out.U=left.U*out.U;
    out.T=out.T*right.T;
    if(!scaleSafeFinite(out.U)||!scaleSafeFinite(out.T))
        throw std::runtime_error("scale-safe UDT*UDT produced nonfinite factors");
    return out;
}

inline MatType onePlusInv(UDT &left, UDT &right)
{
    const int n=right.nDim;
    MatType overlap=right.U.adjoint()*left.U;
    MatType bridge=right.T*left.T.adjoint();
    dVecType rp(n),rm(n),lp(n),lm(n);
    for(int i=0;i<n;++i) {
        rp(i)=right.dLargeInverse(i); rm(i)=right.dSmallPart(i);
        lp(i)=left.dLargeInverse(i); lm(i)=left.dSmallPart(i);
    }
    MatType core=rp.asDiagonal()*overlap*lp.asDiagonal()+rm.asDiagonal()*bridge*lm.asDiagonal();
    MatType rhs=rp.asDiagonal()*right.U.adjoint();
    MatType solved=scaleSafeCheckedSolve(core,rhs,"udtRdivudtL/core");
    MatType out=2.0*left.U*lp.asDiagonal()*solved;
    if(!scaleSafeFinite(out)) throw std::runtime_error("scale-safe udtRdivudtL nonfinite output");
    return out;
}

#else

class UDT
{
public:
    int nDim;
    MatType U;
    dVecType D;
    MatType T;
    UDT() = default;
    UDT(int nDim)
    {
        this->nDim = nDim;
        U = MatType::Identity(nDim, nDim);
        D = dVecType::Ones(nDim);
        T = MatType::Identity(nDim, nDim);
    }

    UDT(const MatType &_U, const dVecType &_D, const MatType &_T)
    {
        this->nDim = _U.rows();
        U = _U;
        D = _D;
        T = _T;
    }

    UDT &operator=(const UDT &other)
    {
        nDim = other.nDim;
        U = other.U;
        D = other.D;
        T = other.T;
        return *this;
    }

    UDT(const UDT &other)
    {
        *this = other;
    }

    UDT &operator=(UDT &&other)
    {
        if (this != &other)
        {
            nDim = other.nDim;
            U = std::move(other.U);
            D = std::move(other.D);
            T = std::move(other.T);
        }
        return (*this);
    }

    UDT(UDT &&other)
    {
        *this = std::move(other);
    }

    // use qr to get UDT decomposition
    explicit UDT(MatType &A)
    {
        nDim = A.rows();
        T = MatType::Zero(nDim, nDim);
        D = dVecType(nDim);
        int jpvt[nDim];
        DataType tau[nDim];
        for (int j = 0; j < nDim; j++)
        {
            jpvt[j] = 0;
        }

        LAPACKE_zgeqp3(LAPACK_COL_MAJOR, nDim, nDim, A.data(), nDim, jpvt, tau);

        double alpha;
        for (int i = 0; i < nDim; i++)
        {
            D(i) = std::abs(A(i, i).real()); // A(i, i)'s are all real
            alpha = 1.0 / D(i);
            for (int j = i; j < nDim; j++)
            {
                A(i, j) = A(i, j) * alpha;
            }
        }

        int j;
        for (int i = 0; i < nDim; i++)
        {
            j = jpvt[i] - 1;
            T(i, j) = 1;
        }

        cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, nDim, nDim, &one, A.data(), nDim, T.data(), nDim);

        LAPACKE_zungqr(LAPACK_COL_MAJOR, nDim, nDim, nDim, A.data(), nDim, tau);

        U = A;
    }

    // // F = (*this) * F
    // inline void factorizedMultUpdate(UDT &F)
    // {
    //     MatType mat = T * F.U;
    //     mat = D.asDiagonal() * mat;
    //     mat = mat * F.D.asDiagonal();
    //     UDT tmp(mat);
    //     F.U = U * tmp.U;
    //     F.T = tmp.T * F.T;
    //     F.D = tmp.D;
    // }

    // // (*this) = B * (*this)
    // inline void bMultUpdate(const MatType &B)
    // {
    //     // MatType r = ((B * U) * D.asDiagonal()) * T;
    //     MatType tmp = B * U;
    //     MatType tmp2 = tmp * D.asDiagonal();
    //     UDT F(tmp2);
    //     U = F.U;
    //     D = F.D;
    //     // std::cout << " r - UDTT " << (r - U * D.asDiagonal() * F.T * T).squaredNorm() << "\n";
    //     tmp = (F.T) * T;
    //     T = tmp;
    // }

    // // g = 2 * (1 + UDT)^{-1}
    // inline void onePlusInv(MatType &g) const
    // {
    //     MatType Xinv = T;
    //     MatType tmp1, tmp2;
    //     // std::cout << Tn << "\n === \n";
    //     int ipiv[nDim];
    //     LAPACKE_zgetrf(LAPACK_COL_MAJOR, nDim, nDim, Xinv.data(), nDim, ipiv);
    //     // std::cout << "ipiv= " << ipiv << "\n";
    //     // std::cout << Tn << "\n";
    //     LAPACKE_zgetri(LAPACK_COL_MAJOR, nDim, Xinv.data(), nDim, ipiv);

    //     // std::cout << "Xinv * X - I = " << (Xinv*T - MatType::Identity(nDim, nDim)).squaredNorm() << "\n";

    //     dVecType Dpinv(nDim);
    //     dVecType Dm(nDim);
    //     for (int i = 0; i < nDim; i++)
    //     {
    //         Dpinv(i) = 1.0 / std::max(D(i), 1.0);
    //         Dm(i) = std::min(D(i), 1.0);
    //     }

    //     tmp1 = Xinv * Dpinv.asDiagonal();
    //     tmp2 = U * Dm.asDiagonal();
    //     tmp1 = tmp1 + tmp2;

    //     UDT f = UDT(tmp1);
    //     LAPACKE_zgetrf(LAPACK_COL_MAJOR, nDim, nDim, f.T.data(), nDim, ipiv);
    //     LAPACKE_zgetri(LAPACK_COL_MAJOR, nDim, f.T.data(), nDim, ipiv);

    //     // std::cout << "Uinv * U - I = " << ((f.U)*(f.U.adjoint()) - MatType::Identity(nDim, nDim)).squaredNorm() << "\n";

    //     // MatType identity = MatType::Identity(nDim, nDim);
    //     // g = Xinv * Dpinv.asDiagonal() * f.T * (f.U * f.D.asDiagonal()).inverse();
    //     // return;

    //     tmp1 = (f.T) * (f.D.cwiseInverse()).asDiagonal();
    //     tmp2 = tmp1 * f.U.adjoint();
    //     tmp1 = Dpinv.asDiagonal() * tmp2;

    //     f = UDT(tmp1);

    //     f.D *= 2.0;

    //     tmp2 = Xinv * f.U;
    //     tmp1 = tmp2 * f.D.asDiagonal();
    //     g = tmp1 * f.T;
    // }


    // g = 2 * (1 + UDT)^{-1}
    inline void onePlusInv(MatType &g) const
    {
        MatType Xinv = T.inverse();
        dVecType Dpinv(nDim);
        dVecType Dm(nDim);
        for (int i = 0; i < nDim; i++)
        {
            Dpinv(i) = 1.0 / std::max(D(i), 1.0);
            Dm(i) = std::min(D(i), 1.0);
        }
        MatType tmp1 = Xinv * Dpinv.asDiagonal();
        MatType tmp2 = (tmp1 + U * Dm.asDiagonal()).inverse();
        g = 2.0 * tmp1 * tmp2;
    }
};

inline UDT operator*(const UDT &udtL, const UDT &udtR)
{
    MatType mat = udtL.T * udtR.U;
    mat = udtL.D.asDiagonal() * mat;
    mat = mat * udtR.D.asDiagonal();
    UDT tmp(mat);
    tmp.U = udtL.U * tmp.U;
    tmp.T = tmp.T * udtR.T;
    return tmp;
}

// TO DO: lazy evaluation, check udta = B * udta
inline UDT operator*(const MatType &B, const UDT &udtR)
{
    MatType mat = B * udtR.U;
    mat = mat * udtR.D.asDiagonal();
    UDT tmp(mat);
    tmp.T = tmp.T * udtR.T;
    return tmp;
}

// return (1+udtR@(udtL).adjoint)^{-1}
inline MatType onePlusInv(UDT &udtL, UDT &udtR)
{
    int n = udtR.U.cols();
    MatType tem1 = udtR.U.adjoint() * udtL.U;
    MatType tem2 = udtR.T * udtL.T.adjoint();
    auto DrPinv = dVecType(n);
    auto DrM = dVecType(n);
    auto DlPinv = dVecType(n);
    auto DlM = dVecType(n);
    for (int i = 0; i < n; i++)
    {
        if (std::abs(udtR.D(i)) > 1.0)
        {
            DrPinv(i) = 1.0 / udtR.D(i);
            DrM(i) = 1.0;
        }
        else
        {
            DrPinv(i) = 1.0;
            DrM(i) = udtR.D(i);
        }
        if (std::abs(udtL.D(i)) > 1.0)
        {
            DlPinv(i) = 1.0 / udtL.D(i);
            DlM(i) = 1.0;
        }
        else
        {
            DlPinv(i) = 1.0;
            DlM(i) = udtL.D(i);
        }
    }

    tem2 = DrPinv.asDiagonal() * tem1 * DlPinv.asDiagonal() + DrM.asDiagonal() * tem2 * DlM.asDiagonal();
    tem1 = DlPinv.asDiagonal() * tem2.inverse() * DrPinv.asDiagonal();

    return 2. * udtL.U * tem1 * (udtR.U.adjoint());
}
#endif // PFQMC_SCALE_SAFE_UDT
#endif
