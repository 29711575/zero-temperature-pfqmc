#ifndef Spinless_tV_H
#define Spinless_tV_H

#include "operator.h"
#include "types.h"
#include <algorithm>
#include <limits>

class SpinlessTvUtils {
   public:
    int Lx, Ly;
    // int nsites;
    int nDim;
    // std::vector<int> nBond;
    double dt;
    double V;
    int l;
    double lambdaV, chlV, shlV, thlV, etaM;

    const int hsScheme;  

    const bool singleMaj;

    SpinlessTvUtils(int _Lx, int _Ly, double _dt, double _V, int _l, int _nDim,
                    bool _singleMaj = false, int _hsScheme = 0)
        : singleMaj(_singleMaj), hsScheme(_hsScheme) {
        // model configuration
        Lx = _Lx;
        Ly = _Ly;
        dt = _dt;
        V = _V;
        l = _l;  // imaginary time slices
        nDim = _nDim;

        lambdaV = acosh(exp(0.5 * V * dt));
        chlV = cosh(lambdaV);
        shlV = sinh(lambdaV);
        thlV = tanh(lambdaV);
        etaM = chlV * chlV;
    }

    inline int unitCellCoord2Idx(int ix, int iy) const { return ix * Ly + iy; }

    // virtual inline void KineticGenerator(MatType &H, DataType t) const {};

    // virtual inline DataType energyFromGreensFunc(const MatType &g) {};
    virtual inline void aux2MajoranaIdx(int idxAux, int imaj, int bType,
                                        int &idx1, int &idx2) const = 0;

    // Generate the Greens function for single slice
    inline void InteractionTanhGenerator(MatType &H, const iVecType &s,
                                         const int bondType,
                                         bool inv = false) const {
        DataType tmp = (std::complex<double>(0.0, 1.0)) * tanh(0.5 * lambdaV);
        if (inv) {
            tmp = (-1.0) / tmp;
        }

        if (hsScheme == 0) {
            int idx1, idx2;
            if (!singleMaj) {
                for (int i = 0; i < s.size(); i++) {
                    for (int k = 0; k < 2; k++) {
                        aux2MajoranaIdx(i, k, bondType, idx1, idx2);
                        H(idx1, idx2) += -tmp * double(s(i));
                        H(idx2, idx1) += +tmp * double(s(i));
                    }
                }
            } else {
                for (int i = 0; i < s.size(); i++) {
                    aux2MajoranaIdx(i, 0, bondType, idx1, idx2);
                    H(idx1, idx2) += -tmp * double(s(i));
                    H(idx2, idx1) += +tmp * double(s(i));
                }
            }
        } else if (hsScheme == 1) {
            int idxi1, idxi2, idxj1, idxj2;
            assert(!singleMaj);
            for (int i = 0; i < s.size(); i++) {
                aux2MajoranaIdx(i, 0, bondType, idxi1, idxj1);
                aux2MajoranaIdx(i, 1, bondType, idxi2, idxj2);
                H(idxi1, idxi2) += -tmp * double(s(i));
                H(idxi2, idxi1) += +tmp * double(s(i));
                H(idxj1, idxj2) += +tmp * double(s(i));
                H(idxj2, idxj1) += -tmp * double(s(i));
            }
        }
    }

    // Directly generate B by directly writing each 2*2 block
    // B should be initialized as Identity
    inline void InteractionBGenerator(MatType &B, const iVecType &s,
                                      const int bondType,
                                      bool inv = false) const {
        DataType ch = chlV;
        DataType ish = (std::complex<double>(0.0, 1.0)) * shlV;
        if (inv) ish = -ish;

        if (hsScheme == 0) {
            // B = exp(λ/2 * s * (i γi1 γj1 + i γi2 γj2))
            int idx1, idx2;
            if (!singleMaj) {
                for (int i = 0; i < s.size(); i++) {
                    for (int k = 0; k < 2; k++) {
                        aux2MajoranaIdx(i, k, bondType, idx1, idx2);
                        B(idx1, idx1) = ch;
                        B(idx2, idx2) = ch;
                        B(idx1, idx2) = +ish * double(s(i));
                        B(idx2, idx1) = -ish * double(s(i));
                    }
                }
            } else {
                for (int i = 0; i < s.size(); i++) {
                    aux2MajoranaIdx(i, 0, bondType, idx1, idx2);
                    B(idx1, idx1) = ch;
                    B(idx2, idx2) = ch;
                    B(idx1, idx2) = +ish * double(s(i));
                    B(idx2, idx1) = -ish * double(s(i));
                }
            }
        } else if (hsScheme == 1) {
            // B = exp(λ/2 * s * (i γi1 γi2 - i γj1 γj2))
            int idxi1, idxi2, idxj1, idxj2;
            assert(!singleMaj);
            for (int i = 0; i < s.size(); i++) {
                aux2MajoranaIdx(i, 0, bondType, idxi1, idxj1);
                aux2MajoranaIdx(i, 1, bondType, idxi2, idxj2);
                B(idxi1, idxi1) = ch;
                B(idxi2, idxi2) = ch;
                B(idxi1, idxi2) = +ish * double(s(i));
                B(idxi2, idxi1) = -ish * double(s(i));

                B(idxj1, idxj1) = ch;
                B(idxj2, idxj2) = ch;
                B(idxj1, idxj2) = -ish * double(s(i));
                B(idxj2, idxj1) = +ish * double(s(i));
            }
        }
    }
};

class SpinlessVOperator : public Operator {
   protected:
    const SpinlessTvUtils *config;
    void interactionBGenerator(MatType &out, bool inv) const {
        DataType ch = chlV;
        DataType ish = std::complex<double>(0.0, 1.0) * shlV;
        if (inv) ish = -ish;
        int a, b, c, d;
        for (int i = 0; i < s->size(); ++i) {
            config->aux2MajoranaIdx(i, 0, bondType, a, b);
            config->aux2MajoranaIdx(i, 1, bondType, c, d);
            if (hsScheme == 0) {
                for (const auto &p : {std::pair<int,int>{a,b}, std::pair<int,int>{c,d}}) {
                    out(p.first,p.first)=ch; out(p.second,p.second)=ch;
                    out(p.first,p.second)=+ish*double((*s)(i));
                    out(p.second,p.first)=-ish*double((*s)(i));
                }
            } else {
                out(a,a)=ch; out(c,c)=ch;
                out(a,c)=+ish*double((*s)(i)); out(c,a)=-ish*double((*s)(i));
                out(b,b)=ch; out(d,d)=ch;
                out(b,d)=-ish*double((*s)(i)); out(d,b)=+ish*double((*s)(i));
            }
        }
    }
    void interactionTanhGenerator(MatType &out, bool inv) const {
        DataType tmp=std::complex<double>(0.0,1.0)*std::tanh(0.5*lambdaV);
        if(inv) tmp=-1.0/tmp;
        int a,b,c,d;
        for(int i=0;i<s->size();++i){
            config->aux2MajoranaIdx(i,0,bondType,a,b);
            config->aux2MajoranaIdx(i,1,bondType,c,d);
            if(hsScheme==0){
                out(a,b)+=-tmp*double((*s)(i)); out(b,a)+=tmp*double((*s)(i));
                out(c,d)+=-tmp*double((*s)(i)); out(d,c)+=tmp*double((*s)(i));
            } else {
                out(a,c)+=-tmp*double((*s)(i)); out(c,a)+=tmp*double((*s)(i));
                out(b,d)+=tmp*double((*s)(i)); out(d,b)+=-tmp*double((*s)(i));
            }
        }
    }
    // delayed update for spinless t-V requires additional diagonalization
    // const int delay_max = 32;
    // cannot be larger than 64 (unless you change the threads value in
    // kernel_linalg.cpp)
   public:
    // Diagnostics for numerical conditioning of the rank-2 updates.  These
    // do not participate in the update math or production decisions.
    DataType last_denom1=DataType(1), last_denom2=DataType(1);
    bool pendingValid=false;
    int pendingAux=-1, pendingOldSigma=0, pendingNewSigma=0;
    DataType pendingRatio=DataType(1), pendingDenom1=DataType(1), pendingDenom2=DataType(1);
    double pendingUniform=0.0;
    int pendingCursor=0;
    const double localV, lambdaV, chlV, shlV, thlV, etaM;
    // const int nUnitcell;
    const int bondType;
    // const int Naux; // number of auxillary fields (live on bonds)
    const int nDim;  // dimension of hamiltonian, number of sites × 2 (number of
                     // majorana species)

    const int hsScheme; // 0: hopping channel, 1: density channel
    iVecType *s;
    MatType B;
    MatType B_inv;
    rdGenerator *rd;

    const bool singleMaj;

    // _s: aux fields, Z_2 variable, length = nUnitcell
    SpinlessVOperator(const SpinlessTvUtils *_config, iVecType *_s,
                      int _bondType, rdGenerator *_rd)
        : config(_config),
          localV(_config->V), lambdaV(_config->lambdaV), chlV(_config->chlV),
          shlV(_config->shlV), thlV(_config->thlV), etaM(_config->etaM),
          bondType(_bondType),
          nDim(config->nDim),
          singleMaj(_config->singleMaj),
          hsScheme(_config->hsScheme) {
        s = _s;
        rd = _rd;
        B = MatType::Identity(nDim, nDim);
        interactionBGenerator(B, false);
    }

    SpinlessVOperator(const SpinlessTvUtils *_config, iVecType *_s,
                      int _bondType, rdGenerator *_rd, double _localV)
        : config(_config), localV(_localV),
          lambdaV(std::acosh(std::exp(0.5*_localV*_config->dt))),
          chlV(std::cosh(lambdaV)), shlV(std::sinh(lambdaV)),
          thlV(std::tanh(lambdaV)), etaM(chlV*chlV),
          bondType(_bondType), nDim(config->nDim),
          singleMaj(_config->singleMaj), hsScheme(_config->hsScheme) {
        s=_s; rd=_rd; B=MatType::Identity(nDim,nDim);
        interactionBGenerator(B,false);
    }

    ~SpinlessVOperator() { delete s; }

    void reCalcInv() {
        B_inv = MatType::Identity(nDim, nDim);
        interactionBGenerator(B_inv, true);
    }

    // virtual inline void aux2MajoranaIdx(int idxAux, int imaj, int& idx1, int&
    // idx2) {};

    void singleFlip(MatType &g, int idxAux, double rand, bool &flag,
                    DataType &signCur) {
        DataType r;
        // auto m = mConfig->idxCell2Coord(idxCell);
        int auxCur = (*s)(idxAux);
        int idx1, idx2, idx3, idx4;
        DataType tmp[2];
        const int inc = 1;
        DataType alpha;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2); // 1j, 1k
        config->aux2MajoranaIdx(idxAux, 1, bondType, idx3, idx4); // 2j, 2k

        if (hsScheme == 0) {
            tmp[0] =
                (1.0 - ((std::complex<double>(0.0, 1.0)) * thlV * double(auxCur) * g(idx1, idx2)));
            tmp[1] =
                (1.0 - ((std::complex<double>(0.0, 1.0)) * thlV * double(auxCur) * g(idx3, idx4)));
            r = tmp[0] * tmp[1];
            // std::cout << "r1" << r << " ";
            r += (thlV * thlV) * ((g(idx1, idx3) * g(idx2, idx4)) -
                                                (g(idx2, idx3) * g(idx1, idx4)));
            // std::cout << " r2" << r << "\n";
            r *= etaM;
        } else if (hsScheme == 1) {
            tmp[0] =
                (1.0 - ((std::complex<double>(0.0, 1.0)) * thlV * double(auxCur) * g(idx1, idx3)));
            tmp[1] =
                (1.0 + ((std::complex<double>(0.0, 1.0)) * thlV * double(auxCur) * g(idx2, idx4)));
            r = tmp[0] * tmp[1];
            // std::cout << "r1" << r << " ";
            r -= (thlV * thlV) * ((g(idx1, idx2) * g(idx3, idx4)) -
                                                (g(idx3, idx2) * g(idx1, idx4)));
            // std::cout << " r2" << r << "\n";
            r *= etaM;
        }
        last_denom1=tmp[0]; last_denom2=tmp[1];
        // for (int imaj = 0; imaj < 2; imaj ++) {
        //     config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
        //     // idx1 = mConfig->majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
        //     // idx2 = mConfig->neighborSiteIdx(m.ix, m.iy, imaj, bondType);
        //     // tmp = [1 + i \sigma_{12} \tanh(\lambda / 2) G_{12}]
        //     tmp[imaj] = ( 1.0 - ( (std::complex<double>(0.0, 1.0)) * (config->thlV) * double(auxCur) *
        //     g(idx1, idx2) ) ); r *= tmp[imaj];
        // }

        flag = rand < std::abs(r);
        // std::cout << rand << "=rand " << "r = " << r << "\n";

        if (flag) {
            // DataType t = (r / std::abs(r));
            // if (std::abs(t.imag()) > 1e-5) {
            //     std::cout << "r= " << r << " sign(r)= " << t << "\n";
            // }
            signCur *= (r / std::abs(r));

            if (hsScheme == 0) {
                for (int imaj = 0; imaj < 2; imaj++) {
                    config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
                    if (imaj == 1) {
                        tmp[1] = (1.0 - ((std::complex<double>(0.0, 1.0)) * thlV * double(auxCur) *
                                        g(idx1, idx2)));
                    }
                    // update aux field and B matrix
                    (*s)(idxAux) = -auxCur;
                    B(idx1, idx2) = -B(idx1, idx2);
                    B(idx2, idx1) = -B(idx2, idx1);

                    // update Green's function
                    cVecType x1 = -g.col(idx1);
                    cVecType x2 = -g.col(idx2);
                    x1(idx1) += 2;
                    x2(idx2) += 2;
                    alpha = (+std::complex<double>(0.0, 1.0)) * double(auxCur) * thlV / tmp[imaj];
                    zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                        g.data(), &nDim);
                    alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                        g.data(), &nDim);
                }
            } else if (hsScheme == 1) {
                int idxj1, idxk1, idxj2, idxk2;
                config->aux2MajoranaIdx(idxAux, 0, bondType, idxj1, idxk1);
                config->aux2MajoranaIdx(idxAux, 1, bondType, idxj2, idxk2);
                for (int iaux = 0; iaux < 2; iaux++) {
                    if (iaux == 0) {
                        idx1 = idxj1;
                        idx2 = idxj2;
                    } else {
                        idx1 = idxk1;
                        idx2 = idxk2;
                        tmp[1] = (1.0 + ((std::complex<double>(0.0, 1.0)) * thlV * double(auxCur) *
                                        g(idx1, idx2)));
                    }

                    // update aux field and B matrix
                    (*s)(idxAux) = -auxCur;
                    B(idx1, idx2) = -B(idx1, idx2);
                    B(idx2, idx1) = -B(idx2, idx1);

                    // update Green's function
                    cVecType x1 = -g.col(idx1);
                    cVecType x2 = -g.col(idx2);
                    x1(idx1) += 2;
                    x2(idx2) += 2;
                    alpha = (+std::complex<double>(0.0, 1.0)) * double(auxCur) * thlV / tmp[iaux];
                    if (iaux == 1) alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                        g.data(), &nDim);
                    alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                        g.data(), &nDim);
                }
            }
        }
    };

    // Prepare one proposal without changing the auxiliary field, B, or G.
    // The random number is deliberately supplied by the caller so a rebuild
    // can recompute this exact proposal without advancing the RNG.
    bool prepareSingleFlip(MatType &g, double *u = nullptr) override {
        if (singleMaj || s->size()==0) return false;
        // Legacy update() traverses auxiliary fields in fixed order; it does
        // not draw a random auxiliary index.  Preserve that proposal order.
        const int aux = pendingCursor % s->size();
        pendingCursor = (pendingCursor + 1) % s->size();
        return prepareSingleFlip(g, aux, u);
    }

    bool prepareSingleFlip(MatType &g, int idxAux, double *u = nullptr) override {
        if (singleMaj || idxAux < 0 || idxAux >= s->size()) return false;
        pendingAux = idxAux;
        pendingOldSigma = (*s)(pendingAux);
        pendingNewSigma = -pendingOldSigma;
        int i1,i2,i3,i4;
        config->aux2MajoranaIdx(pendingAux,0,bondType,i1,i2);
        config->aux2MajoranaIdx(pendingAux,1,bondType,i3,i4);
        const DataType I(0.0,1.0);
        if (hsScheme==0) {
            pendingDenom1 = 1.0 - I*thlV*double(pendingOldSigma)*g(i1,i2);
            pendingDenom2 = 1.0 - I*thlV*double(pendingOldSigma)*g(i3,i4);
            pendingRatio = (pendingDenom1*pendingDenom2 + thlV*thlV*(g(i1,i3)*g(i2,i4)-g(i2,i3)*g(i1,i4)))*etaM;
        } else {
            pendingDenom1 = 1.0 - I*thlV*double(pendingOldSigma)*g(i1,i3);
            pendingDenom2 = 1.0 + I*thlV*double(pendingOldSigma)*g(i2,i4);
            pendingRatio = (pendingDenom1*pendingDenom2 - thlV*thlV*(g(i1,i2)*g(i3,i4)-g(i3,i2)*g(i1,i4)))*etaM;
        }
        last_denom1=pendingDenom1; last_denom2=pendingDenom2;
        pendingValid=true;
        pendingUniform=rd->rdUniform01();
        if (u) *u=pendingUniform;
        return true;
    }

    int singleFlipProposalCount() const override {
        return singleMaj ? 0 : s->size();
    }
    DataType preparedRatio() const override { return pendingRatio; }
    double preparedUniform() const override { return pendingUniform; }
    double preparedMinDenominator() const override {
        return std::min(std::abs(pendingDenom1), std::abs(pendingDenom2));
    }

    DataType recomputePreparedRatio(MatType &g) override {
        if (!pendingValid) return DataType(1);
        int saveAux=pendingAux; int saveOld=pendingOldSigma;
        // Re-evaluate using the saved identity; no RNG is consumed.
        int i1,i2,i3,i4; config->aux2MajoranaIdx(saveAux,0,bondType,i1,i2); config->aux2MajoranaIdx(saveAux,1,bondType,i3,i4);
        const DataType I(0.0,1.0); double q=double(saveOld);
        if (hsScheme==0) { pendingDenom1=1.0-I*thlV*q*g(i1,i2); pendingDenom2=1.0-I*thlV*q*g(i3,i4); pendingRatio=(pendingDenom1*pendingDenom2+thlV*thlV*(g(i1,i3)*g(i2,i4)-g(i2,i3)*g(i1,i4)))*etaM; }
        else { pendingDenom1=1.0-I*thlV*q*g(i1,i3); pendingDenom2=1.0+I*thlV*q*g(i2,i4); pendingRatio=(pendingDenom1*pendingDenom2-thlV*thlV*(g(i1,i2)*g(i3,i4)-g(i3,i2)*g(i1,i4)))*etaM; }
        last_denom1=pendingDenom1; last_denom2=pendingDenom2; return pendingRatio;
    }

    DataType finishSingleFlip(MatType &g, bool accept, bool use_fast_update=true) override {
        if (!pendingValid) return DataType(1);
        bool flag=false; DataType sg=1.0;
        // Force the saved proposal's decision; singleFlip performs the
        // established HS/B/rank-update algebra unchanged.
        if (use_fast_update) {
        singleFlip(g,pendingAux, accept ? 0.0 : std::numeric_limits<double>::infinity(), flag, sg);
        } else {
            // Apply the established field/B algebra while deliberately
            // discarding the rank-updated Green matrix.  The caller can then
            // install a stabilized rebuild for the accepted configuration.
            MatType scratch=g;
            singleFlip(scratch,pendingAux, accept ? 0.0 : std::numeric_limits<double>::infinity(), flag, sg);
        }
        pendingValid=false;
        return flag ? sg : DataType(1);
    }

    void singleFlipSingleMajorana(MatType &g, int idxAux, double rand,
                                  bool &flag, DataType &signCur) {
        DataType r;
        // auto m = mConfig->idxCell2Coord(idxCell);
        int auxCur = (*s)(idxAux);
        int idx1, idx2;
        DataType tmp;
        const int inc = 1;
        DataType alpha;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
        tmp =
            (1.0 - ((std::complex<double>(0.0, 1.0)) * thlV * double(auxCur) * g(idx1, idx2)));
        r = tmp * tmp * etaM;

        flag = rand < std::abs(r);
        // std::cout << rand << "=rand " << "r = " << r << "\n";

        if (flag) {
            // std::cout << "tmp = " << tmp << " r = " << r << "\n";
            signCur *= (tmp / std::abs(tmp));
            config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
            // update aux field and B matrix
            (*s)(idxAux) = -auxCur;
            B(idx1, idx2) = -B(idx1, idx2);
            B(idx2, idx1) = -B(idx2, idx1);

            // update Green's function
            cVecType x1 = -g.col(idx1);
            cVecType x2 = -g.col(idx2);
            x1(idx1) += 2;
            x2(idx2) += 2;
            alpha = (+std::complex<double>(0.0, 1.0)) * double(auxCur) * thlV / tmp;
            zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                  g.data(), &nDim);
            alpha = -alpha;
            zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                  g.data(), &nDim);
        }
    };

    DataType update(MatType &g) override {
        double rand;
        bool flag;
        DataType signCur = 1.0;
        if (singleMaj) {
            for (int i = 0; i < s->size(); i++) {
                // random real between (0, 1)
                rand = rd->rdUniform01();
                singleFlipSingleMajorana(g, i, rand, flag, signCur);
            }
        } else {
            for (int i = 0; i < s->size(); i++) {
                // random real between (0, 1)
                rand = rd->rdUniform01();
                singleFlip(g, i, rand, flag, signCur);
            }
        }
        return signCur;
    }

    void right_multiply(const MatType &AIn, MatType &AOut) override {
        AOut = AIn * B;
    }

    void left_multiply(const MatType &AIn, MatType &Aout) override {
        Aout = B * AIn;
    }

    void left_propagate(MatType &g, MatType &gTmp) override {
        reCalcInv();
        gTmp = B * g;
        g = gTmp * B_inv;
    }
    void right_propagate(MatType &g, MatType &gTmp) override {
        reCalcInv();
        gTmp = B_inv * g;
        g = gTmp * B;
    }

    iVecType *getAuxField() override { return s; }

    int getType() override { return bondType; }

    void getGreensMat(MatType &g) override {
        g = MatType::Zero(nDim, nDim);
        interactionTanhGenerator(g, false);
    }

    void getGreensMatInv(MatType &g) override {
        g = MatType::Zero(nDim, nDim);
        interactionTanhGenerator(g, true);
    }

    // inline DataType getSignPfGInv() override;

    void stabilizedLeftMultiply(UDT &F) override {
        // F.bMultUpdate(B);
        F = B * F;
    }
};

class Spinless_tV {
   public:
    std::vector<Operator *> op_array;
    int nDim;
    ~Spinless_tV() {
        for (int i = 0; i < op_array.size(); i++) {
            delete op_array[i];
        }
    }
};

#endif
