#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"

namespace {
struct Args { int L, boundary, hs_scheme, seed, burn, measurements, threads; double theta, beta_trial, dt, V, delta, mu; };

Args parse_args(int argc, char **argv) {
    if (argc != 14) throw std::runtime_error("usage: projector_driver L theta beta_trial dt V delta mu boundary hs_scheme seed burn measurements threads");
    Args a{std::stoi(argv[1]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),
           std::stoi(argv[11]), std::stoi(argv[12]), std::stoi(argv[13]), std::stod(argv[2]),
           std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7])};
    const long long nt = std::llround(a.theta / a.dt);
    if (a.L < 2 || a.theta <= 0 || a.beta_trial <= 0 || a.dt <= 0 || a.V < 0 ||
        (a.boundary != 0 && a.boundary != 1) || (a.hs_scheme != 0 && a.hs_scheme != 1) ||
        a.burn < 0 || a.measurements <= 0 || a.threads <= 0 ||
        std::abs(nt * a.dt - a.theta) > 1e-10)
        throw std::runtime_error("invalid parameters (theta must be an integer multiple of dt)");
    return a;
}

class ProjectorKitaevWalker : public Spinless_tV {
  public:
    int center_boundary = -1, n_trial_slices = 0, n_physical_slices = 0;
    ProjectorKitaevWalker(const SpinlessTvChainUtils *c, rdGenerator *rd,
                          double theta, double beta_trial) {
        build_projector_static_contour(*this,c,rd,theta,beta_trial,center_boundary,n_trial_slices,n_physical_slices);
    }
};

std::vector<int> fields(const Spinless_tV &w) {
    std::vector<int> r;
    for (Operator *op : w.op_array) if (iVecType *s = op->getAuxField())
        for (int i = 0; i < s->size(); ++i) r.push_back((*s)(i));
    return r;
}
long long changes(const std::vector<int>& a, const std::vector<int>& b) {
    long long n=0; for (size_t i=0;i<a.size();++i) n += a[i]!=b[i]; return n;
}
double blocking_stderr(const std::vector<double>& x) {
    if (x.size() < 2) return 0.0;
    double mean=0.0; for(double v:x) mean+=v; mean/=x.size();
    double ss=0.0; for(double v:x) ss+=(v-mean)*(v-mean);
    return std::sqrt(ss/(double(x.size())*double(x.size()-1)));
}
}

int main(int argc, char **argv) try {
    const Args a = parse_args(argc, argv); mkl_set_num_threads(a.threads);
    const auto started = std::chrono::steady_clock::now();
    const int physical_slices = 2 * int(std::llround(a.theta / a.dt));
    SpinlessTvChainUtils config(a.L,a.dt,a.V,physical_slices,a.boundary,a.delta,a.mu,a.hs_scheme);
    rdGenerator random(a.seed); ProjectorKitaevWalker walker(&config,&random,a.theta,a.beta_trial);
    PfQMC qmc(&walker,10);
    for(int i=0;i<a.burn;++i){qmc.rightSweep();qmc.leftSweep();}
    const int n_bins=std::min(15,a.measurements);
    std::vector<double> bsign(n_bins,0.0),bspi(n_bins,0.0),bspidq(n_bins,0.0);
    std::vector<int> bcount(n_bins,0);
    double sign_sum=0, spi_sum=0, spidq_sum=0, max_si=0, max_oi=0;
    long long accepted=0,attempted=0; int negatives=0,recomputes=0,corrections=0;
    for(int sample=0;sample<a.measurements;++sample){
        if(sample%20==0){ DataType raw=qmc.getSignRaw(); ++recomputes;
            if(std::abs(qmc.sign-raw)>1e-2){qmc.sign=raw.real()>=0?DataType(1,0):DataType(-1,0);++corrections;}}
        auto before=fields(walker); MatType gc; DataType sc;
        qmc.rightSweep(walker.center_boundary,&gc,&sc); auto after=fields(walker);
        accepted+=changes(before,after); attempted+=before.size();
        DataType spi=config.StructureFactorCDW(gc), spidq=config.StructureFactorCDWOffset(gc);
        if(!std::isfinite(sc.real())||!std::isfinite(sc.imag())||!std::isfinite(spi.real())||!std::isfinite(spidq.real())) return 3;
        double s=sc.real()>=0?1.0:-1.0; sign_sum+=s; spi_sum+=s*spi.real(); spidq_sum+=s*spidq.real(); negatives+=s<0;
        const int bin=std::min(n_bins-1,int((static_cast<long long>(sample)*n_bins)/a.measurements));
        bsign[bin]+=s; bspi[bin]+=s*spi.real(); bspidq[bin]+=s*spidq.real(); ++bcount[bin];
        max_si=std::max(max_si,std::abs(sc.imag())); max_oi=std::max({max_oi,std::abs(spi.imag()),std::abs(spidq.imag())});
        before.swap(after); qmc.leftSweep(); after=fields(walker);
        accepted+=changes(before,after); attempted+=before.size();
    }
    if(std::abs(sign_sum)<1e-12) return 4;
    // PfQMC::g has unit diagonal, so legacy already contains the negative
    // onsite contact.  Only flip the overall sign; do not add contact again.
    const double onsite_contact=1.0/(4.0*a.L);
    const double spi=-spi_sum/sign_sum, spidq=-spidq_sum/sign_sum;
    const double spi_offsite=spi-onsite_contact, spidq_offsite=spidq-onsite_contact;
    if(std::abs(spi)<1e-14) return 4;
    std::vector<double> spi_offsite_bins,spidq_offsite_bins,spi_bins,spidq_bins,r_bins,sign_bins;
    for(int b=0;b<n_bins;++b){
        if(bcount[b]>0) sign_bins.push_back(bsign[b]/bcount[b]);
        if(std::abs(bsign[b])<1e-12) continue;
        const double x=-bspi[b]/bsign[b], y=-bspidq[b]/bsign[b];
        const double xoff=x-onsite_contact, yoff=y-onsite_contact;
        if(std::abs(x)<1e-14) continue;
        spi_offsite_bins.push_back(xoff); spidq_offsite_bins.push_back(yoff);
        spi_bins.push_back(x); spidq_bins.push_back(y); r_bins.push_back(1.0-y/x);
    }
    if(spi_bins.size()<2) return 4;
    const double spi_offsite_err=blocking_stderr(spi_offsite_bins),spidq_offsite_err=blocking_stderr(spidq_offsite_bins);
    const double spi_err=blocking_stderr(spi_bins),spidq_err=blocking_stderr(spidq_bins);
    const double r=1.0-spidq/spi,r_err=blocking_stderr(r_bins);
    const double sign_mean=sign_sum/a.measurements,sign_err=blocking_stderr(sign_bins);
    double runtime=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::cout<<std::setprecision(17)<<"{\"mode\":\"projector\",\"L\":"<<a.L<<",\"theta\":"<<a.theta
      <<",\"beta_trial\":"<<a.beta_trial<<",\"dt\":"<<a.dt<<",\"V\":"<<a.V<<",\"delta\":"<<a.delta
      <<",\"mu\":"<<a.mu<<",\"boundary\":"<<a.boundary<<",\"hs_scheme\":"<<a.hs_scheme<<",\"seed\":"<<a.seed
      <<",\"burn\":"<<a.burn<<",\"measurements\":"<<a.measurements<<",\"threads\":"<<a.threads
      <<",\"physical_slices\":"<<physical_slices<<",\"trial_slices\":"<<walker.n_trial_slices
      <<",\"center_operator_boundary\":"<<walker.center_boundary
      <<",\"observable_convention\":\"physical_density_SQ_equals_minus_legacy_contact_already_included\""
      <<",\"error_method\":\"contiguous_sign_reweighted_bins\",\"n_bins\":"<<n_bins
      <<",\"onsite_contact\":"<<onsite_contact<<",\"onsite_contact_is_diagnostic_only\":true"
      <<",\"S_pi_offsite\":"<<spi_offsite<<",\"S_pi_offsite_mean\":"<<spi_offsite<<",\"S_pi_offsite_err\":"<<spi_offsite_err
      <<",\"S_pi_dq_offsite\":"<<spidq_offsite<<",\"S_pi_dq_offsite_mean\":"<<spidq_offsite<<",\"S_pi_dq_offsite_err\":"<<spidq_offsite_err
      <<",\"S_pi\":"<<spi<<",\"S_pi_mean\":"<<spi<<",\"S_pi_err\":"<<spi_err
      <<",\"S_pi_dq\":"<<spidq<<",\"S_pi_dq_mean\":"<<spidq<<",\"S_pi_dq_err\":"<<spidq_err
      <<",\"R_cdw\":"<<r<<",\"R_cdw_mean\":"<<r<<",\"R_cdw_err\":"<<r_err
      <<",\"average_sign\":"<<sign_mean<<",\"average_sign_mean\":"<<sign_mean<<",\"average_sign_err\":"<<sign_err
      <<",\"acceptance\":"<<double(accepted)/attempted<<",\"runtime_seconds\":"<<runtime
      <<",\"negative_signs\":"<<negatives<<",\"sign_recomputes\":"<<recomputes<<",\"sign_corrections\":"<<corrections
      <<",\"max_sign_imag\":"<<max_si<<",\"max_observable_imag\":"<<max_oi<<"}\n";
    return 0;
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
