#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"
#include "projector_json.h"

namespace {
struct Args { int L,boundary,hs,seed,burn,measurements,threads,diag_stride; double theta,beta,dt,V,delta,mu; std::string csv; };
Args parse(int n,char**v){
 if(n!=16) throw std::runtime_error("usage: static_note_driver L theta beta_trial dt V delta mu boundary hs seed burn measurements threads raw.csv diagnostic_stride");
 Args a{std::stoi(v[1]),std::stoi(v[8]),std::stoi(v[9]),std::stoi(v[10]),std::stoi(v[11]),std::stoi(v[12]),std::stoi(v[13]),std::stoi(v[15]),std::stod(v[2]),std::stod(v[3]),std::stod(v[4]),std::stod(v[5]),std::stod(v[6]),std::stod(v[7]),v[14]};
 long long nt=std::llround(a.theta/a.dt);
 if(a.L<2||a.theta<=0||a.beta<=0||a.dt<=0||a.V<0||(a.boundary!=0&&a.boundary!=1)||(a.hs!=0&&a.hs!=1)||a.burn<0||a.measurements<=0||a.threads<=0||a.diag_stride<0||a.csv.empty()||std::abs(nt*a.dt-a.theta)>1e-10) throw std::runtime_error("invalid parameters");
 return a;
}
class Walker:public Spinless_tV{public:int center=-1,ntrial=0,nphys=0;Walker(const SpinlessTvChainUtils*c,rdGenerator*r,double t,double b){build_projector_static_contour(*this,c,r,t,b,center,ntrial,nphys);}};
std::vector<int> fields(const Spinless_tV&w){std::vector<int>r;for(Operator*o:w.op_array)if(auto*s=o->getAuxField())for(int i=0;i<s->size();++i)r.push_back((*s)(i));return r;}
long long changes(const std::vector<int>&a,const std::vector<int>&b){long long n=0;for(size_t i=0;i<a.size();++i)n+=a[i]!=b[i];return n;}
double se(const std::vector<double>&x){if(x.size()<2)return 0;double m=0,q=0;for(double y:x)m+=y;m/=x.size();for(double y:x)q+=(y-m)*(y-m);return std::sqrt(q/(x.size()*(x.size()-1.)));}
double phys_spi(const SpinlessTvChainUtils&c,const MatType&g){return -c.StructureFactorCDW(g).real();}
double phys_sdq(const SpinlessTvChainUtils&c,const MatType&g){return -c.StructureFactorCDWOffset(g).real();}
}

int main(int argc,char**argv)try{
 Args a=parse(argc,argv);mkl_set_num_threads(a.threads);auto started=std::chrono::steady_clock::now();
 int ns=2*int(std::llround(a.theta/a.dt));SpinlessTvChainUtils cfg(a.L,a.dt,a.V,ns,a.boundary,a.delta,a.mu,a.hs);rdGenerator rd(a.seed);Walker w(&cfg,&rd,a.theta,a.beta);PfQMC q(&w,10);
 constexpr double guard_threshold=0.1;q.configureAdaptiveGuard(false,guard_threshold,100.0);
 for(int i=0;i<a.burn;++i){q.rightSweep();q.leftSweep();}
 std::ofstream raw(a.csv);if(!raw)throw std::runtime_error("cannot open raw CSV");
 raw<<"measurement,sign,sign_S_pi_numerator,sign_S_pi_dq_numerator,S_pi,S_pi_dq,R_cdw,acceptance,diag_relative_frobenius,diag_S_pi_abs_diff,diag_R_cdw_abs_diff,diag_sign_mismatch\n"<<std::setprecision(17);
 int nb=std::min(20,a.measurements);std::vector<double>bs(nb),bp(nb),bd(nb),brs,brp,brd,brr;std::vector<int>bn(nb);double ss=0,sp=0,sd=0,maxsi=0,maxoi=0;long long accepted=0,attempted=0,neg=0,recomp=0,corr=0,diag_n=0,diag_sm=0;double maxgf=0,maxdp=0,maxdr=0;
 for(int k=0;k<a.measurements;++k){
  if(k%20==0){DataType z=q.getSignRaw();++recomp;if(std::abs(q.sign-z)>1e-2){q.sign=z.real()>=0?DataType(1,0):DataType(-1,0);++corr;}}
  auto before=fields(w);MatType g;DataType sg;q.rightSweep(w.center,&g,&sg);auto after=fields(w);accepted+=changes(before,after);attempted+=before.size();
  DataType lp=cfg.StructureFactorCDW(g),ld=cfg.StructureFactorCDWOffset(g);double p=-lp.real(),d=-ld.real(),r=1-d/p,s=sg.real()>=0?1.:-1.;ss+=s;sp+=s*p;sd+=s*d;neg+=s<0;int b=std::min(nb-1,int((static_cast<long long>(k)*nb)/a.measurements));bs[b]+=s;bp[b]+=s*p;bd[b]+=s*d;++bn[b];maxsi=std::max(maxsi,std::abs(sg.imag()));maxoi=std::max({maxoi,std::abs(lp.imag()),std::abs(ld.imag())});
  double rel=std::numeric_limits<double>::quiet_NaN(),dp=rel,dr=rel;int sm=-1;
  if(a.diag_stride>0&&k%a.diag_stride==0){MatType full;q.rebuildGreenFromFullContourAtBoundary(w.center,full);double pf=phys_spi(cfg,full),df=phys_sdq(cfg,full),rf=1-df/pf;rel=(g-full).norm()/std::max(full.norm(),std::numeric_limits<double>::min());dp=std::abs(p-pf);dr=std::abs(r-rf);/* q.sign and getSignRaw() now refer to the same live HS configuration at the completed right-sweep boundary. */DataType sr=q.getSignRaw();sm=((q.sign.real()>=0)!=(sr.real()>=0));++diag_n;diag_sm+=sm;maxgf=std::max(maxgf,rel);maxdp=std::max(maxdp,dp);maxdr=std::max(maxdr,dr);}
  raw<<k<<','<<s<<','<<s*p<<','<<s*d<<','<<p<<','<<d<<','<<r<<','<<(attempted?double(accepted)/attempted:0)<<','<<rel<<','<<dp<<','<<dr<<','<<sm<<'\n';
  before.swap(after);q.leftSweep();after=fields(w);accepted+=changes(before,after);attempted+=before.size();
 }
#ifdef PFQMC_TEST_FORCE_ZERO_AVERAGE_SIGN
 ss=0.0;
#endif
 const bool sign_resolved=std::abs(ss)>=1e-12;
 const double unresolved=std::numeric_limits<double>::quiet_NaN();
 for(int b=0;b<nb;++b)if(bn[b]){
  brs.push_back(bs[b]/bn[b]);
  if(std::abs(bs[b])<=1e-12)continue;
  double pb=bp[b]/bs[b],db=bd[b]/bs[b];
  if(!std::isfinite(pb)||!std::isfinite(db)||std::abs(pb)<1e-14)continue;
  brp.push_back(pb);brd.push_back(db);brr.push_back(1-db/pb);
 }
 const double p=sign_resolved?sp/ss:unresolved,d=sign_resolved?sd/ss:unresolved;
 const bool observable_resolved=sign_resolved&&std::isfinite(p)&&std::isfinite(d)&&
                                std::abs(p)>=1e-14&&brp.size()>=2;
 const double r=observable_resolved?1-d/p:unresolved;
 const char* observable_status=observable_resolved?"resolved":
  (sign_resolved?"unresolved_insufficient_sign_reweighted_bins":
                 "unresolved_zero_average_sign");
 const double run=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
 const double gf=q.proposal_attempt_count?double(q.pre_decision_rebuild_count)/q.proposal_attempt_count:0;
 std::cout<<std::setprecision(17)<<"{\"status\":\"complete\",\"mode\":\"static_note_projector\",\"L\":"<<a.L<<",\"theta\":"<<a.theta<<",\"beta_trial\":"<<a.beta<<",\"dt\":"<<a.dt<<",\"V\":"<<a.V<<",\"delta\":"<<a.delta<<",\"mu\":"<<a.mu<<",\"boundary\":"<<a.boundary<<",\"hs_scheme\":"<<a.hs<<",\"seed\":"<<a.seed<<",\"burn\":"<<a.burn<<",\"measurements\":"<<a.measurements<<",\"physical_slices\":"<<ns<<",\"trial_slices\":"<<w.ntrial<<",\"observable_convention\":\"physical_density_SQ_equals_minus_legacy_contact_already_included\",\"sign_reweighted_observables_status\":\""<<observable_status<<"\",\"S_pi\":";
 projectorJsonNumber(std::cout,p,observable_resolved);
 std::cout<<",\"S_pi_err\":";projectorJsonNumber(std::cout,se(brp),observable_resolved);
 std::cout<<",\"S_pi_dq\":";projectorJsonNumber(std::cout,d,observable_resolved);
 std::cout<<",\"S_pi_dq_err\":";projectorJsonNumber(std::cout,se(brd),observable_resolved);
 std::cout<<",\"R_cdw\":";projectorJsonNumber(std::cout,r,observable_resolved);
 std::cout<<",\"R_cdw_err\":";projectorJsonNumber(std::cout,se(brr),observable_resolved);
 std::cout<<",\"average_sign\":";projectorJsonNumber(std::cout,ss/a.measurements);
 std::cout<<",\"average_sign_err\":";projectorJsonNumber(std::cout,se(brs));
 std::cout<<",\"acceptance\":";projectorJsonNumber(std::cout,double(accepted)/attempted);
 std::cout<<",\"runtime_seconds\":";projectorJsonNumber(std::cout,run);
 std::cout<<",\"negative_signs\":"<<neg<<",\"sign_recomputes\":"<<recomp<<",\"sign_corrections\":"<<corr<<",\"max_sign_imag\":";
 projectorJsonNumber(std::cout,maxsi);
 std::cout<<",\"max_observable_imag\":";projectorJsonNumber(std::cout,maxoi);
 std::cout<<",\"adaptive_guard\":false,\"guard_threshold\":"<<guard_threshold<<",\"proposal_attempt_count\":"<<q.proposal_attempt_count<<",\"min_update_denominator\":";
 projectorJsonNumber(std::cout,q.min_update_denominator);
 std::cout<<",\"guard_trigger_fraction\":";projectorJsonNumber(std::cout,gf);
 std::cout<<",\"adaptive_rebuild_count\":"<<q.adaptive_rebuild_count<<",\"pre_decision_rebuild_count\":"<<q.pre_decision_rebuild_count<<",\"post_accept_rebuild_count\":"<<q.post_accept_rebuild_count<<",\"full_rebuild_count\":"<<q.adaptive_rebuild_count<<",\"multiprecision_fallback\":false,\"multiprecision_fallback_count\":"<<q.multiprecision_fallback_count<<",\"diagnostic_stride\":"<<a.diag_stride<<",\"diagnostic_comparisons\":"<<diag_n<<",\"diagnostic_relative_frobenius_max\":";
 projectorJsonNumber(std::cout,maxgf,diag_n>0);
 std::cout<<",\"diagnostic_S_pi_abs_diff_max\":";projectorJsonNumber(std::cout,maxdp,diag_n>0);
 std::cout<<",\"diagnostic_R_cdw_abs_diff_max\":";projectorJsonNumber(std::cout,maxdr,diag_n>0);
 std::cout<<",\"diagnostic_sign_mismatch_count\":"<<diag_sm;
 projectorJsonBuildProvenance(std::cout,q);
 std::cout<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
