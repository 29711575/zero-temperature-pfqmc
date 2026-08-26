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
struct Args { int L,boundary,hs,seed,burn,measurements,threads,diag_stride,stb,guard,sign_stride; double theta,beta,dt,V,delta,mu; std::string csv; };
Args parse(int n,char**v){
 if(n!=19) throw std::runtime_error("usage: regression_driver L theta beta_trial dt V delta mu boundary hs seed burn measurements threads raw.csv diagnostic_stride stabilization_interval guard sign_stride");
 Args a{std::stoi(v[1]),std::stoi(v[8]),std::stoi(v[9]),std::stoi(v[10]),std::stoi(v[11]),std::stoi(v[12]),std::stoi(v[13]),std::stoi(v[15]),std::stoi(v[16]),std::stoi(v[17]),std::stoi(v[18]),std::stod(v[2]),std::stod(v[3]),std::stod(v[4]),std::stod(v[5]),std::stod(v[6]),std::stod(v[7]),v[14]};
 long long nt=std::llround(a.theta/a.dt);
 if(a.L<2||a.theta<=0||a.beta<=0||a.dt<=0||a.V<0||(a.boundary!=0&&a.boundary!=1)||(a.hs!=0&&a.hs!=1)||a.burn<0||a.measurements<=0||a.threads<=0||a.diag_stride<0||a.stb<=0||(a.guard!=0&&a.guard!=1)||a.sign_stride<=0||a.csv.empty()||std::abs(nt*a.dt-a.theta)>1e-10) throw std::runtime_error("invalid parameters");
 return a;
}
class Walker:public Spinless_tV{public:int center=-1,ntrial=0,nphys=0;Walker(const SpinlessTvChainUtils*c,rdGenerator*r,double t,double b){build_projector_static_contour(*this,c,r,t,b,center,ntrial,nphys);}};
std::vector<int> fields(const Spinless_tV&w){std::vector<int>r;for(Operator*o:w.op_array)if(auto*s=o->getAuxField())for(int i=0;i<s->size();++i)r.push_back((*s)(i));return r;}
long long changes(const std::vector<int>&a,const std::vector<int>&b){long long n=0;for(size_t i=0;i<a.size();++i)n+=a[i]!=b[i];return n;}
double se(const std::vector<double>&x){if(x.size()<2)return 0;double m=0,q=0;for(double y:x)m+=y;m/=x.size();for(double y:x)q+=(y-m)*(y-m);return std::sqrt(q/(x.size()*(x.size()-1.)));}
double phys_spi(const SpinlessTvChainUtils&c,const MatType&g){return -c.StructureFactorCDW(g).real();}
double phys_sdq(const SpinlessTvChainUtils&c,const MatType&g){return -c.StructureFactorCDWOffset(g).real();}
double percentile(std::vector<double> x,double p){if(x.empty())return 0;std::sort(x.begin(),x.end());double at=p*(x.size()-1),lo=std::floor(at),hi=std::ceil(at),f=at-lo;return x[size_t(lo)]*(1-f)+x[size_t(hi)]*f;}
}

int main(int argc,char**argv)try{
 Args a=parse(argc,argv);mkl_set_num_threads(a.threads);resetScaleSafeQRGuardDiagnostics();auto started=std::chrono::steady_clock::now();
 int ns=2*int(std::llround(a.theta/a.dt));SpinlessTvChainUtils cfg(a.L,a.dt,a.V,ns,a.boundary,a.delta,a.mu,a.hs);rdGenerator rd(a.seed);Walker w(&cfg,&rd,a.theta,a.beta);PfQMC q(&w,a.stb);
 constexpr double guard_threshold=0.1;q.configureAdaptiveGuard(a.guard!=0,guard_threshold,100.0);
 for(int i=0;i<a.burn;++i){q.rightSweep();q.leftSweep();}
 std::ofstream raw(a.csv);if(!raw)throw std::runtime_error("cannot open raw CSV");
 raw<<"measurement,sign,sign_S_pi_numerator,sign_S_pi_dq_numerator,S_pi,S_pi_dq,R_cdw,acceptance,diag_relative_frobenius,diag_S_pi_abs_diff,diag_R_cdw_abs_diff,diag_sign_mismatch\n"<<std::setprecision(17);
 int nb=std::min(20,a.measurements);std::vector<double>bs(nb),bp(nb),bd(nb),brs,brp,brd,brr;std::vector<int>bn(nb);double ss=0,sp=0,sd=0,maxsi=0,maxoi=0,center_green_norm=0;long long accepted=0,attempted=0,neg=0,recomp=0,corr=0,diag_n=0,diag_sm=0;double maxgf=0,maxdp=0,maxdr=0;ProjectorRawSignChecks rawChecks;
 for(int k=0;k<a.measurements;++k){
  if(k%a.sign_stride==0){const PfaffianResult z=q.getSignRawWithStatus();++recomp;rawChecks.record(q.sign,z);}
  auto before=fields(w);MatType g;DataType sg;q.rightSweep(w.center,&g,&sg);auto after=fields(w);accepted+=changes(before,after);attempted+=before.size();
  DataType lp=cfg.StructureFactorCDW(g),ld=cfg.StructureFactorCDWOffset(g);double p=-lp.real(),d=-ld.real(),r=1-d/p,s=sg.real()>=0?1.:-1.;if(k==0)center_green_norm=g.norm();ss+=s;sp+=s*p;sd+=s*d;neg+=s<0;int b=std::min(nb-1,int((static_cast<long long>(k)*nb)/a.measurements));bs[b]+=s;bp[b]+=s*p;bd[b]+=s*d;++bn[b];maxsi=std::max(maxsi,std::abs(sg.imag()));maxoi=std::max({maxoi,std::abs(lp.imag()),std::abs(ld.imag())});
  double rel=std::numeric_limits<double>::quiet_NaN(),dp=rel,dr=rel;int sm=-1;
  if(a.diag_stride>0&&k%a.diag_stride==0){MatType full;q.rebuildGreenFromFullContourAtBoundary(0,full);double pfast=phys_spi(cfg,q.g),dfast=phys_sdq(cfg,q.g),rfast=1-dfast/pfast,pf=phys_spi(cfg,full),df=phys_sdq(cfg,full),rf=1-df/pf;rel=(q.g-full).norm()/std::max(full.norm(),std::numeric_limits<double>::min());dp=std::abs(pfast-pf);dr=std::abs(rfast-rf);const PfaffianResult sr=q.getSignRawWithStatus();if(sr.ok()){sm=((q.sign.real()>=0)!=(sr.value.real()>=0));diag_sm+=sm;}++diag_n;maxgf=std::max(maxgf,rel);maxdp=std::max(maxdp,dp);maxdr=std::max(maxdr,dr);}
  raw<<k<<','<<s<<','<<s*p<<','<<s*d<<','<<p<<','<<d<<','<<r<<','<<(attempted?double(accepted)/attempted:0)<<','<<rel<<','<<dp<<','<<dr<<','<<sm<<'\n';
  before.swap(after);q.leftSweep();after=fields(w);accepted+=changes(before,after);attempted+=before.size();
 }
 if(std::abs(ss)<1e-12)throw std::runtime_error("zero average sign");
 for(int b=0;b<nb;++b)if(bn[b]&&std::abs(bs[b])>1e-12){double p=bp[b]/bs[b],d=bd[b]/bs[b];brs.push_back(bs[b]/bn[b]);brp.push_back(p);brd.push_back(d);brr.push_back(1-d/p);}
 double p=sp/ss,d=sd/ss,r=1-d/p,run=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();double gf=q.proposal_attempt_count?double(q.pre_decision_rebuild_count)/q.proposal_attempt_count:0;
 const ScaleSafeQRGuardDiagnostics &udt=scaleSafeQRGuardDiagnostics();
 std::cout<<std::setprecision(17)<<"{\"status\":\"complete\",\"mode\":\"static_regression_projector\",\"L\":"<<a.L<<",\"theta\":"<<a.theta<<",\"beta_trial\":"<<a.beta<<",\"dt\":"<<a.dt<<",\"V\":"<<a.V<<",\"delta\":"<<a.delta<<",\"mu\":"<<a.mu<<",\"boundary\":"<<a.boundary<<",\"hs_scheme\":"<<a.hs<<",\"seed\":"<<a.seed<<",\"burn\":"<<a.burn<<",\"measurements\":"<<a.measurements<<",\"stabilization_interval\":"<<a.stb<<",\"sign_recompute_stride\":"<<a.sign_stride<<",\"physical_slices\":"<<ns<<",\"trial_slices\":"<<w.ntrial<<",\"center_green_norm_first\":";
 projectorJsonNumber(std::cout,center_green_norm);
 std::cout<<",\"S_pi\":";projectorJsonNumber(std::cout,p);
 std::cout<<",\"S_pi_dq\":";projectorJsonNumber(std::cout,d);
 std::cout<<",\"R_cdw\":";projectorJsonNumber(std::cout,r);
 std::cout<<",\"average_sign\":";projectorJsonNumber(std::cout,ss/a.measurements);
 std::cout<<",\"acceptance\":";projectorJsonNumber(std::cout,double(accepted)/attempted);
 std::cout<<",\"runtime_seconds\":";projectorJsonNumber(std::cout,run);
 std::cout<<",\"sign_corrections\":"<<corr<<",\"max_sign_imag\":";projectorJsonNumber(std::cout,maxsi);
 std::cout<<",\"max_observable_imag\":";projectorJsonNumber(std::cout,maxoi);
 std::cout<<",\"diagnostic_comparisons\":"<<diag_n<<",\"diagnostic_relative_frobenius_max\":";
 projectorJsonNumber(std::cout,maxgf,diag_n>0);
 std::cout<<",\"diagnostic_sign_mismatch_count\":"<<diag_sm<<",\"udt_guard_bits\":"<<scaleSafeUDTRankLossGuardBits<<",\"udt_orth_precheck_bits\":"<<scaleSafeUDTOrthogonalityPrecheckBits<<",\"udt_guard_triggers\":"<<udt.trigger_count<<",\"udt_max_lost_bits\":";
 projectorJsonNumber(std::cout,udt.max_lost_bits);
 std::cout<<",\"udt_min_guard_margin\":";projectorJsonNumber(std::cout,udt.min_guard_margin);
 projectorJsonRawSignChecks(std::cout,rawChecks);
 projectorJsonBuildProvenance(std::cout,q);
 std::cout<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
