#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "multiprecision_driven_rebuild.h"

struct DebugWalker:Spinless_tV{
 int center=-1,ntrial=0,ninit=0,nramp=0;std::vector<std::string> region;
 DebugWalker(const SpinlessTvChainUtils*c,rdGenerator*rd,double V0,double Vf,double rate,double theta,double beta){
  nDim=c->nDim;ninit=llround(theta/c->dt);nramp=llround((Vf-V0)/(rate*c->dt));int nb[2]={(c->Lx+1)/2,c->Lx/2};
  MatType h(nDim,nDim);c->KineticGenerator(h);
  for(double rem=beta;rem>1e-12;){double st=std::min(c->dt,rem);MatType x=h,b=expm(x,-st);x=st*h;op_array.push_back(new DenseOperator(b,signOfHamiltonian(x)));region.push_back("trial");++ntrial;rem-=st;}
  MatType x=h,kh=expm(x,-c->dt/2);x=c->dt*h/2;DataType sk=signOfHamiltonian(x);
  auto add=[&](double V,bool dag,const std::string&r){op_array.push_back(new DenseOperator(kh,sk));region.push_back(r);for(int z=0;z<2;++z){int bond=dag?1-z:z;auto*s=new iVecType(nb[bond]);for(int j=0;j<s->size();++j)(*s)(j)=rd->rdZ2();op_array.push_back(new SpinlessVOperator(c,s,bond,rd,V));region.push_back(r);}op_array.push_back(new DenseOperator(kh,sk));region.push_back(r);};
  for(int l=0;l<ninit;++l)add(V0,false,"ket_init");
  for(int l=0;l<nramp;++l)add(V0+rate*(l+.5)*c->dt,false,l<nramp/3?"ket_early":l<2*nramp/3?"ket_middle":"ket_late");
  center=op_array.size();
  for(int l=nramp-1;l>=0;--l)add(V0+rate*(l+.5)*c->dt,true,l>=2*nramp/3?"bra_late":l>=nramp/3?"bra_middle":"bra_early");
  for(int l=0;l<ninit;++l)add(V0,true,"bra_init");
 }
};

MatType green_at(DebugWalker&w,int target){UDT f(w.nDim);for(int j=0;j<(int)w.op_array.size();++j)w.op_array[(target+j)%w.op_array.size()]->stabilizedLeftMultiply(f);MatType g;f.onePlusInv(g);return g;}

struct UltraProxy { double d_spread, core_condition, solve_residual; };
UltraProxy ultra_proxy_at(DebugWalker&w,int target){
 UDT f(w.nDim);for(int j=0;j<(int)w.op_array.size();++j)w.op_array[(target+j)%w.op_array.size()]->stabilizedLeftMultiply(f);
 double dmin=f.D.minCoeff(),dmax=f.D.maxCoeff();MatType xinv=f.T.inverse();dVecType dpi(f.nDim),dm(f.nDim);
 for(int i=0;i<f.nDim;++i){dpi(i)=1./std::max(f.D(i),1.);dm(i)=std::min(f.D(i),1.);}
 MatType lhs=xinv*dpi.asDiagonal(),core=lhs+f.U*dm.asDiagonal(),g;f.onePlusInv(g);
 Eigen::JacobiSVD<MatType> svd(core);auto sv=svd.singularValues();
 double residual=(g*core-2.*lhs).norm()/std::max(2.*lhs.norm(),std::numeric_limits<double>::min());
 return {dmax/std::max(dmin,std::numeric_limits<double>::min()),sv(0)/std::max(sv(sv.size()-1),std::numeric_limits<double>::min()),residual};
}

std::complex<double> logdet_full(DebugWalker&w){UDT f(w.nDim);for(auto*o:w.op_array)o->stabilizedLeftMultiply(f);dVecType dp(f.nDim),di(f.nDim),dm(f.nDim);for(int i=0;i<f.nDim;++i){dp(i)=std::max(f.D(i),1.);di(i)=1./dp(i);dm(i)=std::min(f.D(i),1.);}MatType m=f.T.inverse()*di.asDiagonal()+f.U*dm.asDiagonal();Eigen::PartialPivLU<MatType>lu(m);auto d=lu.matrixLU().diagonal();std::complex<double> z=lu.permutationP().determinant()<0?std::complex<double>(0,M_PI):0.;for(int i=0;i<d.size();++i)z+=std::log(d(i))+std::log(dp(i));return z;}

DataType ratio_fast(const SpinlessTvChainUtils&c,const SpinlessVOperator&o,const MatType&g,int ia){int a,b,d,e;c.aux2MajoranaIdx(ia,0,o.bondType,a,b);c.aux2MajoranaIdx(ia,1,o.bondType,d,e);double s=(*o.s)(ia);DataType t0=1.-DataType(0,1)*o.thlV*s*g(a,b),t1=1.-DataType(0,1)*o.thlV*s*g(d,e);return o.etaM*(t0*t1+o.thlV*o.thlV*(g(a,d)*g(b,e)-g(b,d)*g(a,e)));}

void complex_json(std::ostream&f,DataType z){f<<'['<<z.real()<<','<<z.imag()<<']';}
void matrix_json(std::ostream&f,const MatType&m){f<<'[';for(int i=0;i<m.rows();++i)for(int j=0;j<m.cols();++j){if(i||j)f<<',';complex_json(f,m(i,j));}f<<']';}
void dump_snapshot(const std::string&dir,int k,int oi,int ia,const DebugWalker&w,const MatType&gf,const MatType&gfull,DataType rf,DataType rfull){mkdir(dir.c_str(),0777);std::ofstream f(dir+"/flip_"+std::to_string(k)+".json");f<<std::setprecision(17)<<"{\"flip\":"<<k<<",\"operator_index\":"<<oi<<",\"aux_index\":"<<ia<<",\"region\":\""<<w.region[oi]<<"\",\"R_fast\":";complex_json(f,rf);f<<",\"R_full\":";complex_json(f,rfull);f<<",\"G_fast\":";matrix_json(f,gf);f<<",\"G_full\":";matrix_json(f,gfull);f<<",\"operators\":[";bool first=true;for(int j=0;j<(int)w.op_array.size();++j)if(auto*o=dynamic_cast<SpinlessVOperator*>(w.op_array[j])){if(!first)f<<',';first=false;f<<"{\"index\":"<<j<<",\"local_V\":"<<o->localV<<",\"bond\":"<<o->bondType<<",\"s\":[";for(int q=0;q<o->s->size();++q){if(q)f<<',';f<<(*o->s)(q);}f<<"]}";}f<<"]}\n";}
void rebuild_B(const SpinlessTvChainUtils&c,SpinlessVOperator&o){o.B=MatType::Identity(o.nDim,o.nDim);for(int q=0;q<o.s->size();++q){int a,b,c1,d;c.aux2MajoranaIdx(q,0,o.bondType,a,b);c.aux2MajoranaIdx(q,1,o.bondType,c1,d);for(auto p:{std::pair<int,int>{a,b},std::pair<int,int>{c1,d}}){o.B(p.first,p.first)=o.chlV;o.B(p.second,p.second)=o.chlV;o.B(p.first,p.second)=DataType(0,1)*o.shlV*double((*o.s)(q));o.B(p.second,p.first)=-DataType(0,1)*o.shlV*double((*o.s)(q));}}}

int main(int argc,char**argv){
 if(argc!=6&&argc!=7&&argc!=8){std::cerr<<"usage: driven_fastupdate_check Vf dt flips stb trace.csv [guard threshold | snapshot_dir]\n";return 2;}
 double Vf=std::stod(argv[1]),dt=std::stod(argv[2]);int flips=std::stoi(argv[3]),stb=std::stoi(argv[4]);
 int L=6;if(const char* value=std::getenv("PFQMC_FASTCHECK_L"))L=std::stoi(value);
 double rate=1.;if(const char* value=std::getenv("PFQMC_FASTCHECK_RATE"))rate=std::stod(value);
 bool guard=argc==8&&std::string(argv[6])=="guard";double threshold=guard?std::stod(argv[7]):0.1;
 std::string snap=argc==7?argv[6]:"";
 SpinlessTvChainUtils c(L,dt,0,1,1,1,0,0);rdGenerator rd(424242);DebugWalker w(&c,&rd,0,Vf,rate,6,8);PfQMC rebuilder(&w,stb);
 std::ofstream tr(argv[5]);tr<<"flip,operator_index,region,local_V,bond,aux,min_denominator,R_before,R_used,ratio_abs_error,ratio_complex_error,green_error,guard_triggered,udt_d_spread,udt_core_condition,udt_solve_residual,stabilized,error_before_reset,error_after_reset\n"<<std::setprecision(17);
 double maxr=0,maxabs=0,maxg=0,maxpre=0,maxpost=0,minDen=INFINITY;long long triggers=0;int k=0,oi=0,ops_since_reset=0;MatType g=green_at(w,0),tmp;
 while(k<flips){
  auto*o=dynamic_cast<SpinlessVOperator*>(w.op_array[oi]);
  if(o){
   int ia=k%o->s->size();DataType rbefore=ratio_fast(c,*o,g,ia),rused=rbefore;MatType gpre=g;
   double ignored=0;o->prepareSingleFlip(g,ia,&ignored);double md=o->preparedMinDenominator();minDen=std::min(minDen,md);
   bool dangerous=guard&&(md<threshold||std::abs(rbefore)>100.0);auto ld0=logdet_full(w);PfQMC q0(&w,10);DataType s0=q0.getSignRaw();
   if(dangerous){++triggers;rebuilder.rebuildGreenFromFullContourAtBoundary(oi,g);rused=o->recomputePreparedRatio(g);o->finishSingleFlip(g,true,false);rebuilder.rebuildGreenFromFullContourAtBoundary(oi,g);double cp=rebuilder.fullContourCoreConditionAtBoundary(oi);if(cp>1e4){MatType mg;if(driven_multiprecision::rebuild(c,w.op_array,oi,w.ntrial,mg))g.swap(mg);if(k==254){std::ofstream f(std::string(argv[5])+".flip254_mp.json");f<<std::setprecision(17)<<"{\"green\":";matrix_json(f,g);f<<",\"core_condition\":"<<cp<<"}\n";}}}
   else {bool accepted=false;DataType phase=1;o->singleFlip(g,ia,0.,accepted,phase);}
   auto ld1=logdet_full(w);PfQMC q1(&w,10);DataType s1=q1.getSignRaw();double mag=std::exp(.5*(ld1.real()-ld0.real()));DataType rfull=s1/s0;rfull=rfull/std::abs(rfull)*mag;MatType gf=green_at(w,oi);
   if(!snap.empty()&&((k>=164&&k<=170)||(k>=251&&k<=257)||(Vf<4.5&&k>=160&&k<=166))){(*o->s)(ia)=-(*o->s)(ia);rebuild_B(c,*o);dump_snapshot(snap,k,oi,ia,w,gpre,green_at(w,oi),rbefore,rfull);(*o->s)(ia)=-(*o->s)(ia);rebuild_B(c,*o);}
   double er=std::abs(rused-rfull),ea=std::abs(std::abs(rused)-std::abs(rfull)),eg=(g-gf).cwiseAbs().maxCoeff();maxr=std::max(maxr,er);maxabs=std::max(maxabs,ea);maxg=std::max(maxg,eg);UltraProxy up=ultra_proxy_at(w,oi);
   tr<<k<<','<<oi<<','<<w.region[oi]<<','<<o->localV<<','<<o->bondType<<','<<ia<<','<<md<<','<<std::abs(rbefore)<<','<<std::abs(rused)<<','<<ea<<','<<er<<','<<eg<<','<<dangerous<<','<<up.d_spread<<','<<up.core_condition<<','<<up.solve_residual<<",0,0,0\n";++k;
  }
  w.op_array[oi]->left_propagate(g,tmp);oi=(oi+1)%w.op_array.size();++ops_since_reset;
  if(ops_since_reset>=stb){MatType gf=green_at(w,oi);double pre=(g-gf).cwiseAbs().maxCoeff();maxpre=std::max(maxpre,pre);g=gf;double post=(g-gf).cwiseAbs().maxCoeff();maxpost=std::max(maxpost,post);tr<<k<<','<<oi<<",reset,0,-1,-1,0,0,0,0,0,0,0,0,0,0,1,"<<pre<<','<<post<<'\n';ops_since_reset=0;}
 }
 std::cout<<std::setprecision(17)<<"{\"Vf\":"<<Vf<<",\"dt\":"<<dt<<",\"flips\":"<<flips<<",\"stb\":"<<stb<<",\"guard\":"<<(guard?"true":"false")<<",\"guard_threshold\":"<<threshold<<",\"guard_triggers\":"<<triggers<<",\"guard_frequency\":"<<double(triggers)/flips<<",\"min_denominator\":"<<minDen<<",\"max_ratio_complex_abs_error\":"<<maxr<<",\"max_ratio_magnitude_abs_error\":"<<maxabs<<",\"max_green_abs_error\":"<<maxg<<",\"max_error_before_reset\":"<<maxpre<<",\"max_error_after_reset\":"<<maxpost<<"}\n";return 0;
}
