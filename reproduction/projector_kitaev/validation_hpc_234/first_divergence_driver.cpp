#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"

struct Walker: public Spinless_tV { int center=-1,ntrial=0,nphys=0; Walker(const SpinlessTvChainUtils*c,rdGenerator*r,double t,double b){build_projector_static_contour(*this,c,r,t,b,center,ntrial,nphys);} };
struct Stop { bool hit=false; std::string why; };
static double relerr(const MatType&a,const MatType&b){return (a-b).norm()/std::max(b.norm(),std::numeric_limits<double>::min());}
static int sgn(DataType z){return z.real()>=0?1:-1;}
static void emit(std::ofstream&o,const char*dir,int sweep,int l,Operator*op,int aux,int boundary,double pre,double post,DataType tracked0,DataType direct0,DataType fast,DataType full,double den,double u,bool af,bool afull,DataType tracked1,DataType direct1,const std::string&event){
 o<<dir<<','<<sweep<<','<<l<<','<<op->getType()<<','<<aux<<','<<boundary<<','<<std::setprecision(17)<<pre<<','<<post<<','<<tracked0.real()<<','<<tracked0.imag()<<','<<direct0.real()<<','<<direct0.imag()<<','<<fast.real()<<','<<fast.imag()<<','<<full.real()<<','<<full.imag()<<','<<den<<','<<u<<','<<af<<','<<afull<<','<<tracked1.real()<<','<<tracked1.imag()<<','<<direct1.real()<<','<<direct1.imag()<<','<<event<<'\n';
}
static bool flips(PfQMC&q,Operator*op,int l,int boundary,const char*dir,int sweep,std::ofstream&o,Stop&stop){
 int n=op->singleFlipProposalCount(); if(n<=0){q.sign*=op->update(q.g);return false;}
 for(int aux=0;aux<n;++aux){
  MatType fullPre; q.rebuildGreenFromFullContourAtBoundary(boundary,fullPre); double pre=relerr(q.g,fullPre); DataType t0=q.sign,d0=q.getSignRaw();
  double u=0; if(!op->prepareSingleFlip(q.g,aux,&u)) throw std::runtime_error("proposal unavailable"); DataType fast=op->preparedRatio(); double den=op->preparedMinDenominator(); DataType full=op->recomputePreparedRatio(fullPre); bool af=u<std::abs(fast), afull=u<std::abs(full);
  std::string event; if(pre>1e-6)event="GREEN_PRE"; else if(sgn(t0)!=sgn(d0))event="SIGN_PRE"; else if(std::abs(fast-full)>1e-6*std::max(1.0,std::abs(full)))event="RATIO"; else if(af!=afull)event="ACCEPTANCE";
  if(!event.empty()){emit(o,dir,sweep,l,op,aux,boundary,pre,std::numeric_limits<double>::quiet_NaN(),t0,d0,fast,full,den,u,af,afull,t0,d0,event);stop.hit=true;stop.why=event;return true;}
  DataType delta=op->finishSingleFlip(q.g,af,true); q.sign*=delta; MatType fullPost; q.rebuildGreenFromFullContourAtBoundary(boundary,fullPost); double post=relerr(q.g,fullPost); DataType d1=q.getSignRaw(); std::string after; if(post>1e-6)after="GREEN_POST"; else if(sgn(q.sign)!=sgn(d1))after="SIGN_POST";
  emit(o,dir,sweep,l,op,aux,boundary,pre,post,t0,d0,fast,full,den,u,af,afull,q.sign,d1,after);
  if(!after.empty()){stop.hit=true;stop.why=after;return true;}
 }
 return false;
}
static bool right(PfQMC&q,int sweep,std::ofstream&o,Stop&stop){MatType tmp=MatType::Identity(q.nDim,q.nDim),A=MatType::Identity(q.nDim,q.nDim);int seg=0;for(int l=0;l<q.op_length;++l){Operator*op=q.op_array[l];if(flips(q,op,l,l,"right",sweep,o,stop))return true;op->left_multiply(A,tmp);std::swap(A,tmp);if(q.need_stabilization[(l+1)%q.op_length]){if(seg==0)q.udtR[seg]=UDT(A);else q.udtR[seg]=A*q.udtR[seg-1];A=MatType::Identity(q.nDim,q.nDim);if(seg==q.checkpoints-1)q.udtR[seg].onePlusInv(q.g);else q.g=onePlusInv(q.udtL[seg+1],q.udtR[seg]);++seg;}else op->left_propagate(q.g,tmp);}return false;}
static bool left(PfQMC&q,int sweep,std::ofstream&o,Stop&stop){MatType tmp=MatType::Identity(q.nDim,q.nDim),A=MatType::Identity(q.nDim,q.nDim);int seg=q.checkpoints-1;for(int l=q.op_length-1;l>=0;--l){Operator*op=q.op_array[l];op->right_propagate(q.g,tmp);if(flips(q,op,l,(l+1)%q.op_length,"left",sweep,o,stop))return true;op->right_multiply(A,tmp);std::swap(A,tmp);if(q.need_stabilization[l]){A.adjointInPlace();if(seg==q.checkpoints-1)q.udtL[seg]=UDT(A);else q.udtL[seg]=A*q.udtL[seg+1];A=MatType::Identity(q.nDim,q.nDim);if(seg==0){q.udtL[seg].onePlusInv(q.g);q.g.adjointInPlace();}else q.g=onePlusInv(q.udtL[seg],q.udtR[seg-1]);--seg;}}return false;}
int main(int argc,char**argv)try{
 if(argc!=3)throw std::runtime_error("usage: driver seed csv");int seed=std::stoi(argv[1]);std::ofstream o(argv[2]);o<<"direction,sweep,operator,bond,aux,boundary,green_rel_pre,green_rel_post,tracked_pre_re,tracked_pre_im,direct_pre_re,direct_pre_im,fast_ratio_re,fast_ratio_im,full_ratio_re,full_ratio_im,min_denominator,uniform,fast_accept,full_accept,tracked_post_re,tracked_post_im,direct_post_re,direct_post_im,event\n";
 SpinlessTvChainUtils cfg(10,.1,4,200,0,1,0,0);rdGenerator rd(seed);Walker w(&cfg,&rd,10,8);PfQMC q(&w,10);q.configureAdaptiveGuard(false,.1,100.);Stop s;
 for(int sw=0;sw<3&&!s.hit;++sw){right(q,sw,o,s);if(!s.hit)left(q,sw,o,s);} std::cout<<"{\"seed\":"<<seed<<",\"event\":\""<<s.why<<"\",\"stopped\":"<<(s.hit?"true":"false")<<"}\n";return s.hit?0:1;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
