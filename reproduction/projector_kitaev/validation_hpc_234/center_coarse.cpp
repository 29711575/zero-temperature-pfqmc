#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"
struct W:public Spinless_tV{int center=-1,a=0,b=0;W(const SpinlessTvChainUtils*c,rdGenerator*r){build_projector_static_contour(*this,c,r,10,8,center,a,b);}};
static double re(const MatType&a,const MatType&b){return(a-b).norm()/std::max(b.norm(),std::numeric_limits<double>::min());}static int ss(DataType z){return z.real()>=0?1:-1;}
static bool right_center(PfQMC&q,int center,int iter,std::ofstream&o){MatType tmp=MatType::Identity(q.nDim,q.nDim),A=MatType::Identity(q.nDim,q.nDim);int seg=0;for(int l=0;l<q.op_length;++l){Operator*op=q.op_array[l];q.sign*=op->update(q.g);op->left_multiply(A,tmp);std::swap(A,tmp);if(q.need_stabilization[(l+1)%q.op_length]){if(seg==0)q.udtR[seg]=UDT(A);else q.udtR[seg]=A*q.udtR[seg-1];A.setIdentity();if(seg==q.checkpoints-1)q.udtR[seg].onePlusInv(q.g);else q.g=onePlusInv(q.udtL[seg+1],q.udtR[seg]);++seg;}else op->left_propagate(q.g,tmp);if(l+1==center){MatType f;q.rebuildGreenFromFullContourAtBoundary(center,f);DataType d=q.getSignRaw();double e=re(q.g,f);bool sm=ss(q.sign)!=ss(d);o<<iter<<','<<center<<','<<std::setprecision(17)<<q.sign.real()<<','<<q.sign.imag()<<','<<d.real()<<','<<d.imag()<<','<<sm<<','<<e<<'\n';if(sm||e>1e-6)return true;}}return false;}
int main(int ac,char**av){int seed=std::stoi(av[1]);std::ofstream o(av[2]);o<<"measurement,boundary,tracked_re,tracked_im,direct_re,direct_im,sign_mismatch,green_rel\n";SpinlessTvChainUtils c(6,.1,4,200,0,1,0,0);rdGenerator r(seed);W w(&c,&r);PfQMC q(&w,10);q.configureAdaptiveGuard(false,.1,100.);for(int i=0;i<500;++i){q.rightSweep();q.leftSweep();}int hit=-1;for(int i=0;i<5000;++i){if(right_center(q,w.center,i,o)){hit=i;break;}q.leftSweep();}std::ofstream j(av[3]);j<<"{\"seed\":"<<seed<<",\"first_center_event\":"<<hit<<"}\n";}
