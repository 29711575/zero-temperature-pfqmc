#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "kitaevChain.h"
#include "pfqmc.h"
#include "projector_contour.h"
struct Walker:public Spinless_tV{Walker(const SpinlessTvChainUtils*c,rdGenerator*r){int x,y,z;build_projector_static_contour(*this,c,r,10,8,x,y,z);}};
static double rel(const MatType&a,const MatType&b){return(a-b).norm()/std::max(b.norm(),std::numeric_limits<double>::min());}
static int s(DataType z){return z.real()>=0?1:-1;}
static bool check(PfQMC&q,std::ofstream&o,const char*phase,int sw,int l,int b,int aux=-1){MatType f;q.rebuildGreenFromFullContourAtBoundary(b,f);DataType d=q.getSignRaw();double e=rel(q.g,f);o<<phase<<','<<sw<<','<<l<<','<<b<<','<<aux<<','<<std::setprecision(17)<<e<<','<<q.sign.real()<<','<<q.sign.imag()<<','<<d.real()<<','<<d.imag()<<"\n";return e>1e-6||s(q.sign)!=s(d);}
int main(int ac,char**av)try{
 if(ac!=3)throw std::runtime_error("usage seed csv");int seed=std::stoi(av[1]);std::ofstream o(av[2]);o<<"phase,sweep,operator,boundary,aux,green_rel,tracked_re,tracked_im,direct_re,direct_im\n";
 SpinlessTvChainUtils c(10,.1,4,200,0,1,0,0);rdGenerator r(seed);Walker w(&c,&r);PfQMC q(&w,10);q.configureAdaptiveGuard(false,.1,100.);q.rightSweep();
 MatType tmp=MatType::Identity(q.nDim,q.nDim),A=MatType::Identity(q.nDim,q.nDim);int seg=q.checkpoints-1;std::string why="";
 for(int l=q.op_length-1;l>=0;--l){int before=(l+1)%q.op_length;if(check(q,o,"left_pre_propagation",0,l,before)){why="PRE_PROPAGATION";break;}Operator*op=q.op_array[l];op->right_propagate(q.g,tmp);if(check(q,o,"left_post_propagation",0,l,l)){why="DIRECT_PROPAGATION";break;}
  int n=op->singleFlipProposalCount();for(int aux=0;aux<n;++aux){double u=0;op->prepareSingleFlip(q.g,aux,&u);DataType ratio=op->preparedRatio();bool accept=u<std::abs(ratio);DataType ds=op->finishSingleFlip(q.g,accept,true);q.sign*=ds;if(check(q,o,accept?"left_accepted_update":"left_rejected_update",0,l,l,aux)){why=accept?"ACCEPTED_LOCAL_UPDATE":"REJECTED_LOCAL_UPDATE";break;}} if(!why.empty())break;
  op->right_multiply(A,tmp);std::swap(A,tmp);if(q.need_stabilization[l]){A.adjointInPlace();if(seg==q.checkpoints-1)q.udtL[seg]=UDT(A);else q.udtL[seg]=A*q.udtL[seg+1];A=MatType::Identity(q.nDim,q.nDim);if(seg==0){q.udtL[seg].onePlusInv(q.g);q.g.adjointInPlace();}else q.g=onePlusInv(q.udtL[seg],q.udtR[seg-1]);if(check(q,o,"left_stabilization",0,l,l)){why="STABILIZATION";break;}--seg;}}
 std::cout<<"{\"seed\":"<<seed<<",\"event\":\""<<why<<"\"}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
