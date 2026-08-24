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

struct Walker:public Spinless_tV{int center=-1,nt=0,np=0;Walker(const SpinlessTvChainUtils*c,rdGenerator*r){build_projector_static_contour(*this,c,r,10,8,center,nt,np);}};
static double rel(const MatType&a,const MatType&b){return(a-b).norm()/std::max(b.norm(),std::numeric_limits<double>::min());}
static int sg(DataType z){return z.real()>=0?1:-1;}
int main(int ac,char**av)try{
 if(ac!=4)throw std::runtime_error("usage seed csv json");int seed=std::stoi(av[1]);std::ofstream out(av[2]);
 out<<"phase_index,stage,direction,iteration,boundary,tracked_re,tracked_im,direct_re,direct_im,sign_mismatch,green_rel,event\n"<<std::setprecision(17);
 SpinlessTvChainUtils c(6,.1,4,200,0,1,0,0);rdGenerator r(seed);Walker w(&c,&r);PfQMC q(&w,10);q.configureAdaptiveGuard(false,.1,100.);
 long long phase=0;std::string event;int hit_iter=-1;std::string hit_dir;
 auto probe=[&](const char*stage,const char*dir,int iter){MatType full;q.rebuildGreenFromFullContourAtBoundary(0,full);DataType direct=q.getSignRaw();double e=rel(q.g,full);bool sm=sg(q.sign)!=sg(direct);std::string ev=sm?"SIGN_MISMATCH":(e>1e-6?"GREEN_DRIFT":"");out<<phase++<<','<<stage<<','<<dir<<','<<iter<<",0,"<<q.sign.real()<<','<<q.sign.imag()<<','<<direct.real()<<','<<direct.imag()<<','<<sm<<','<<e<<','<<ev<<'\n';if(!ev.empty()){event=ev;hit_iter=iter;hit_dir=dir;return true;}return false;};
 if(probe("initial","none",-1))goto done;
 for(int i=0;i<500;++i){q.rightSweep();if(probe("burn","right",i))goto done;q.leftSweep();if(probe("burn","left",i))goto done;}
 for(int i=0;i<5000;++i){MatType center;DataType center_sign;q.rightSweep(w.center,&center,&center_sign);if(probe("measurement","right",i))goto done;q.leftSweep();if(probe("measurement","left",i))goto done;}
 done: std::ofstream js(av[3]);js<<"{\"seed\":"<<seed<<",\"event\":\""<<event<<"\",\"direction\":\""<<hit_dir<<"\",\"iteration\":"<<hit_iter<<",\"phase_index\":"<<(phase-1)<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
