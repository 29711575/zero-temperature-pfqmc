#include <iomanip>
#include <iostream>
#include "kitaevChain.h"
#include "pfqmc.h"
int main(){int L=4;SpinlessTvChainUtils c(L,.1,0,40,1,1,0,0);rdGenerator r(1);Chain_tV w(&c,&r);PfQMC q(&w,10);MatType g=q.g;double p=-c.StructureFactorCDW(g).real(),d=-c.StructureFactorCDWOffset(g).real();std::cout<<std::setprecision(17)<<p<<' '<<d<<' '<<1-d/p<<'\n';}
