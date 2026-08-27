#include <cmath>
#include <iomanip>
#include <iostream>
#include "pfqmc.h"

class PhaseOperator:public Operator{
    DataType phase;
public:
    explicit PhaseOperator(double angle):phase(std::polar(1.0,angle)){}
    DataType update(MatType&)override{return phase;}
    void left_multiply(const MatType&a,MatType&b)override{b=a;}
    void right_multiply(const MatType&a,MatType&b)override{b=a;}
    void left_propagate(MatType&,MatType&)override{}
    void right_propagate(MatType&,MatType&)override{}
    void stabilizedLeftMultiply(UDT&)override{}
    void getGreensMat(MatType&g)override{g=MatType::Zero(2,2);}
    DataType getSignOfWeight()override{return phase;}
};
class Walker:public Spinless_tV{public:Walker(){nDim=2;op_array.push_back(new PhaseOperator(.23));op_array.push_back(new PhaseOperator(.23));}};
int main(){Walker w;PfQMC q(&w,1);DataType before=q.sign;q.rightSweep();q.leftSweep();DataType expected=before*std::polar(1.0,.92);double error=std::abs(q.sign-expected);std::cout<<std::setprecision(17)<<"{\"mode\":\"generic_complex\",\"before_re\":"<<before.real()<<",\"before_im\":"<<before.imag()<<",\"after_re\":"<<q.sign.real()<<",\"after_im\":"<<q.sign.imag()<<",\"expected_error\":"<<error<<",\"real_z2_mode\":"<<(q.realZ2Mode()?"true":"false")<<"}\n";return error<1e-12&&!q.realZ2Mode()&&std::abs(q.sign.imag())>1e-3?0:2;}
