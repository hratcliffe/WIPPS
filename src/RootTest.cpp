#include <iostream>
#include <sstream>
#include <string>
 
 
#include <boost/math/tools/roots.hpp>
 
class Test {
public:
    double operator()(const double x) {
        //return x * cos(x);
        return x*x*x*x*x - 6.0*x*x*x*x + 20.0*x*x*x - 60.0*x*x + 99.0*x - 54.0;
    }
};
 
// see: 
// http://www.boost.org/doc/libs/1_47_0/libs/math/doc/sf_and_dist/html/math_toolkit/toolkit/internals1/roots2.html
 
int main(int argc, char** argv) {
    Test t;
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=500;
    boost::math::tools::eps_tolerance<double> tol(30);
 
    Result r1 = boost::math::tools::toms748_solve(t, 0.0, 2.5, tol, max_iter);
    std::cout << "root bracketed: [ " << r1.first << " , " << r1.second <<  " ]" << std::endl;
    std::cout << "f("<< r1.first << ")=" << t(r1.first) << std::endl;
    std::cout << "f("<< r1.second << ")=" << t(r1.second) << std::endl;
    std::cout << "max_iter=" << max_iter << std::endl;
    return 0;
}
