#include <iostream>
#include <sstream>
#include <string>
#include "support.h"
//#include "main_support.h"

#include <boost/math/tools/roots.hpp>
double * coeff;
double g_om_ce = 50000;
double g_om_pe = 10.0*g_om_ce;

#define calc_type double
std::vector<calc_type> cubic_result(calc_type x, calc_type v_par, calc_type n){
/**Get resonant frequency for particular x, v_parallel, n
*
*Solve high density approx to get omega. for pure electron proton plasma....
* Calls cubic_solve Note for slowly changing v_par, suggests Newtons method might be more efficient. Although much of this could be precomputed for given grids.
Return empty vector if no valid solutions
*/

  std::vector<calc_type> ret_vec;
  calc_type wc = g_om_ce;
  calc_type omega_pe_loc = g_om_pe;
  calc_type om_ref_ce = wc;

  calc_type a, b, c, d;
  calc_type an, bn, cn;
  calc_type cos_th = std::cos(std::atan(x));
  calc_type vel = v_par / v0;
  calc_type vel_cos = std::pow(vel * cos_th, 2) ;

  calc_type gamma, gamma2;
  gamma2 = 1.0/( 1.0 - std::pow(vel, 2));
  gamma = std::sqrt(gamma2);


  a = (vel_cos - 1.0) * gamma2;
  b = (vel_cos*cos_th*gamma2 + 2.0*gamma*n - gamma2*cos_th)*wc/om_ref_ce;
  c = ((2.0*gamma*n*cos_th - n*n)* std::pow(wc/om_ref_ce, 2) - std::pow(omega_pe_loc/om_ref_ce, 2)*vel_cos*gamma2);
  d = -n*n*std::pow(wc/om_ref_ce, 3)*cos_th;
  
  //To maintain best precision we solve for omega/ reference omega_ce
    
  an = b/a;
  bn = c/a;
  cn = d/a;
  
  ret_vec = cubic_solve(an, bn, cn);

  for(size_t i=0; i<ret_vec.size(); ++i) ret_vec[i] *= om_ref_ce;
  //restore factor
  for(size_t i=0; i<ret_vec.size(); ++i){
    if(std::abs(ret_vec[i]) > std::abs(wc)){
      ret_vec.erase(ret_vec.begin() + i);
      --i;
    }
  }
  //Discard any larger than omega_ce because they can't be whistlers
  return ret_vec;

}

class BuildPoly {
private:
  double om_pe_divide_ce_sq;
public:
  BuildPoly(double om_ce, double om_pe){ this->om_pe_divide_ce_sq = std::pow(om_pe/om_ce, 2);}
  void allocate(int n=11){ coeff = (double *) malloc(n*sizeof(double));}
  void calculate10(double v_par, double psi ,int n);
  void deallocate(){free(coeff);}
};

void BuildPoly::calculate10(double psi ,double v_par, int n){
//Calculate the coefficients

  double A6, A4, A2, A0;
  double mu44, mu34, mu24, mu14, mu04;
  double mu22, mu12, mu02;
  double B6, B4, B2, B0;
  double C8, C6, C4, C2;
  double L_sq, cos_psi2, sin_psi2, gamma, gamma_sq, L;

///Coefficents in limit of electron freqs >> ion

  cos_psi2 = std::cos(psi);
  gamma_sq = 1.0/(1.0 - std::pow(v_par/cos_psi2/v0, 2));
  gamma = std::sqrt(gamma_sq);
  L = std::pow(v0/ (v_par*cos_psi2 * gamma), 2);
  L_sq = std::pow(L, 2);

  cos_psi2 = cos_psi2 * cos_psi2;
  sin_psi2 = std::sin(psi);
  sin_psi2 = sin_psi2*sin_psi2;
  
  A6 = 1.0;
  A2 = - om_pe_divide_ce_sq * cos_psi2;
  A4 = A2 - 1.0;
  A0 = 0.0;
  
  B6 = 2.0;
  //B4 = -2.0 - om_pe_divide_ce_sq * (3.0 + sin_psi2);
  //B2 = 2.0 * om_pe_divide_ce_sq*om_pe_divide_ce_sq + om_pe_divide_ce_sq;
  B4 = -2.0 - om_pe_divide_ce_sq * 4.0;
  B2 = 2.0 * om_pe_divide_ce_sq*om_pe_divide_ce_sq;
  B0 = 0.0;
  
  C8 = 1.0;
  C6 = -1.0 - 3.0*om_pe_divide_ce_sq;
  C4 = 3.0*om_pe_divide_ce_sq*om_pe_divide_ce_sq + om_pe_divide_ce_sq;
  C2 = - om_pe_divide_ce_sq*om_pe_divide_ce_sq*om_pe_divide_ce_sq;
  
  mu44 = L_sq * gamma_sq*gamma_sq;
  mu34 = -4.0 *L_sq*n*gamma*gamma_sq;
  mu24 = 6.0*L_sq*n*n*gamma_sq;
  mu14 = -4.0 * L_sq*n*n*n*gamma;
  mu04 = L_sq*std::pow(n, 4);
  
  mu22 = L*gamma_sq;
  mu12 = -2.0*L*n*gamma;
  mu02 = L*n*n;
  
  
  coeff[0] = A6 *mu44 + B6 *mu22 + C8;
  coeff[1] = A6 *mu34 + B6 *mu12;
  coeff[2] = A6 *mu24 + A4 *mu44 + B6 *mu02 + B4 *mu22 + C6;
  coeff[3] = A6 *mu14 + A4 *mu34 + B4 *mu12;
  coeff[4] = A6 *mu04 + A4 *mu24 + A2 *mu44 + B4 *mu02 + B2 *mu22 + C4;
  coeff[5] = A4 *mu14 + A2 *mu34 + B2 *mu12;
  coeff[6] = A4 *mu04 + A2 *mu24 + A0 *mu44 + B2 *mu02 + B0 *mu22 + C2;
  coeff[7] = A2 *mu14 + A0 *mu34 + B0 *mu12;
  coeff[8] = A2 *mu04 + A0 *mu24 + B0 *mu02;
  coeff[9] = A0 *mu14;
 coeff[10] = A0 *mu04;
  
//   (((((((((coeff[0]*x + coeff[1])*x + coeff[2])*x + coeff[3])*x + coeff[4])*x + coeff[5])*x + coeff[6])*x + coeff[7])*x + coeff[8])*x + coeff[9])*x + coeff[10]

}

class Test {
public:
    double operator()(const double x) {
        //return x * cos(x);
        //return x*x*x*x*x - 6.0*x*x*x*x + 20.0*x*x*x - 60.0*x*x + 99.0*x - 54.0;
        return (((((((((coeff[0]*x + coeff[1])*x + coeff[2])*x + coeff[3])*x + coeff[4])*x + coeff[5])*x + coeff[6])*x + coeff[7])*x + coeff[8])*x + coeff[9])*x + coeff[10];
    }
};
 
// see: 
// http://www.boost.org/doc/libs/1_47_0/libs/math/doc/sf_and_dist/html/math_toolkit/toolkit/internals1/roots2.html
 
int main(int argc, char** argv) {
    Test t;
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=500;
    boost::math::tools::eps_tolerance<double> tol(30);
 
/*    Result r1 = boost::math::tools::toms748_solve(t, 0.0, 2.5, tol, max_iter);
    std::cout << "root bracketed: [ " << r1.first << " , " << r1.second <<  " ]" << std::endl;
    std::cout << "f("<< r1.first << ")=" << t(r1.first) << std::endl;
    std::cout << "f("<< r1.second << ")=" << t(r1.second) << std::endl;
    std::cout << "max_iter=" << max_iter << std::endl;
    return 0;*/
  
  calc_type x = 0.5, omega = g_om_ce/2.0, v_par = 0.1*v0, psi;
  int n=1;
  psi = std::atan(x);

  std::vector<calc_type> cubic = cubic_result(x, v_par, n);

  BuildPoly * tenth_order = new BuildPoly(g_om_ce, g_om_pe);
  tenth_order->allocate();
  tenth_order->calculate10(psi, v_par, n);

  for(int i=0; i< cubic.size(); ++i){
    std::cout<<cubic[i]/g_om_ce<<'\n';
    std::cout<<t(cubic[i]/g_om_ce)<<'\n';
  }
  double min, max;
  if(cubic.size()> 0){
    if(cubic[0] < 0){
      min = 1.1 * cubic[0]/g_om_ce;
      max= 0.9 * cubic[0]/g_om_ce;
    }else{
      min = 0.9 * cubic[0]/g_om_ce;
      max= 1.1 * cubic[0]/g_om_ce;
    
    }
  }
  
  min = -1.0;
  max = -0.01;
  
  Result r1 = boost::math::tools::toms748_solve(t, min, max, tol, max_iter);
    std::cout << "root bracketed: [ " << r1.first << " , " << r1.second <<  " ]" << std::endl;

  
}

std::vector<calc_type> cubic_solve(calc_type an, calc_type bn, calc_type cn){
/** \brief Finds roots of cubic x^3 + an x^2 + bn x + cn = 0
*
* Uses Num. Rec. equations, which are optimised for precision. Note that if x >>1 precision errors may result. Returns real solutions only
*/

  calc_type Q, R, bigA, bigB, Q3, R2, bigTheta;
  std::vector<calc_type> ret_vec;

  Q = (std::pow(an, 2) - 3.0 * bn)/9.0;
  R = (2.0* std::pow(an, 3) - 9.0 * an *bn + 27.0*cn)/54.0;
  
  R2 = std::pow(R, 2);
  Q3 = std::pow(Q, 3);
  
  if( R2 < Q3){
    
    bigTheta = std::acos(R/sqrt(Q3));
    calc_type minus2sqrtQ = -2.0*std::sqrt(Q);
    
    ret_vec.push_back(minus2sqrtQ*std::cos(bigTheta/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta + 2.0*pi)/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta - 2.0*pi)/3.0) - an/3.0);

  }else{
    calc_type ret_root;
    bigA = - boost::math::sign(R)*std::pow((std::abs(R) + std::sqrt(R2 - Q3)), 1.0/3.0 );

    (bigA != 0.0) ? (bigB = Q / bigA) : (bigB = 0.0);
    ret_root = (bigA + bigB) - an/3.0;

    ret_vec.push_back(ret_root);
  }

/** Used to test when writing
  calc_type tmp;
  for(int i=0; i<ret_vec.size(); ++i){
    
    tmp = std::pow(ret_vec[i], 3) + an*std::pow(ret_vec[i], 2) + bn*ret_vec[i] + cn;
    std::cout<<"solution gives "<<tmp<<std::endl;
  
  }
*/
  return ret_vec;

}

