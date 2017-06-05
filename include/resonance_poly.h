//
//  resonance_poly.h
//  
//
//  Created by Heather Ratcliffe on 01/06/2017.
//
//

#ifndef ____resonance_poly__
#define ____resonance_poly__

#include <stdio.h>
#include "support.h"

void check_stix_r(double psi, double n, double gamma, double k_par, double omega, double om_ce, double om_pe);//Check the given omega solves the equation for mu for given parameters

class resonance_poly {
protected:
  double om_ce;
  double om_pe;
  double om_ce_ref;
  double om_ci;
  double om_pi;
  std::vector<double> coeff;
  bool first_calc;
public:
  explicit resonance_poly();
  explicit resonance_poly(double om_ce, double om_pe, double om_ce_ref);
  explicit resonance_poly(double om_ce, double om_pe, double om_ce_ref, double om_ci, double om_pi);
  void calculate_coeffs_no_ion(double psi, double v_par, int n, double gamma);
  void calculate_coeffs_full(double psi, double v_par, int n, double gamma);

};

//Boost's Newton-Raphson needs a pair of value and deriv. on calling
typedef std::pair<double, double> NRVals;
class NR_poly: public resonance_poly{
  public:
  explicit NR_poly():resonance_poly(){;}
  explicit NR_poly(double om_ce, double om_pe, double om_ce_ref, double om_ci, double om_pi):resonance_poly(om_ce, om_pe, om_ce_ref, om_ci, om_pi){;}
  explicit NR_poly(double om_ce, double om_pe, double om_ce_ref):resonance_poly(om_ce, om_pe, om_ce_ref, om_ci, om_pi){;}

  NRVals operator()(const double x) { return std::make_pair(
    (((((((((coeff[10]*x + coeff[9])*x + coeff[8])*x + coeff[7])*x + coeff[6])*x + coeff[5])*x + coeff[4])*x + coeff[3])*x + coeff[2])*x + coeff[1])*x + coeff[0],
    ((((((((10.0*coeff[10]*x + 9.0*coeff[9])*x + 8.0*coeff[8])*x + 7.0*coeff[7])*x + 6.0*coeff[6])*x + 5.0*coeff[5])*x + 4.0*coeff[4])*x + 3.0*coeff[3])*x + 2.0*coeff[2])*x + coeff[1]);}

  bool sign_changes(double min, double max);
  bool is_root(double val, double tol = 1e-4);
  std::pair<double, double> refine_interval(std::pair<double, double> init, int n_steps);
  void dump_vals();
  
};


#endif /* defined(____resonance_poly__) */
