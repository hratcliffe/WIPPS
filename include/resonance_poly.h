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

class resonance_poly {
protected:
  double om_ce;
  double om_pe;
  double om_ce_ref;
  std::vector<double> coeff;
  bool first_calc;
public:
  explicit resonance_poly();
  explicit resonance_poly(double om_ce, double om_pe, double om_ce_ref);
  void calculate_coeffs_no_ion(double v_par, double psi, int n, double gamma);
  void calculate_coeffs_full(double v_par, double psi, int n, double gamma);

};

//Boost's Newton-Raphson needs a pair of value and deriv. on calling
typedef std::pair<double, double> NRVals;
class NR_poly: public resonance_poly{
  public:
  NR_poly(double om_ce, double om_pe, double om_ce_ref):resonance_poly(om_ce, om_pe, om_ce_ref){;}
  NRVals operator()(const double x) { return std::make_pair(
    (((((((((coeff[10]*x + coeff[9])*x + coeff[8])*x + coeff[7])*x + coeff[6])*x + coeff[5])*x + coeff[4])*x + coeff[3])*x + coeff[2])*x + coeff[1])*x + coeff[0],
    ((((((((10.0*coeff[10]*x + 9.0*coeff[9])*x + 8.0*coeff[8])*x + 7.0*coeff[7])*x + 6.0*coeff[6])*x + 5.0*coeff[5])*x + 4.0*coeff[4])*x + 3.0*coeff[3])*x + 2.0*coeff[2])*x + coeff[1]);}
};


#endif /* defined(____resonance_poly__) */
