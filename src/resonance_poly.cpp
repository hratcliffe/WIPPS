//
//  resonance_poly.cpp
//  
//
//  Created by Heather Ratcliffe on 01/06/2017.
//
//

#include "resonance_poly.h"

resonance_poly::resonance_poly(){

  this->om_ce = 1;
  this->om_pe = 1;
  this->om_ce_ref = 1;
  first_calc = true;
}

resonance_poly::resonance_poly(double om_ce, double om_pe, double om_ce_ref){

  this->om_ce = om_ce;
  this->om_pe = om_pe;
  this->om_ce_ref = om_ce_ref;
  coeff.resize(11);
  first_calc = true;

}


void resonance_poly::calculate_coeffs_no_ion(double psi, double v_par, int n, double gamma){
/** \brief Calculate coeffs for negligible ion frequencies
*
*Calculate the 10th order coefficients in the limit of omega, omega_ce, Omega_e >> omega_ci, Omega_i. These are normalised, so the polynomial returns solutions for omega/omega_ce_ref
@param psi Wave normal angle
@param v_par Parallel particle velocity
@param n Resonant number
@param gamma Particle gamma factor
\todo Pre-calc and stash as much as possible, in particular consider psi in A, B, C
*/

  double A6, A4, A2, A0;
  double mu44, mu34, mu24, mu14, mu04;
  double mu22, mu12, mu02;
  double B6, B4, B2, B0;
  double C8, C6, C4, C2;
  double L_sq, cos_psi2, sin_psi2, gamma_sq, L;
  double om_pe_tilde = om_pe/om_ce_ref, om_ce_tilde = om_ce/om_ce_ref;

///Coefficents in limit of electron freqs >> ion
  gamma_sq = gamma*gamma;
  L = std::pow(v0/ (v_par* std::cos(psi) * gamma), 2);
  L_sq = std::pow(L, 2);

  cos_psi2 = std::cos(psi)*std::cos(psi);
  sin_psi2 = std::sin(psi)*std::sin(psi);

  A6 = 1.0;
  A4 = - om_pe_tilde*om_pe_tilde - om_ce_tilde*om_ce_tilde;
  A2 = std::pow(om_pe_tilde*om_ce_tilde, 2)*cos_psi2;
  A0 = 0.0;

  B6 = 2.0;
  B4 = -2.0*(om_ce_tilde*om_ce_tilde + 2.0*om_pe_tilde*om_pe_tilde);
  B2 = 2.0*std::pow(om_pe_tilde, 4) + std::pow(om_ce_tilde*om_pe_tilde, 2)*(1.0 + cos_psi2);
  B0 = 0.0;
 
  C8 = 1.0;
  C6 = -om_ce_tilde*om_ce_tilde - 3.0* om_pe_tilde*om_pe_tilde;
  C4 = 3.0*std::pow(om_pe_tilde, 4) + std::pow(om_pe_tilde*om_ce_tilde, 2);
  C2 = -std::pow(om_pe_tilde, 6);

  mu44 = L_sq * gamma_sq*gamma_sq;
  mu34 = -4.0 *L_sq*n*gamma*gamma_sq *om_ce_tilde;
  mu24 = 6.0*L_sq*n*n*gamma_sq *om_ce_tilde*om_ce_tilde;
  mu14 = -4.0 * L_sq*n*n*n*gamma*std::pow(om_ce_tilde, 3);
  mu04 = L_sq*std::pow(n*om_ce_tilde, 4);
  
  mu22 = L*gamma_sq;
  mu12 = -2.0*L*n*gamma*om_ce_tilde;
  mu02 = L*n*n*om_ce_tilde*om_ce_tilde;
  
  //Set coefficients in vector
 coeff[10] = A6 *mu44 - B6 *mu22 + C8;
  coeff[9] = A6 *mu34 - B6 *mu12;
  coeff[8] = A6 *mu24 + A4 *mu44 - B6 *mu02 - B4 *mu22 + C6;
  coeff[7] = A6 *mu14 + A4 *mu34 - B4 *mu12;
  coeff[6] = A6 *mu04 + A4 *mu24 + A2 *mu44 - B4 *mu02 - B2 *mu22 + C4;
  coeff[5] = A4 *mu14 + A2 *mu34 - B2 *mu12;
  coeff[4] = A4 *mu04 + A2 *mu24 + A0 *mu44 - B2 *mu02 - B0 *mu22 + C2;
  coeff[3] = A2 *mu14 + A0 *mu34 - B0 *mu12;
  coeff[2] = A2 *mu04 + A0 *mu24 - B0 *mu02;
  coeff[1] = A0 *mu14;
  coeff[0] = A0 *mu04;
 
}


