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
  coeff.resize(11);
  first_calc = true;
}

resonance_poly::resonance_poly(double om_ce, double om_pe, double om_ce_ref){

  this->om_ce = om_ce;
  this->om_pe = om_pe;
  this->om_ci = 0;
  this->om_pi = 0;
  this->om_ce_ref = om_ce_ref;
  coeff.resize(11);
  first_calc = true;

}

resonance_poly::resonance_poly(double om_ce, double om_pe, double om_ce_ref, double om_ci, double om_pi){

  this->om_ce = om_ce;
  this->om_pe = om_pe;
  this->om_ce_ref = om_ce_ref;
  this->om_ci = om_ci;
  this->om_pi = om_pi;

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

  mu44 = L_sq * gamma_sq * gamma_sq;
  mu34 = -4.0 * L_sq * n * std::pow(gamma, 3) * om_ce_tilde;
  mu24 = 6.0 * L_sq * n * n * gamma_sq * om_ce_tilde *om_ce_tilde;
  mu14 = -4.0 * L_sq *  gamma * std::pow(n*om_ce_tilde, 3);
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

void resonance_poly::calculate_coeffs_full(double psi, double v_par, int n, double gamma){
/** \brief Calculate coeffs for general case
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
  double om_pe_tilde = om_pe/om_ce_ref, om_ce_tilde = om_ce/om_ce_ref, om_pi_tilde = om_pi/om_ce_ref, om_ci_tilde = om_ci/om_ce_ref, om_p_tilde_sq = (om_pe_tilde*om_pe_tilde + om_pi_tilde*om_pi_tilde);

///Coefficents in limit of electron freqs >> ion
  gamma_sq = gamma*gamma;
  L = std::pow(v0/ (v_par* std::cos(psi) * gamma), 2);
  L_sq = std::pow(L, 2);

  cos_psi2 = std::cos(psi)*std::cos(psi);
  sin_psi2 = std::sin(psi)*std::sin(psi);

  A6 = 1.0;
  A4 = - om_p_tilde_sq - om_ce_tilde*om_ce_tilde - om_ci_tilde*om_ci_tilde;
  A2 = std::pow(om_ci_tilde*om_ce_tilde, 2) - sin_psi2*(std::pow(om_pe_tilde*om_ci_tilde, 2) + std::pow(om_pi_tilde*om_ce_tilde, 2)) + om_p_tilde_sq*(om_ce_tilde*om_ce_tilde + om_ci_tilde*om_ci_tilde)*cos_psi2;
  A0 = -om_p_tilde_sq *std::pow(om_ce_tilde*om_ci_tilde, 2)*cos_psi2;

  B6 = 2.0;
  B4 = -2.0*(om_ce_tilde*om_ce_tilde +om_ci_tilde*om_ci_tilde + 2.0*om_p_tilde_sq);
  B2 = 2.0*std::pow(om_ce_tilde*om_ci_tilde, 2) + (2.0+sin_psi2)*(std::pow(om_pe_tilde*om_ci_tilde, 2) + std::pow(om_pi_tilde*om_ce_tilde, 2)) + 2.0*(std::pow(om_pe_tilde, 4) + std::pow(om_pi_tilde, 4)) + 4.0*std::pow(om_pe_tilde*om_pi_tilde, 2) + om_p_tilde_sq*(om_ce_tilde*om_ce_tilde + om_ci_tilde*om_ci_tilde)*(1.0 + cos_psi2);
  B0 = -2.0*(std::pow(om_pe_tilde, 4)*om_ci_tilde*om_ci_tilde + std::pow(om_pi_tilde, 4)*om_ce_tilde*om_ce_tilde) - om_p_tilde_sq*std::pow(om_ce_tilde*om_ci_tilde, 2)*(1.0+cos_psi2) + std::pow(om_pe_tilde*om_pi_tilde, 2)*( (1.0+ cos_psi2)*std::pow(om_ce_tilde - om_ci_tilde, 2) - 4.0*om_ce_tilde*om_ci_tilde);
  ;
 
  C8 = 1.0;
  C6 = -om_ce_tilde*om_ce_tilde -om_ci_tilde*om_ci_tilde - 3.0* om_p_tilde_sq;
  C4 = 3.0*std::pow(om_p_tilde_sq, 2) + om_p_tilde_sq*std::pow(om_ce_tilde - om_ci_tilde, 2) + std::pow(om_ce_tilde*om_ci_tilde, 2);
  C2 = -om_p_tilde_sq*std::pow(om_ce_tilde*om_ci_tilde - om_p_tilde_sq, 2);

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

bool NR_poly::sign_changes(double min, double max){
/** \brief Check if sign changes on interval
*
*Checks for sign change between min and max. Does not assume min <= max
@param min Minimum of interval
@param max Maximum of interval
@return True if sign changes between min and max, false else
*/

  return (this->operator()(min).first/this->operator()(max).first < 0);
}

bool NR_poly::is_root(double val, double tol){
/** \brief Check if val is a root
*
*Checks if value is a root, by checking for a sign change in the interval val*(1-tol) to val*(1+tol)
@param min Value to check
@param tol (optional) Fractional width of interval to check, default 1e-4
@return True if val is (probably) a root, false else
*/
  return sign_changes(val*(1.0-tol), val*(1.0+tol));
}

void NR_poly::dump_vals(){

  std::ofstream outfile;
  outfile.open("./polyvals.dat", std::ios::out);
  double val = 0;
  for(int i = -200; i< 200; i++){
    val = (float) i /(float) 200;
    outfile<<val<<' '<<this->operator()(val).first<<'\n';
  }
  outfile.close();

}

std::pair<double, double> NR_poly::refine_interval(std::pair<double, double> init, int n_steps){

  double min, max, midpt;
  min = init.first;
  max = init.second;
  if(!sign_changes(min, max)) return std::make_pair(min, max);

  for(int i = 0; i < n_steps; i++){
    midpt = (min + max)/2.0;
    if(sign_changes(min, midpt)){
      max = midpt;
    }else if(sign_changes(midpt, max)){
      min = midpt;
    }else{
      min = midpt*0.99;
      max = midpt*1.01;
    }
  }

  return std::make_pair(min, max);

}


void check_stix_r(double psi, double n, double gamma, double k_par, double omega, double om_ce, double om_pe){

  double P, R, L, S;
  
  P = 1.0 - std::pow(om_pe/omega, 2);
  R = 1.0 - std::pow(om_pe, 2)/(omega*(omega + om_ce));
  L = 1.0 - std::pow(om_pe, 2)/(omega*(omega - om_ce));
  S = 0.5*(R + L);
  
  double sin_psi2 = std::pow(std::sin(psi), 2), cos_psi = std::cos(psi), cos_psi2 = cos_psi*cos_psi;
  double mu = v0 * k_par/cos_psi/omega;
  double A = S*sin_psi2 + P*cos_psi2;
  double B = R*L*sin_psi2 + P*S*(1.0 + cos_psi2);
  double C = P*R*L;
  std::cout<<"Stix gives "<<A*std::pow(mu, 4) - B* std::pow(mu, 2) + C <<'\n';
  
}

