//
//  plasma.cpp
//  
//
//  Created by Heather Ratcliffe on 07/10/2015.

#include <stdio.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <boost/math/special_functions/sign.hpp>
#include "plasma.h"
#include "resonance_poly.h"

/********Basic setup and allocation functions ****/
plasma::plasma(std::string file_prefix, my_type Bx_local){
/** \brief Set up plasma
*
*Sets up components from {file_prefix}plasma.conf. If a Bx_local is given, store and calc local cyclotron frequency from this. Else use the cyclotron frequency from deck constants. IMPORTANT: Make sure my_consts is defined (read_deck and share_consts) before creating plasma!
@param file_prefix File prefix prepended to all files read
@param Bx_local Local x-component of magnetic field
*/

  if(my_const.omega_ce == 0){
    my_error_print("!!! Reference om_ce is zero!", mpi_info.rank);
    is_setup = p_state::p_none;
    return;
  }
  //Set up plasma components
  plasma_state state = configure_from_file(file_prefix);
  
  //Set B0 and om_ce values
  if(Bx_local == -1){
    //Use reference B0
    B0 = my_const.omega_ce * me/std::abs(q0);
  }else{
    B0 = Bx_local;
  }
  this->om_ce_ref = my_const.omega_ce;
  this->om_ce_local = std::abs(q0) * B0 / me;

  is_setup = state;
  this->full_poly = new NR_poly(om_ce_local, my_const.omega_pe, om_ce_ref);
}

plasma_state plasma::configure_from_file(std::string file_prefix){
/** \brief Setup plasma from file
*
*Reads {file_prefix}plasma.conf and parses component mass, charge and density
*Sample block looks like
* electron:
*mass = 1.0*me
*charge = -1.0
*dens = 1.0
*v_th = 0.01
*Masses can be given relative to me or mp by ending the line with *mX. Charges are assumed relative to q0, densities to the reference_density given by omega_pe in deck constants, and v_th to c
On error we continue using defaults set below
@param file_prefix File prefix prepended to "plasma.conf" to get file to read
@return plasma_state object showing success and any caveats
*/

  calc_type ref_dens = my_const.omega_pe * my_const.omega_pe * eps0 * me / q0/q0;

//Default values -----------------------------------------
  pmass[0] = me;
  pmass[1] = mp;
  pmass[2] = mp;
  pmass[3] = mp;
  
  pcharge[0] = -1.0*q0;
  pcharge[1] = 1.0*q0;
  pcharge[2] = 0.0*q0;
  pcharge[3] = 0.0*q0;
  

  pdens[0] = 1.0 * ref_dens;
  pdens[1] = 1.0 * ref_dens;
  pdens[2] = 0.0 * ref_dens;
  pdens[3] = 0.0*ref_dens;

  pvth[0] = 0.01*v0;
  pvth[1] = 0.001*v0;
  pvth[2] = 0.001*v0;
  pvth[3] = 0.001*v0;

//End default values -----------------------------------------

  std::ifstream infile;
  infile.open(file_prefix+"plasma.conf");
  std::string line, name, val, head, tail;
  int block_num = -1;
  bool parse_err;
  size_t pos;

  //Very naive parsing. We spin through until we find a ":" and read the next lines until we find another
  //If we don't find n_comps such blocks, we report and continue
  
  while(getline(infile, line)){
    
    if(line.find(':') != std::string::npos){
      block_num ++;
      if(block_num >= ncomps) break;
      continue;
      //Found next block, so skip this header line
    }
    if(block_num >= 0){
      //This line might be a valid input one
      parse_err = parse_name_val(line, name, val);
      if(!parse_err){
        name = str_to_lower(name);
        if(name == "mass" && val.find('*') == std::string::npos){
          pmass[block_num] = checked_strtof(val.c_str());
        }else if(name == "mass"){
          //Find which mass is relative to
          pos = val.find('*');
          if(pos != std::string::npos){
            tail = val.substr(pos+1, val.size());
            head = val.substr(0, pos);
            trim_string(tail);
            tail = str_to_lower(tail);
            trim_string(head);
            val[pos] = 0;
          }
          if(tail == "me") pmass[block_num] = checked_strtof(val.c_str())* me;
          else if(tail == "mp") pmass[block_num] = checked_strtof(val.c_str())* mp;
          else pmass[block_num] = checked_strtof(val.c_str());
          
        }else if(name == "charge"){
          pcharge[block_num] = checked_strtof(val.c_str()) * q0;
        }else if(name == "dens"){
          pdens[block_num] = checked_strtof(val.c_str()) * ref_dens;
        }else if(name == "v_th"){
          pvth[block_num] = checked_strtof(val.c_str())*v0;
        }
      }
    }
  }

  if(block_num == -1){
    return p_default;
  }
  else if(block_num >= ncomps){
    my_print("Too many blocks in plasma file, truncating!", mpi_info.rank);
    return p_overflow;
  }
  else if(block_num < ncomps-1){
    my_print(mk_str(block_num), mpi_info.rank);
  
    my_print("Insufficient blocks in config file, using defaults for others", mpi_info.rank);
    return p_underflow;
  }else{
    return p_good;
  }

}

/********Get/set functions ****/
calc_type plasma::get_omega_ref(std::string code)const{
/** \brief Reference plasma and cyclotron frequencies
*
*Get value of omega at local position
@param code two character code string. ce is actual Cyclotron freq. c0 is a reference value. pe is plasma frequency
@return Value of reference omega
*/

  code = str_to_lower(code);
#ifdef DEBUG_ALL
  std::string codes = "c0 pe ce";
  //Argument preconditions. Check only in debug mode for speed
  if(codes.find(code) == std::string::npos) my_error_print("!!!!!!!!Error in get_omega_ref, unknown code (code="+code+")!!!!!!", mpi_info.rank);
#endif
  if(code == "c0") return this->om_ce_ref;
  if(code == "pe") return my_const.omega_pe;
  if(code == "ce") return this->om_ce_local;
  else return 0.0;

}

void plasma::set_B0(my_type B0){
/** \brief Set B0
*
*Sets the local reference B field value and thus om_ce_local value
@param B0 Input B0 value
*/
  this->B0 = B0;
  this->om_ce_local = std::abs(q0) * B0 / me;
}

/********Dispersion solvers ****/
mu_dmudom plasma::get_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type gamma_particle, bool skip_phi, bool Righthand)const{
/** \brief Solve plasma dispersion and extensions
*
*Solves Appleton-Hartree plasma dispersion and returns struct containing mu, its derivatives and error code. Also returns the Phi defined by Lyons \cite Lyons1974B. I.e. the set of values needed to calculate D See \ref mu_dmudom
*
*Duplicated from mufunctions by CEJ Watt
@param w Wave frequency 
@param psi Wave normal angle
@param alpha particle pitch angle (for phi) 
@param n Resonance number 
@param gamma_particle Relativistic gamma for resonant particle 
@param skip_phi Omit phi calculation 
@param Righthand True for Righthand wave mode, false for left
@return mu_dmudom object containing mu info
*
*On notation: within this routine we use notation as from mufunctions3.f90. In the return values as defined in support.h we match with Lyons \cite Lyons1974B and Albert \cite Albert2005 . Thus in my_mu, we have lat, r, theta, omega for polar coordinate, r, wave normal angle and wave frequency
 */
 
#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(psi < 0 || psi >= pi) my_error_print("!!!!!!!!Error in get_phi_mu_om, pitch angle (psi="+mk_str(psi)+") out of range!!!!!!", mpi_info.rank);
  if(alpha < 0 || alpha >= pi) my_error_print("!!!!!!!!Error in get_phi_mu_om, particle pitch angle (alpha="+mk_str(alpha)+") out of range!!!!!!", mpi_info.rank);
  if(gamma_particle < 1) my_error_print("!!!!!!!!Error in get_phi_mu_om, particle gamma (gamma_particle="+mk_str(gamma_particle)+") out of range!!!!!!", mpi_info.rank);
  //I don't think there's an upper or lower bound on w we need to enforce
  //Nor any actual bounds on n
#endif

  calc_type w2 = w*w;
  
  calc_type R=1.0, L=1.0, P=1.0;
  calc_type wp[ncomps], wp2[ncomps], wc[ncomps];

  for(int i = 0; i < ncomps; ++i){
    wp[i] = sqrt(pdens[i] * pcharge[i]*pcharge[i]/(eps0 * pmass[i]));
    wp2[i] = wp[i]*wp[i];
    wc[i] =  (pcharge[i]) * this->B0 / pmass[i];
  }
  
  for(int i=0; i<ncomps; ++i){
    R = R - wp2[i]/(w*(w + wc[i]));
    L = L - wp2[i]/(w*(w - wc[i]));
    P = P - wp2[i]/w2;
  }

  return get_phi_mu_from_stix(w, psi, alpha, n, gamma_particle, R, L, P, skip_phi, Righthand);
}

mu_dmudom plasma::get_high_dens_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type gamma_particle, bool skip_phi, bool Righthand)const{
  /** \brief Solve plasma dispersion and extensions
*
*Duplicates plasma::get_phi_mu_om but using reduced form of Stix parameters corresponding to a high-density assumption assuming the first species is the electrons. This is mainly for comparison with the exact solution to validate this assumption.
@param w Wave frequency 
@param psi Wave normal angle 
@param alpha particle pitch angle (for phi)
@param n Resonance number 
@param gamma_particle Relativistic gamma for resonant particle
@param skip_phi Omit phi calculation 
@param Righthand True for Righthand wave mode, false for left
@return mu_dmudom object containing mu info
\caveat Unsurprisingly this routine uses a high density approximation to the dispersion, which assumes \f$ \omega_{pe} >> \omega_{ce} \f$
\caveat This routine assumes that electrons are the first species of the plasma, i.e. the first species in the plasma.conf file
*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(psi < 0 || psi >= pi) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, pitch angle (psi="+mk_str(psi)+") out of range!!!!!!", mpi_info.rank);
  if(alpha < 0 || alpha >= pi) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, particle pitch angle (alpha="+mk_str(alpha)+") out of range!!!!!!", mpi_info.rank);
  if(gamma_particle < 1) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, particle gamma (gamma_particle="+mk_str(gamma_particle)+") out of range!!!!!!", mpi_info.rank);
  //I don't think there's an upper or lower bound on w we need to enforce.
  //Nor any bounds on n
#endif

  calc_type w2 = w*w;
  
  calc_type R=1.0, L=1.0, P=1.0;
  calc_type wp[ncomps], wp2[ncomps], wc[ncomps];

  for(int i=0; i<ncomps; ++i){
    wp[i] = sqrt(pdens[i] * pcharge[i]*pcharge[i]/(eps0 * pmass[i]));
    wp2[i] = wp[i]*wp[i];
    wc[i] =  (pcharge[i]) * this->B0 / pmass[i];

  }
  
  for(int i=0; i<1; ++i){
    //consider electron component only...
    R = 0.0 - wp2[i]/(w*(w + wc[i]));
    L = 0.0 - wp2[i]/(w*(w - wc[i]));
    P = 0.0 - wp2[i]/w2;
  }
  
  return get_phi_mu_from_stix(w, psi, alpha, n, gamma_particle, R, L, P, skip_phi, Righthand);
}

mu_dmudom plasma::get_phi_mu_from_stix(calc_type w, calc_type psi, calc_type alpha, int n, calc_type gamma_particle, calc_type R, calc_type L, calc_type P, bool skip_phi, bool Righthand)const{
/** Goes from the Stix RLP to the final mu. This part is the same for high dens and normal, only the calculations of the Stix params differ
@param w Wave frequency 
@param psi Wave normal angle
@param alpha particle pitch angle (for phi) 
@param n Resonance number 
@param gamma_particle Relativistic gamma for resonant particle 
@param R Stix param
@param L Stix param
@param P Stix param
@param skip_phi Omit phi calculation
@param Righthand True for Righthand wave mode, false for left
@return mu_dmudom object containing mu info
*/

  calc_type w2, w3;
  w2 = w*w;
  w3 = w2*w;
  
  calc_type J, S, D, A, B, C, s2psi, c2psi, mua2, mub2, mu2, scpsi;
  calc_type F, G, smu;
  calc_type dHdF, dHdG, dmudw, dFdpsi,dGdpsi, dAdpsi,dBdpsi;
  calc_type wp[ncomps], wp2[ncomps], wc[ncomps], X[ncomps], Y[ncomps];
  
  calc_type dmudX[ncomps], dmudY[ncomps];
  calc_type dPdX[ncomps], dLdX[ncomps], dRdX[ncomps], dSdX[ncomps];
  calc_type dAdX[ncomps], dBdX[ncomps], dCdX[ncomps], dFdX[ncomps], dGdX[ncomps];
  calc_type dLdY[ncomps], dRdY[ncomps], dSdY[ncomps];
  calc_type dAdY[ncomps], dBdY[ncomps], dCdY[ncomps], dFdY[ncomps], dGdY[ncomps], dXdw[ncomps], dYdw[ncomps];
  
  //These will hold suitably calc'd plasma frequency, square and cyclotron freq. If we have to derive from position, we do...
  for(int i=0; i<ncomps; ++i){
    wp[i] = sqrt(pdens[i] * pcharge[i]*pcharge[i]/(eps0 * pmass[i]));
    wp2[i] = wp[i]*wp[i];
    wc[i] =  (pcharge[i]) * this->B0 / pmass[i];

    X[i] = wp2[i]/w2;
    Y[i] = wc[i]/w;
  }
  
//We loop over components and trust compiler to unroll for us :) most of these will trivially vectorise anyway.

  scpsi = std::sin(psi);
  s2psi = std::pow(scpsi, 2);
  c2psi = 1.0 - s2psi;
  scpsi *= std::sqrt(c2psi);
  //To make it smokin'

  S = 0.5*(R + L);
  D = 0.5*(R - L);
  A = S*s2psi + P*c2psi;
  B = R*L*s2psi + P*S*(1.0+c2psi);
  C = P*R*L;
  J = std::sqrt(B*B - 4.0*A*C);

  mua2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B + J);
  mub2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B - J);

  mu_dmudom my_mu;
  //placeholder values if we can't fill...
  my_mu.mu = 1.0;
  my_mu.dmudom = 0.0;
  my_mu.dmudtheta = 0.0;
  my_mu.err = 1;
  my_mu.phi = 0;
  my_mu.cone_ang = pi;
  smu = 0;
  mu2 = 0;

  if( (mua2 > 0.0) || (mub2 > 0.0) ){
    //see Albert [2005] or Stix [1972] for mode selection
    if(Righthand){//Select Mode
      if(w > 0){
        if(D < 0.0){ smu = 1.0; mu2 = mua2;}
        else{smu = -1.0; mu2 = mub2;}
      }else{
        if(D < 0.0){smu = -1.0; mu2 = mub2;}
        else{smu = 1.0; mu2 = mua2;}
      }
    }else{
      if(w > 0){
        if(D > 0.0){ smu = 1.0; mu2 = mua2;}
        else{smu = -1.0; mu2 = mub2;}
      }else{
        if(D > 0.0){smu = -1.0; mu2 = mub2;}
        else{smu = 1.0; mu2 = mua2;}
      }
    }

    my_mu.mu = std::sqrt(mu2);
    my_mu.err = 0;

    F = 2.0*(A - B + C);
    G = 2.0*A - B + smu*J;
    dHdF = -1.0/G;
    dHdG = F/(G*G);

    for(int i=0; i<ncomps; i++){
      dPdX[i] = -1.0;
      dLdX[i] = -1.0/(1.0 - Y[i]);

      dRdX[i] = -1.0/(1.0 + Y[i]);
      dSdX[i] = 0.5*(dRdX[i] + dLdX[i]);
      dAdX[i] = 0.5*(dRdX[i] + dLdX[i])*s2psi + dPdX[i]*c2psi;
      dBdX[i] = (L*dRdX[i] + R*dLdX[i])*s2psi + (P*dSdX[i] + S*dPdX[i])*(1.0 + c2psi);
      dCdX[i] = P*R*dLdX[i] + P*L*dRdX[i] + R*L*dPdX[i];
      dFdX[i] = 2.0*(dAdX[i] - dBdX[i] + dCdX[i]);
      dGdX[i] = 2.0*dAdX[i] - dBdX[i] + (smu/J)*(B*dBdX[i] - 2.0*(A*dCdX[i] + C*dAdX[i]));
      dmudX[i] = (0.5/my_mu.mu)*(dHdF*dFdX[i] + dHdG*dGdX[i]);

      dRdY[i] = X[i]/( pow(1.0 + Y[i], 2) );
      dLdY[i] = -X[i]/( pow(1.0 - Y[i], 2) );
      dSdY[i] = 0.5*(dRdY[i] + dLdY[i]);

      dAdY[i] = dSdY[i]*s2psi;
      dBdY[i] = (L*dRdY[i] + R*dLdY[i])*s2psi + P*dSdY[i]*(1.0 + c2psi);
      dCdY[i] = P*(L*dRdY[i] + R*dLdY[i]);

      dFdY[i] = 2.0*(dAdY[i] - dBdY[i] + dCdY[i]);
      dGdY[i] = 2.0*dAdY[i] - dBdY[i] + (smu/J)*(B*dBdY[i] - 2.0*(A*dCdY[i] + C*dAdY[i]));

      dmudY[i] = (0.5/my_mu.mu)*(dHdF*dFdY[i] + dHdG*dGdY[i]);


      dXdw[i] = -2.0*wp2[i]/w3;
      dYdw[i] = -wc[i]/w2;

    }

    dAdpsi = sin(2.0*psi)*(S - P);
    dBdpsi = sin(2.0*psi)*(R*L - P*S);
    dFdpsi = 2.0*(dAdpsi - dBdpsi);
    dGdpsi = 2.0*dAdpsi - dBdpsi + (smu/J)*(B*dBdpsi - 2.0*C*dAdpsi);

    my_mu.dmudtheta = (0.5/my_mu.mu)*(dHdF*dFdpsi + dHdG*dGdpsi);

    dmudw = 0.0;
    for(int i=0; i<ncomps; i++) dmudw += dmudX[i]*dXdw[i] + dmudY[i]*dYdw[i];
    
    my_mu.dmudom = dmudw;
    my_mu.cone_ang = std::atan(std::sqrt(-P/S));

    if(!skip_phi){
      calc_type term1, term2, term3, denom, tmp_bes, tmp_besp, tmp_besm, bessel_arg, D_mu2S, omega_n;

      D_mu2S = D / (mu2 - S);
      calc_type calc_n = (calc_type) n;
      omega_n = -1.0 * calc_n * std::abs(get_omega_ref("ce"))/gamma_particle;
      calc_type omega_n_slash_n = -1.0 * std::abs(get_omega_ref("ce"))/gamma_particle;
      //temporaries for simplicity

      term1 = (mu2* s2psi - P)/(mu2);
      
      denom = std::pow(D_mu2S*term1, 2) + c2psi*std::pow((P/mu2), 2);
      
      bessel_arg = std::tan(psi)*std::tan(alpha) * (w - omega_n)/omega_n_slash_n;// n x tan alpha (om - om_n)/om_n, but allowing n/n = 1 even when n = 0

      tmp_besp = boost::math::cyl_bessel_j(n+1, bessel_arg);
      tmp_besm = boost::math::cyl_bessel_j(n-1, bessel_arg);

      term2  = (1.0 + D_mu2S)*tmp_besp;
      term2 += (1.0 - D_mu2S)*tmp_besm;
      
      if(n != 0) tmp_bes = 0.5*bessel_arg *(tmp_besp + tmp_besm)/calc_n;
      else tmp_bes = boost::math::cyl_bessel_j(0, bessel_arg);
      //Use bessel identity to save time.
      // J_(n-1) + J_(n+1) = (2 n / arg) J_n, n!=0
      term3 = scpsi*tmp_bes/std::tan(alpha);

      my_mu.phi = std::pow((0.5*term1*term2 + term3), 2)/denom;
      if(alpha == 0) my_mu.phi = 0;
    }
  }
  return my_mu;
}

std::vector<calc_type> plasma::get_resonant_omega(calc_type theta, calc_type v_par, calc_type gamma_particle, int n)const{
/** \brief Solve plasma dispersion and doppler resonance simultaneously
*
*Obtains solutions of the Doppler resonance condition omega - k_par v_par = -n Omega_ce and a high-density approximation to the Whistler mode dispersion relation simultaneously. Assumes pure electron-proton plasma and uses cubic_solve. ONLY solutions between -om_ce_local and om_ce_local, excluding omega = 0, are considered. "Zero" solutions are those less than the GEN_PRECISION constant in support.h.
*
*Note that since k_parallel and v_parallel in resonant condition are signed, we will get multiple entries of ± omega for the corresponding ±k and ±n. These should be handled by the calling code, as k may or may not be handled with both signs
@param theta Wave normal angle 
@param v_par Particle velocity to solve with 
@param gamma_particle Relativistic gamma for resonant particle
@param n Resonance number
@return Vector of solutions for resonant omega, or empty vector if no solutions are found
\ext Extend this to use full solution rather than the high density approx
*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(std::abs(v_par) >= v0) my_error_print("!!!!!!!!Error in get_resonant_omega, velocity (v_par="+mk_str(v_par)+") out of range!!!!!!", mpi_info.rank);
  if(gamma_particle < 1.0) my_error_print("!!!!!!!!Error in get_resonant_omega, particle gamma (gamma_particle="+mk_str(gamma_particle)+") out of range!!!!!!", mpi_info.rank);
#endif

  std::vector<calc_type> ret_vec;
  calc_type om_ce = this->get_omega_ref("ce");
  calc_type om_pe_loc = this->get_omega_ref("pe");
  calc_type om_ce_ref = this->get_omega_ref("c0");
  calc_type signed_om_ce_slash_ref = - std::abs(om_ce/om_ce_ref);//Contains sign of electron charge
  //Special case when v=0 and we can save time
  if(std::abs(v_par) < tiny_calc_type){
    if( n == 1 || n == -1 ) ret_vec.push_back(om_ce*n);
    return ret_vec;
  }
  //We let n=0 fall through even though we can handle it quicker because it's simpler code. We do exclude omega "=" zero solutions below though.

  //Equation to solve is cubic of the form
  //a x^3 + b x^2 + c x + d = 0
  //or equivalently
  //x^3 + an x^2 + bn x + cn = 0

  calc_type a, b, c, d;
  calc_type an, bn, cn;
  calc_type cos_th = std::cos(theta);
  calc_type vel = v_par / v0;
  calc_type calc_n = (calc_type) n;
  //For clarity
  calc_type vel_cos = std::pow(vel * cos_th, 2);//Product of the two, for simplicity in expressions. But note is not meaningful, theta is the wave angle

  calc_type gamma_sq = gamma_particle * gamma_particle;

  //Calculate coefficients
  //To maintain best precision we solve for x = omega/omega_ce_ref so x ~ 1
  a = (vel_cos - 1.0) * gamma_sq;

  b = (vel_cos*cos_th*gamma_sq + 2.0*gamma_particle* calc_n - gamma_sq*cos_th)*signed_om_ce_slash_ref;

  c = (2.0 * gamma_particle * calc_n * cos_th - calc_n * calc_n)* std::pow(signed_om_ce_slash_ref, 2) - std::pow(om_pe_loc/om_ce_ref, 2)*vel_cos*gamma_sq;

  d = -calc_n*calc_n*std::pow(signed_om_ce_slash_ref, 3)*cos_th;
  //Note: for n=0, d is 0. We covered v_par being zero above, so either omega = 0 and k_par = 0 or we have a normal solution, but with redundancy. We assume omega = 0 is an unhelpful solution to return. And for simplicity we don't do a special quadratic solution but just continue anyway

  an = b/a;
  bn = c/a;
  cn = d/a;

  ret_vec = cubic_solve(an, bn, cn);

  //Restore om_ce_ref factor and delete any entries > om_ce as these aren't whistler modes. Also delete answers which are "zero"
  for(size_t i=0; i<ret_vec.size(); ++i) ret_vec[i] *= om_ce_ref;
  for(size_t i=0; i<ret_vec.size(); ++i){
    if(std::abs(ret_vec[i]) > std::abs(om_ce) || std::abs(ret_vec[i]) < GEN_PRECISION){
      ret_vec.erase(ret_vec.begin() + i);
      --i;
    }
  }

  return ret_vec;

}

std::vector<calc_type> plasma::get_resonant_omega_full(calc_type theta, calc_type v_par, calc_type gamma_particle, int n)const{
/** \brief Solve plasma dispersion and doppler resonance simultaneously
*
*Obtains solutions of the Doppler resonance condition omega - k_par v_par = -n Omega_ce and the Whistler mode dispersion relation simultaneously. Assumes pure electron-proton plasma and uses cubic_solve. ONLY solutions between -om_ce_local and om_ce_local, excluding omega = 0, are considered. "Zero" solutions are those less than the GEN_PRECISION constant in support.h.
*
*Note that since k_parallel and v_parallel in resonant condition are signed, we will get multiple entries of ± omega for the corresponding ±k and ±n. These should be handled by the calling code, as k may or may not be handled with both signs
@param theta Wave normal angle 
@param v_par Particle velocity to solve with 
@param gamma_particle Relativistic gamma for resonant particle
@param n Resonance number
@return Vector of solutions for resonant omega, or empty vector if no solutions are found
\caveat I am assuming the solutions move but no more appear in the full equation. This may not be accurate, but we need a first guess for the solver
\caveat We set a hard minimum for resonant frequencies of interest, NR_min_om
*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(std::abs(v_par) >= v0) my_error_print("!!!!!!!!Error in get_resonant_omega, velocity (v_par="+mk_str(v_par)+") out of range!!!!!!", mpi_info.rank);
  if(gamma_particle < 1.0) my_error_print("!!!!!!!!Error in get_resonant_omega, particle gamma (gamma_particle="+mk_str(gamma_particle)+") out of range!!!!!!", mpi_info.rank);
#endif

  std::vector<calc_type> ret_vec, cubic_guesses;
  calc_type om_ce = this->get_omega_ref("ce");
  calc_type om_pe_loc = this->get_omega_ref("pe");
  calc_type om_ce_ref = this->get_omega_ref("c0");

  if(std::abs(v_par) < tiny_calc_type){
    if( n == 1 || n == -1 ) ret_vec.push_back(om_ce*n);
    return ret_vec;
  }
  //We let n=0 fall through even though we can handle it quicker because it's simpler code. We do exclude omega "=" zero solutions below though.

  //First attempt we just recalculate every time

  full_poly->calculate_coeffs_no_ion(theta, v_par, n, gamma_particle);
  cubic_guesses = this->get_resonant_omega(theta, v_par, gamma_particle, n);
  NR_poly callable_poly = *full_poly;
  double nrResult;
  for(size_t i = 0; i < cubic_guesses.size(); i++){
    nrResult = boost::math::tools::newton_raphson_iterate(callable_poly, cubic_guesses[i]*0.9/om_ce_ref, 0.6*cubic_guesses[i]/om_ce_ref, 1.0, 20);
    ret_vec.push_back(nrResult);
//    std::cout<<"NR returns "<<nrResult*om_ce_ref<<'\n';
//    std::cout<<"NR method found root matching to "<<(nrResult*om_ce_ref/cubic_guesses[i] - 1.0)*100.0<<" %"<<'\n';
  }

  //Restore om_ce_ref factor and delete any entries > om_ce as these aren't whistler modes. Also delete answers which are "zero"
  for(size_t i=0; i<ret_vec.size(); ++i) ret_vec[i] *= om_ce_ref;
  for(size_t i=0; i<ret_vec.size(); ++i){
    if(std::abs(ret_vec[i]) > std::abs(om_ce) || std::abs(ret_vec[i]) < GEN_PRECISION){
      ret_vec.erase(ret_vec.begin() + i);
      --i;
    }
  }

  return ret_vec;

}

calc_type plasma::get_dispersion(my_type in, int wave_type, bool reverse, bool deriv, my_type theta)const{
/** \brief Solve analytic dispersion (approx)
*
* By default returns omega for a given k (see reverse and deriv params param). Uses local reference cyclotron and plasma frequencies and works with UNNORMALISED quantitites. NB: parameters out of range will silently return 0. 
@param k Wavenumber 
@param wave_type wave species (see support.h) 
@param reverse Return k for input omega 
@param deriv Whether to instead return anayltic v_g 
@param theta Wavenormal angle, default 0.0  
@return Value of omega, or k if reverse is set
\todo Complete Xmode?
\caveat For Whistler modes this is an approximation and intended to be perfectly reversible.
*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(wave_type < WAVE_WHISTLER || wave_type > WAVE_X_LOW) my_error_print("!!!!!!!!Error in get_dispersion, wave_type unknown (wave_type="+mk_str(wave_type)+")!!!!!!", mpi_info.rank);
  //Theta is not conditioned because we reduce it below. This may change though
  //In is either a frequency or a wavenumber and thus has no limits, although may wish to add sanity limits
  if((wave_type == WAVE_X_LOW || wave_type == WAVE_X_UP) && deriv) my_error_print("!!!!This path not yet implemented", mpi_info.rank);
#endif


  calc_type ret = 0.0;
  calc_type om_ce_loc, om_pe_loc;
  calc_type cos_th;
  if(theta != 0.0){
  //Expect theta=0 to be common case so shortcut it
  //Else we step by step get theta into 0-pi/2 range
    if(theta < 0) theta = -theta;
    if(theta > 2.0*pi) theta = theta - 2.0*pi*(int)(theta/(2.0*pi));
    if(theta > pi) theta = 2.0*pi - theta;
    if(theta > pi/2.0) theta = pi/2.0 - theta;
  }
  cos_th = std::cos(theta);

  om_ce_loc = this->get_omega_ref("ce");
  om_pe_loc = this->get_omega_ref("pe");
  
  switch(wave_type){

    case WAVE_WHISTLER :
      if(!reverse){
        //This is omega(k) case
        calc_type csq_ksq = v0*v0*in*in;
        ret = v0*v0*std::abs(om_ce_loc)*cos_th/(csq_ksq + om_pe_loc*om_pe_loc);
        if(!deriv) ret *=std::pow(in, 2);
        else ret *= (2.0*in *std::pow(om_pe_loc, 2))/(csq_ksq + om_pe_loc*om_pe_loc);
      }else{
        //k(omega)
        //NB! We have omitted the one so that our k and omega approxes are inverses!!!
        in = std::abs(in);
        if(in >= om_ce_loc*cos_th || in == 0.0) break;
        
        if(!deriv){
          ret = in/v0*std::sqrt(0.0 - std::pow(om_pe_loc, 2)/(in*(in - std::abs(om_ce_loc)*cos_th)));
        }
        else{
          calc_type om_pe_sq = std::pow(om_pe_loc, 2);
          calc_type in_minus_om = in-std::abs(om_ce_loc);
          ret = 0.5*om_pe_sq*std::abs(om_ce_loc)/std::pow(in_minus_om, 2)/std::sqrt(-om_pe_sq*in/in_minus_om)  /v0;
        }
      }
      break;
    case WAVE_PLASMA:
      //Grab v_t as the first electron species in the plasma file
      {
      my_type v_t = 0.0;
      my_type om_ce_sin_sq = std::pow(om_ce_loc*std::sin(theta), 2);
      for(int i = 0; i < ncomps; i++){
        if(pmass[i] == me){
          v_t = pvth[i];
          break;
        }
      }
      if(!reverse){
      //This is "Z" mode, to get pure plasma use theta=0. Note also that it requires a single dominant electron thermal speed
        my_type kv_t_3_sq = 3.0*in*std::pow(v_t, 2);
        ret = std::sqrt(om_pe_loc*om_pe_loc + kv_t_3_sq*in + om_ce_sin_sq);
        if(deriv) ret = kv_t_3_sq/ret;
      }else{
        my_type tmp = in*in - om_pe_loc*om_pe_loc - om_ce_sin_sq;
        if(tmp > 0){
          if(deriv) ret = in/(sqrt(3.0*tmp) *v_t);
          else ret = std::sqrt(tmp/3.0)/v_t;
        }
      }
      break;
      }
    case WAVE_O :
      if(!reverse){
        ret = std::sqrt(om_pe_loc*om_pe_loc + std::pow(v0*in, 2));
        if(deriv) ret = v0*v0 * in/ret;
      }else{
        if(in < om_pe_loc) break;
        ret = std::sqrt(in*in -om_pe_loc*om_pe_loc);
        if(!deriv) ret = ret/v0;
        else ret = in/v0/ret;
      }
      break;
    case WAVE_X_UP:
    case WAVE_X_LOW:
      if(!reverse){
        int sgn = 1;
        if(wave_type == WAVE_X_LOW) sgn = -1;
        my_type c_sq_k_sq = std::pow(v0*in, 2), om_pe_sq = om_pe_loc*om_pe_loc, om_ce_sq = om_ce_loc*om_ce_loc, tmp;
        tmp = (sgn * std::sqrt(std::pow(om_ce_sq - c_sq_k_sq, 2) + 4.0 *om_pe_sq*om_ce_sq) + c_sq_k_sq + 2.0*om_pe_sq + om_ce_sq)/2.0;
        if(tmp != tmp ) break;//no solution...
        if(tmp > 0) ret = std::sqrt(tmp);
        else ret = std::sqrt(- tmp);
      }else{
        //How do we select correct branch??
        my_type om_pe_sq = om_pe_loc*om_pe_loc, om_sq = in * in;
        if(!deriv) ret = in/v0 * std::sqrt(1.0 - om_pe_sq/om_sq * (om_sq - om_pe_sq)/(om_sq - om_pe_sq - om_ce_loc*om_ce_loc));
      }
      break;
  }

  return ret;

}

