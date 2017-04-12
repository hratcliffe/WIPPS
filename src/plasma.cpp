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

/********Basic setup and allocation functions ****/
plasma::plasma(std::string file_prefix, my_type Bx_local){
/** \brief Set up plasma
*
*Sets up components from {file_prefix}plasma.conf. If a Bx_local is given, store and calc local cyclotron frequency from this. Else use the cyclotron frequency from deck constants.
*/

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
          pmass[block_num] = atof(val.c_str());
        }else if(name == "mass"){
          //Find which mass is relative to
          pos = val.find('*');
          tail = val.substr(pos+1, val.size());
          head = val.substr(0, pos);
          trim_string(tail);
          tail = str_to_lower(tail);
          trim_string(head);
          if(tail == "me") pmass[block_num] = atof(val.c_str())* me;
          else if(tail == "mp") pmass[block_num] = atof(val.c_str())* mp;
        
        }else if(name == "charge"){
          pcharge[block_num] = atof(val.c_str()) * q0;
        }else if(name == "dens"){
          pdens[block_num] = atof(val.c_str()) * ref_dens;
        }else if(name == "v_th"){
          pvth[block_num] = atof(val.c_str())*v0;
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
*Takes a two char code string and returns the specified frequency at local position. ce is actual Cyclotron freq. c0 is a reference value
*/

  code = str_to_lower(code);
#ifdef DEBUG_ALL
  std::string codes = "c0 pe ce";
  //Argument preconditions. Check only in debug mode for speed
  if(codes.find(code) == std::string::npos) my_error_print("!!!!!!!!Error in get_omega_ref, unknown code (code="+code+")!!!!!!", 0);
#endif
  if(code == "c0") return this->om_ce_ref;
  if(code == "pe") return my_const.omega_pe;
  if(code == "ce") return this->om_ce_local;
  else return 0.0;

}

void plasma::set_B0(my_type B0){
/** \brief Set B0
*
*Sets the local reference B field value and thus om_ce_local value.
*/
  this->B0 = B0;
  this->om_ce_local = std::abs(q0) * B0 / me;
}

/********Dispersion solvers ****/
mu_dmudom plasma::get_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type gamma_particle, bool skip_phi, bool Righthand)const{
/** \brief Solve plasma dispersion and extensions
*
*Solves Appleton-Hartree plasma dispersion and returns struct containing mu, its derivatives and error code. Also returns the Phi defined by Lyons 1974. I.e. the set of values needed to calculate D See \ref str
*Duplicated from mufunctions by CEJ Watt
*
*On notation: within this routine and plasma::get_root we use notation as from mufunctions3.f90. In the return values as defined in support.h we match with Lyons and Albert. Thus in my_mu, we have lat, r, theta, omega for polar coordinate, r, wave normal angle and wave frequency
*Also needs particle pitch angle alpha WATCH for Clares version which uses a different alpha entirely...
 */
 
#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(psi < 0 || psi >= pi) my_error_print("!!!!!!!!Error in get_phi_mu_om, pitch angle (psi="+mk_str(psi)+") out of range!!!!!!", 0);
  if(alpha < 0 || alpha >= pi) my_error_print("!!!!!!!!Error in get_phi_mu_om, particle pitch angle (alpha="+mk_str(alpha)+") out of range!!!!!!", 0);
  if(gamma_particle < 1) my_error_print("!!!!!!!!Error in get_phi_mu_om, particle gamma (gamma_particle="+mk_str(gamma_particle)+") out of range!!!!!!", 0);
  //I don't think there's an upper or lower bound on w we need to enforce, only positivity. We'll take abs below and succeed, but in debug mode also warn
  if(w < 0) my_error_print("!!!!!!Error in get_phi_mu_om, wave frequency (w="+mk_str(w)+") is negative!!!!!!", 0);
  //Nor any actual bounds on n
#endif

  w = std::abs(w);
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
*Duplicates plasma::get_phi_mu_omega() but using reduced form of Stix parameters corresponding to a high-density assumption assuming the first species is the electrons. This is mainly for comparison with the exact solution to validate this assumption.
*Also needs particle pitch angle alpha */

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(psi < 0 || psi >= pi) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, pitch angle (psi="+mk_str(psi)+") out of range!!!!!!", 0);
  if(alpha < 0 || alpha >= pi) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, particle pitch angle (alpha="+mk_str(alpha)+") out of range!!!!!!", 0);
  if(gamma_particle < 1) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, particle gamma (gamma_particle="+mk_str(gamma_particle)+") out of range!!!!!!", 0);
  //I don't think there's an upper or lower bound on w we need to enforce. Only positivity. We'll take abs below and succeed, but in debug mode also warn
  if(w < 0) my_error_print("!!!!!!Error in get_high_dens_phi_mu_om, wave frequency (w="+mk_str(w)+") is negative!!!!!!", 0);
  //Nor any bounds on n
#endif

  w = std::abs(w);
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
/** Goes from the Stix RLP to the final mu. This part is the same for high dens and normal, only the params differ*/

  calc_type w2, w3;
  w = std::abs(w);
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

  if( (mua2 > 0.0) || (mub2 > 0.0) ){
    if(Righthand){//Select Mode
      if(D < 0.0 ){ smu = 1.0; mu2 = mua2;} //see Albert [2005]
      else{smu = -1.0; mu2 = mub2;}
    }else{
      if(D < 0.0 ){ smu = -1.0; mu2 = mub2;}
      else{smu = 1.0; mu2 = mua2;}
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

    if(!skip_phi){
      calc_type term1, term2, term3, denom, tmp_bes, tmp_besp, tmp_besm, bessel_arg, D_mu2S, omega_n;

      D_mu2S = D / (mu2 - S);
      calc_type calc_n = (calc_type) n;
      omega_n = -1.0 * calc_n * my_const.omega_ce/gamma_particle;
      calc_type omega_n_slash_n = -1.0 * my_const.omega_ce/gamma_particle;
      //temporaries for simplicity

      term1 = (mu2* s2psi - P)/(mu2);
      
      denom = pow(D_mu2S*term1, 2) + c2psi*pow((P/mu2), 2);
      
      bessel_arg = tan(psi)*tan(alpha) * (w - omega_n)/omega_n_slash_n;// n x tan alpha (om - om_n)/om_n, but allowing n/n = 1 even when n = 0

      tmp_besp = boost::math::cyl_bessel_j(abs(n)+1, bessel_arg);
      tmp_besm = boost::math::cyl_bessel_j(abs(n)-1, bessel_arg);

      term2 = (1.0 + D_mu2S)*tmp_besp;
      term2 += (1.0 - D_mu2S)*tmp_besm;
      
      //tmp_bes = boost::math::cyl_bessel_j(abs(n), bessel_arg);
      if(n != 0) tmp_bes = 0.5*bessel_arg *(tmp_besp + tmp_besm)/calc_n;
      else tmp_bes = boost::math::cyl_bessel_j(abs(n), bessel_arg);
      //Use bessel identity to save time.
      // J_(n-1) + J_(n+1) = (2 n / arg) J_n
      term3 = scpsi*tmp_bes/std::tan(alpha);

      my_mu.phi = std::pow((0.5*term1*term2 + term3), 2)/denom;
      if(alpha == 0) my_mu.phi = 0;
    }
  }
  return my_mu;
}

std::vector<calc_type> plasma::get_resonant_omega(calc_type x, calc_type v_par, calc_type gamma_particle, int n)const{
/** \brief Solve plasma dispersion and doppler resonance simultaneously
*
*Obtains solutions of the Doppler resonance condition omega - k_par v_par = -n Omega_ce and a high-density approximation to the Whistler mode dispersion relation simultaneously. Assumes pure electron-proton plasma and uses cubic_solve. ONLY solutions between -om_ce_local and om_ce_local, excluding omega = 0, are considered. "Zero" solutions are those less than the GEN_PRECISION constant in support.h. If solutions are found, they're returned in vector, otherwise empty vector is returned. 
*
*Note that since k_parallel and v_parallel in resonant condition are signed, we will get multiple entries of ± omega for the corresponding ±k and ±n. These should be handled by the calling code, as k may or may not be handled with both signs
@param x Wave normal angle @param v_par Particle velocity to solve with @param n Resonance number.
*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(std::abs(v_par) >= v0) my_error_print("!!!!!!!!Error in get_resonant_omega, velocity (v_par="+mk_str(v_par)+") out of range!!!!!!", 0);
  if(gamma_particle < 1.0) my_error_print("!!!!!!!!Error in get_resonant_omega, particle gamma (gamma_particle="+mk_str(gamma_particle)+") out of range!!!!!!", 0);
#endif

  std::vector<calc_type> ret_vec;
  calc_type om_ce = this->get_omega_ref("ce");
  calc_type om_pe_loc = this->get_omega_ref("pe");
  calc_type om_ce_ref = this->get_omega_ref("c0");

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
  calc_type cos_th = std::cos(std::atan(x));
  calc_type vel = v_par / v0;
  //For clarity
  calc_type vel_cos = std::pow(vel * cos_th, 2);//Product of the two, for simplicity in expressions. But note is not meaningful, theta is the wave angle

  calc_type gamma_sq = gamma_particle * gamma_particle;

  //Calculate coefficients
  //To maintain best precision we solve for x = omega/omega_ce_ref so x ~ 1
  a = (vel_cos - 1.0) * gamma_sq;
  b = (vel_cos*cos_th*gamma_sq + 2.0*gamma_particle*n - gamma_sq*cos_th)*om_ce/om_ce_ref;
  c = ((2.0 * gamma_particle * n * cos_th - n * n)* std::pow(om_ce/om_ce_ref, 2) - std::pow(om_pe_loc/om_ce_ref, 2)*vel_cos*gamma_sq);
  d = -n*n*std::pow(om_ce/om_ce_ref, 3)*cos_th;
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

calc_type plasma::get_dispersion(my_type in, int wave_type, bool reverse, bool deriv, my_type theta)const{
/** \brief Solve analytic dispersion (approx)
*
* By default returns omega for a given k (see reverse and deriv params param). Uses local reference cyclotron and plasma frequencies and works with UNNORMALISED quantitites. NB: parameters out of range will silently return 0. NB: For Whistler modes this is an approximation and intended to be perfectly reversible. @param k Wavenumber @param wave_type wave species (see support.h) @param reverse Return k for input omega @param deriv Whether to instead return anayltic v_g @param theta Wavenormal angle, default 0.0  \todo Complete Xmode?*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(wave_type < WAVE_WHISTLER || wave_type > WAVE_X_LOW) my_error_print("!!!!!!!!Error in get_dispersion, wave_type unknown (wave_type="+mk_str(wave_type)+")!!!!!!", 0);
  //Theta is not conditioned because we reduce it below. This may change though
  //In is either a frequency or a wavenumber and thus has no limits, although may wish to add sanity limits
  if((wave_type == WAVE_X_LOW || wave_type == WAVE_X_UP) && deriv) my_error_print("!!!!This path not yet implemented", 0);
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

