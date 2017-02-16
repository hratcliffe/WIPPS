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

extern deck_constants my_const;
extern const mpi_info_struc mpi_info;

/********Basic setup and allocation functions ****/
plasma::plasma(std::string file_prefix, my_type Bx_local){
/** \brief Set up plasma
*
*Sets up components from <file_prefix>plasma.conf. If a Bx_local is given, store and calc local cyclotron frequency from this. Else use the cyclotron frequency from deck constants. \todo What is is_setup actually doing?
*/

  //Set up plasma components
  configure_from_file(file_prefix);
  
  //Set B0 and om_ce values
  if(Bx_local == -1){
    //Use reference B0
    B0 = my_const.omega_ce * me/std::abs(q0);
  }else{
    B0 = Bx_local;
  }
  this->om_ce_ref = my_const.omega_ce;
  this->om_ce_local = std::abs(q0) * B0 / me;

  is_setup = true;
}

bool plasma::configure_from_file(std::string file_prefix){
/** \brief Setup plasma from file
*
*Reads <file_prefix>plasma.conf and parses component mass, charge and density
*Sample block looks like
* electron:
*mass = 1.0*me
*charge = -1.0
*dens = 1.0
*Masses can be given relative to me or mp by ending the line with *mX. Charges are assumed relative to q0 and densities to the reference_density given by omega_pe in deck constants
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
        if(name == "mass" && val.find('*') == std::string::npos){
          pmass[block_num] = atof(val.c_str());
        }else if(name == "mass"){
          //Find which mass is relative to
          pos = val.find('*');
          tail = val.substr(pos+1, val.size());
          head = val.substr(0, pos);
          trim_string(tail);
          trim_string(head);
          if(tail == "me") pmass[block_num] = atof(val.c_str())* me;
          else if(tail == "mp") pmass[block_num] = atof(val.c_str())* mp;
        
        }else if(name == "charge") pcharge[block_num] = atof(val.c_str()) * q0;
        else if(name == "dens") pdens[block_num] = atof(val.c_str()) * ref_dens;
      }
    }
  }

  if(block_num >= ncomps){
    my_print("Too many blocks in plasma file, truncating!", mpi_info.rank);
    return 1;
  }
  else if(block_num < ncomps-1){
    my_print(mk_str(block_num), mpi_info.rank);
  
    my_print("Insufficient blocks in config file, using defaults for others", mpi_info.rank);
    return 1;
  }else{
    return 0;
  }

}

/********Get/set functions ****/
calc_type plasma::get_omega_ref(std::string code)const{
/** \brief Reference plasma and cyclotron frequencies
*
*Takes a two char code string and returns the specified frequency at local position. ce is actual Cyclotron freq. c0 is a reference value
*/

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
mu plasma::get_root(calc_type th, calc_type w, calc_type psi, bool Righthand){
/** \brief Solve plasma dispersion 
*
*Solves Appleton-Hartree plasma dispersion and returns struct containing mu, its derivatives and error code. See \ref str
*Duplicated from mufunctions by CEJ Watt
*
*@param th Polar coordinate, or latitude @param w Wave frequency @param psi Wave pitch angle k to B0 \todo Check constant types OK \todo Better root selection
*
*On notation: within this routine and plasma::get_phi_mu_om() we use notation as from mufunctions3.f90. In the return values as defined in support.h we match with Lyons and Albert. Thus in my_mu, we have lat, r, theta, omega for polar coordinate, r, wave normal angle and wave frequency
*/
  mu mu_ret;
  calc_type dndr[ncomps], dndth[ncomps];
  calc_type dB0dr, dB0dth;

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(psi < 0 || psi >= pi) my_error_print("!!!!!!!!Error in get_root, pitch angle (psi="+mk_str(psi)+") out of range!!!!!!", 0);
  if(th < 0 || th >= pi/2.0) my_error_print("!!!!!!!!Error in get_root, position angle (th="+mk_str(th)+") out of range!!!!!!", 0);
  //I don't think there's an upper or lower bound on w we need to enforce
  
#endif

  
  /** \todo get or calc these...*/
  // FAKENUMBERS
  for(int i=0;i<ncomps; i++){
    dndr[i] = 0.0;
    dndth[i] = 0.0;
  }
  dB0dr = 1.0;
  dB0dth = 1.0;
  
  calc_type R=1.0, L=1.0, P=1.0, S, D, A, B, C, J, s2psi, c2psi, mua2, mub2, mu2;
  calc_type F, G, smu;
  calc_type dHdF, dHdG, dmudw,dAdpsi,dBdpsi,dFdpsi,dGdpsi,dpsidth;
  calc_type wp[ncomps], wp2[ncomps], wc[ncomps], X[ncomps], Y[ncomps];
  
  calc_type dmudX[ncomps], dXdr[ncomps], dXdth[ncomps], dmudY[ncomps], dYdr[ncomps], dYdth[ncomps];
  calc_type dPdX[ncomps], dLdX[ncomps], dRdX[ncomps], dSdX[ncomps], dDdX[ncomps];
  calc_type dAdX[ncomps], dBdX[ncomps], dCdX[ncomps], dFdX[ncomps], dGdX[ncomps];
  calc_type dLdY[ncomps], dRdY[ncomps], dSdY[ncomps], dDdY[ncomps];
  calc_type dAdY[ncomps], dBdY[ncomps], dCdY[ncomps], dFdY[ncomps], dGdY[ncomps], dXdw[ncomps], dYdw[ncomps];
  
  //These will hold suitably calc'd plasma frequency, square and cyclotron freq. If we have to derive from position, we do...
  for(int i=0; i<ncomps; ++i){
    wp[i] = sqrt(pdens[i] * pcharge[i]*pcharge[i]/(eps0 * pmass[i]));
    wp2[i] = wp[i]*wp[i];
    wc[i] =  (pcharge[i]) * this->B0 / pmass[i];

    X[i] = wp2[i]/(w*w);
    Y[i] = wc[i]/w;

  }
//We loop over components and trust compiler to unroll for us :) most of these will trivially vectorise anyway.

  s2psi = std::pow(sin(psi), 2);
  c2psi = 1.0 - s2psi;
  //To make it smokin'

  for(int i=0; i<ncomps; ++i){
    R = R - wp2[i]/(w*(w + wc[i]));
    L = L - wp2[i]/(w*(w - wc[i]));
    P = P - wp2[i]/(w*w);
  }


  S = 0.5*(R + L);
  D = 0.5*(R - L);
  A = S*s2psi + P*c2psi;
  B = R*L*s2psi + P*S*(1.0+c2psi);
  C = P*R*L;
  J = sqrt(B*B - 4.0*A*C);


  mua2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B + J);
  mub2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B - J);

  //placeholder values if we can't fill...
  mu_ret.mu = 1.0;
  mu_ret.mug = 0.0;
  mu_ret.dmudr = 0.0;
  mu_ret.dmudlat = 0.0;
  mu_ret.dmudom = 0.0;
  mu_ret.dmudtheta = 0.0;
  mu_ret.alpha = 0.0;
  mu_ret.err = 1;
  
  if( (mua2 > 0.0) || (mub2 > 0.0) ){
  
    if(Righthand){//Select Mode
      if(D < 0.0 ){ smu = 1.0; mu2 = mua2;} //see Albert [2005]
      else{smu = -1.0; mu2 = mub2;}
    }else{
      if(D < 0.0 ){ smu = -1.0; mu2 = mub2;}
      else{smu = 1.0; mu2 = mua2;}
    }
    mu_ret.mu = sqrt(mu2);
    
    F = 2.0*(A - B + C);
    G = 2.0*A - B + smu*J;
    dHdF = -1.0/G;
    dHdG = F/(G*G);

    for(int i=0; i<ncomps; i++){
      dPdX[i] = -1.0;
      dLdX[i] = -1.0/(1.0 - Y[i]);
      dRdX[i] = -1.0/(1.0 + Y[i]);
      dSdX[i] = 0.5*(dRdX[i] + dLdX[i]);
      dDdX[i] = 0.5*(dRdX[i] - dLdX[i]);
      dAdX[i] = 0.5*(dRdX[i] + dLdX[i])*s2psi + dPdX[i]*c2psi;
      dBdX[i] = (L*dRdX[i] + R*dLdX[i])*s2psi + (P*dSdX[i] + S*dPdX[i])*(1.0 + c2psi);
      dCdX[i] = P*R*dLdX[i] + P*L*dRdX[i] + R*L*dPdX[i];
      dFdX[i] = 2.0*(dAdX[i] - dBdX[i] + dCdX[i]);
      dGdX[i] = 2.0*dAdX[i] - dBdX[i] + (smu/J)*(B*dBdX[i] - 2.0*(A*dCdX[i] + C*dAdX[i]));
      dmudX[i] = (0.5/mu_ret.mu)*(dHdF*dFdX[i] + dHdG*dGdX[i]);

      dRdY[i] = X[i]/( pow(1.0 + Y[i], 2) );
      dLdY[i] = -X[i]/( pow(1.0 - Y[i], 2) );
      dSdY[i] = 0.5*(dRdY[i] + dLdY[i]);
      dDdY[i] = 0.5*(dRdY[i] - dLdY[i]);
      dAdY[i] = dSdY[i]*s2psi;
      dBdY[i] = (L*dRdY[i] + R*dLdY[i])*s2psi + P*dSdY[i]*(1.0 + c2psi);
      dCdY[i] = P*(L*dRdY[i] + R*dLdY[i]);
      dFdY[i] = 2.0*(dAdY[i] - dBdY[i] + dCdY[i]);
      dGdY[i] = 2.0*dAdY[i] - dBdY[i] + (smu/J)*(B*dBdY[i] - 2.0*(A*dCdY[i] + C*dAdY[i]));
      dmudY[i] = (0.5/mu_ret.mu)*(dHdF*dFdY[i] + dHdG*dGdY[i]);

      dXdr[i] = q0*q0*dndr[i]/(w*w*eps0*pmass[i]);
      dYdr[i] = pcharge[i]*dB0dr/(w*pmass[i]);
      dXdth[i] = q0*q0*dndth[i]/(w*w*eps0*pmass[i]);
      dYdth[i] = pcharge[i]*dB0dth/(w*pmass[i]);
      dXdw[i] = -2.0*wp2[i]/(std::pow(w, 3));
      dYdw[i] = -wc[i]/(std::pow(w, 2));

    }
    
    dAdpsi = sin(2.0*psi)*(S - P);
    dBdpsi = sin(2.0*psi)*(R*L - P*S);
    dFdpsi = 2.0*(dAdpsi - dBdpsi);
    dGdpsi = 2.0*dAdpsi - dBdpsi + (smu/J)*(B*dBdpsi - 2.0*C*dAdpsi);
    mu_ret.dmudtheta = (0.5/mu_ret.mu)*(dHdF*dFdpsi + dHdG*dGdpsi);

    dpsidth = -2.0/(1.0 + 3.0*std::pow(cos(th), 2));
    dmudw = 0.0;
    
    mu_ret.dmudlat = mu_ret.dmudtheta*dpsidth;
    //Even if this one can be folded into above, keep it out as not vectorisable
    for(int i=0; i<ncomps; i++){
       mu_ret.dmudr = mu_ret.dmudr + dmudX[i]*dXdr[i] + dmudY[i]*dYdr[i];
       mu_ret.dmudlat = mu_ret.dmudlat + dmudX[i]*dXdth[i] + dmudY[i]*dYdth[i];
       dmudw = dmudw + dmudX[i]*dXdw[i] + dmudY[i]*dYdw[i];
    }
    mu_ret.mug = mu_ret.mu + w*dmudw;
    mu_ret.alpha = -mu_ret.dmudtheta/mu_ret.mu;
    mu_ret.dmudom = dmudw;

    mu_ret.err = 0;
  
  }

  return mu_ret;
}

mu_dmudom plasma::get_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n, bool Righthand)const{
/** \brief Solve plasma dispersion and extensions
*
*Solves Appleton-Hartree plasma dispersion and returns struct containing mu, its derivatives and error code. Also returns the Phi defined by Lyons 1974. I.e. the set of values needed to calculate D See \ref str
*Duplicated from mufunctions by CEJ Watt
*
*On notation: within this routine and plasma::get_root we use notation as from mufunctions3.f90. In the return values as defined in support.h we match with Lyons and Albert. Thus in my_mu, we have lat, r, theta, omega for polar coordinate, r, wave normal angle and wave frequency
*Also needs particle pitch angle alpha \todo Fix relativistic gamma... WATCH for Clares version which uses a different alpha entirely... \todo Why do we pass omega_n??

 */
 
  calc_type term1, term2, term3, denom, tmp_bes, tmp_besp, tmp_besm, bessel_arg, D_mu2S, gamma, w2, w3;

  mu_dmudom my_mu;

  calc_type dndr[ncomps], dndth[ncomps];
  calc_type dB0dr, dB0dth;

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(psi < 0 || psi >= pi) my_error_print("!!!!!!!!Error in get_phi_mu_om, pitch angle (psi="+mk_str(psi)+") out of range!!!!!!", 0);
  if(alpha < 0 || alpha >= pi) my_error_print("!!!!!!!!Error in get_phi_mu_om, particle pitch angle (alpha="+mk_str(alpha)+") out of range!!!!!!", 0);
  //I don't think there's an upper or lower bound on w we need to enforce
  //Nor any actual bounds on n
#endif
  
  w2 = w*w;
  w3 = w2*w;
  
  /** \todo get or calc these...*/
  // FAKENUMBERS
  for(int i=0;i<ncomps; i++){
    dndr[i] = 0.0;
    dndth[i] = 0.0;
  }
  dB0dr = 0.0;
  dB0dth = 0.0;
  
  calc_type R=1.0, L=1.0, P=1.0, J, S, D, A, B, C, s2psi, c2psi, mua2, mub2, mu2, scpsi;
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

  for(int i=0; i<ncomps; ++i){
    R = R - wp2[i]/(w*(w + wc[i]));
    L = L - wp2[i]/(w*(w - wc[i]));
    P = P - wp2[i]/w2;
  }

  S = 0.5*(R + L);
  D = 0.5*(R - L);
  A = S*s2psi + P*c2psi;
  B = R*L*s2psi + P*S*(1.0+c2psi);
  C = P*R*L;
  J = std::sqrt(B*B - 4.0*A*C);


  mua2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B + J);
  mub2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B - J);

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
//      dDdX[i] = 0.5*(dRdX[i] - dLdX[i]);
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

//    dAdpsi = 2.0*scpsi;
//    dBdpsi = dAdpsi*(R*L - P*S);
//    dAdpsi *=(S - P);

    dFdpsi = 2.0*(dAdpsi - dBdpsi);
    dGdpsi = 2.0*dAdpsi - dBdpsi + (smu/J)*(B*dBdpsi - 2.0*C*dAdpsi);

    my_mu.dmudtheta = (0.5/my_mu.mu)*(dHdF*dFdpsi + dHdG*dGdpsi);

    dmudw = 0.0;
    for(int i=0; i<ncomps; i++) dmudw += dmudX[i]*dXdw[i] + dmudY[i]*dYdw[i];
    
    my_mu.dmudom = dmudw;
  
    D_mu2S = D / (mu2 - S);
    gamma = 1;// FAKENUMBERS /** \todo FIX!!!! */
    calc_type calc_n = (calc_type) n;
    omega_n = -1.0 * calc_n * my_const.omega_ce/gamma;
    //temporaries for simplicity

    term1 = (mu2* s2psi - P)/(mu2);
    
    denom = pow(D_mu2S*term1, 2) + c2psi*pow((P/mu2), 2);
    
    bessel_arg = calc_n* std::tan(psi)*tan(alpha) * (w - omega_n)/omega_n;// n x tan alpha (om - om_n)/om_n
    
    tmp_besp = boost::math::cyl_bessel_j(abs(n)+1, bessel_arg);
    tmp_besm = boost::math::cyl_bessel_j(abs(n)-1, bessel_arg);

    term2 = (1 + D_mu2S)*tmp_besp;
    term2 += (1 - D_mu2S)*tmp_besm;
    
    //tmp_bes = boost::math::cyl_bessel_j(abs(n), bessel_arg);
    tmp_bes = 0.5*bessel_arg *(tmp_besp + tmp_besm)/calc_n;
    //Use bessel identity to save time.
    // J_(n-1) + J_(n+1) = (2 n / arg) J_n
    term3 = scpsi*tmp_bes/std::tan(alpha);

    my_mu.phi = std::pow((0.5*term1*term2 + term3), 2)/denom;

  }

  
  return my_mu;
}

mu_dmudom plasma::get_high_dens_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n, bool Righthand)const{
  /** \brief Solve plasma dispersion and extensions
*
*Duplicates plasma::get_phi_mu_omega() but using reduced form of Stix parameters corresponding to a high-density assumption. This is mainly for comparison with the exact solution to validate this assumption.
*Also needs particle pitch angle alpha \todo Fix relativistic gamma... \todo Multispecies???? */

  calc_type term1, term2, term3, denom, tmp_bes, tmp_besp, tmp_besm, bessel_arg, D_mu2S, gamma, w2, w3;

  mu_dmudom my_mu;

  calc_type dndr[ncomps], dndth[ncomps];
  calc_type dB0dr, dB0dth;
  
#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(psi < 0 || psi >= pi) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, pitch angle (psi="+mk_str(psi)+") out of range!!!!!!", 0);
  if(alpha < 0 || alpha >= pi) my_error_print("!!!!!!!!Error in get_high_dens_phi_mu_om, particle pitch angle (alpha="+mk_str(alpha)+") out of range!!!!!!", 0);
  //I don't think there's an upper or lower bound on w we need to enforce
  //Nor any bounds on n
#endif
  
  w2 = w*w;
  w3 = w2*w;
  
  /** \todo get or calc these...*/
  // FAKENUMBERS
  for(int i=0;i<ncomps; i++){
    dndr[i] = 0.0;
    dndth[i] = 0.0;
  }
  dB0dr = 0.0;
  dB0dth = 0.0;
  
  calc_type R=1.0, L=1.0, P=1.0, J, S, D, A, B, C, s2psi, c2psi, mua2, mub2, mu2, scpsi;
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

  for(int i=0; i<1; ++i){
    //consider electron component only...
    R = 0.0 - wp2[i]/(w*(w + wc[i]));
    L = 0.0 - wp2[i]/(w*(w - wc[i]));
    P = 0.0 - wp2[i]/w2;
  }
  
  S = 0.5*(R + L);
  D = 0.5*(R - L);
  A = S*s2psi + P*c2psi;
  B = R*L*s2psi + P*S*(1.0+c2psi);
  C = P*R*L;
  J = std::sqrt(B*B - 4.0*A*C);


  mua2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B + J);
  mub2 = 1.0 - 2.0*(A - B + C)/(2.0*A - B - J);

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

    for(int i=0; i<1; i++){
      dPdX[i] = -1.0;
      dLdX[i] = -1.0/(1.0 - Y[i]);

      dRdX[i] = -1.0/(1.0 + Y[i]);
      dSdX[i] = 0.5*(dRdX[i] + dLdX[i]);
//      dDdX[i] = 0.5*(dRdX[i] - dLdX[i]);
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

//    dAdpsi = 2.0*scpsi;
//    dBdpsi = dAdpsi*(R*L - P*S);
//    dAdpsi *=(S - P);

    dFdpsi = 2.0*(dAdpsi - dBdpsi);
    dGdpsi = 2.0*dAdpsi - dBdpsi + (smu/J)*(B*dBdpsi - 2.0*C*dAdpsi);

    my_mu.dmudtheta = (0.5/my_mu.mu)*(dHdF*dFdpsi + dHdG*dGdpsi);

    dmudw = 0.0;
    for(int i=0; i<1; i++) dmudw += dmudX[i]*dXdw[i] + dmudY[i]*dYdw[i];
    
    my_mu.dmudom = dmudw;
  
    D_mu2S = D / (mu2 - S);
    gamma = 1;// FAKENUMBERS /** \todo FIX!!!! */
    calc_type calc_n = (calc_type) n;
    omega_n = -1.0 * calc_n * my_const.omega_ce/gamma;
    //temporaries for simplicity

    term1 = (mu2* s2psi - P)/(mu2);
    
    denom = pow(D_mu2S*term1, 2) + c2psi*pow((P/mu2), 2);
    
    bessel_arg = calc_n* std::tan(psi)*tan(alpha) * (w - omega_n)/omega_n;// n x tan alpha (om - om_n)/om_n
    
    tmp_besp = boost::math::cyl_bessel_j(abs(n)+1, bessel_arg);
    tmp_besm = boost::math::cyl_bessel_j(abs(n)-1, bessel_arg);

    term2 = (1 + D_mu2S)*tmp_besp;
    term2 += (1 - D_mu2S)*tmp_besm;
    
    //tmp_bes = boost::math::cyl_bessel_j(abs(n), bessel_arg);
    tmp_bes = 0.5*bessel_arg *(tmp_besp + tmp_besm)/calc_n;
    //Use bessel identity to save time.
    // J_(n-1) + J_(n+1) = (2 n / arg) J_n
    term3 = scpsi*tmp_bes/std::tan(alpha);

    my_mu.phi = std::pow((0.5*term1*term2 + term3), 2)/denom;

  }

  
  return my_mu;
}

std::vector<calc_type> plasma::get_resonant_omega(calc_type x, calc_type v_par, calc_type n)const{
/** \brief Solve plasma dispersion and doppler resonance simultaneously
*
*Obtains solutions of the Doppler resonance condition omega - k_par v_par = -n Omega_ce and a high-density approximation to the Whistler mode dispersion relation simultaneously. Assumes pure electron-proton plasma and uses cubic_solve. If solutions are found, they're returned in vector, otherwise empty vector is returned. @param x Wave normal angle @param v_par Particle velocity to solve with @param n Resonance number.
  \todo Extend to general case?
*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(v_par >= v0) my_error_print("!!!!!!!!Error in get_resonant_omega, velocity (v_par="+mk_str(v_par)+") out of range!!!!!!", 0);
#endif

  std::vector<calc_type> ret_vec;
  calc_type om_ce = this->get_omega_ref("ce");
  calc_type om_pe_loc = this->get_omega_ref("pe");
  calc_type om_ce_ref = this->get_omega_ref("c0");

  if(std::abs(v_par) < tiny_calc_type){
    if( std::abs(n)-1.0 < tiny_calc_type && std::abs(n) == 1.0) ret_vec.push_back(om_ce*n);
    return ret_vec;
  }
  else if(std::abs(n) < tiny_calc_type){
  //also special...
  //coeff d in cubic is 0 and we reduce to quadratic assuming omega != 0
    
    return ret_vec;
    //if omega_pe > omega_ce no solution...???
  /** \todo add case */
  
  }

  //Equation to solve is cubic of the form
  //a x^3 + b x^2 + c x + d = 0
  //or equivalently
  //x^3 + an x^2 + bn x + cn = 0

  calc_type a, b, c, d;
  calc_type an, bn, cn;
  calc_type cos_th = std::cos(std::atan(x));
  calc_type vel = v_par / v0;
  //For clarity
  calc_type vel_cos = std::pow(vel * cos_th, 2);

  calc_type gamma, gamma_sq;
  gamma_sq = 1.0/( 1.0 - std::pow(vel, 2));
  gamma = std::sqrt(gamma_sq);

  //Calculate coefficients
  //To maintain best precision we solve for x = omega/omega_ce_ref so x ~ 1
  a = (vel_cos - 1.0) * gamma_sq;
  b = (vel_cos*cos_th*gamma_sq + 2.0*gamma*n - gamma_sq*cos_th)*om_ce/om_ce_ref;
  c = ((2.0*gamma*n*cos_th - n*n)* std::pow(om_ce/om_ce_ref, 2) - std::pow(om_pe_loc/om_ce_ref, 2)*vel_cos*gamma_sq);
  d = -n*n*std::pow(om_ce/om_ce_ref, 3)*cos_th;
  an = b/a;
  bn = c/a;
  cn = d/a;
  
  ret_vec = cubic_solve(an, bn, cn);

  //Restore om_ce_re factor and delete any entries > om_ce as these aren't whistler modes
  for(size_t i=0; i<ret_vec.size(); ++i) ret_vec[i] *= om_ce_ref;
  for(size_t i=0; i<ret_vec.size(); ++i){
  /** \todo Rather than remove, don't ever add*/
    if(std::abs(ret_vec[i]) > std::abs(om_ce)){
      ret_vec.erase(ret_vec.begin() + i);
      --i;
    }
  }

  return ret_vec;

}

calc_type plasma::get_dispersion(my_type in, int wave_type, bool reverse, bool deriv, my_type theta)const{
/** \brief Solve analytic dispersion (approx)
*
* By default returns omega for a given k (see reverse and deriv params param). Uses local reference cyclotron and plasma frequencies and works with UNNORMALISED quantitites. NB: parameters out of range will silently return 0. NB: For Whistler modes this is an approximation and intended to be perfectly reversible. @param k Wavenumber @param wave_type wave species (see support.h) @param reverse Return k for input omega @param deriv Whether to instead return anayltic v_g  \todo Any other modes?*/

#ifdef DEBUG_ALL
  //Argument preconditions. Check only in debug mode for speed
  if(wave_type < WAVE_WHISTLER || wave_type > WAVE_X) my_error_print("!!!!!!!!Error in get_dispersion, wave_type unknown (wave_type="+mk_str(wave_type)+")!!!!!!", 0);
  //Theta is not conditioned because we reduce it below. This may change though
  //In is either a frequency or a wavenumber and thus has no limits, although may wish to add sanity limits
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
      break;//Case break

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
  }

  return ret;

}

