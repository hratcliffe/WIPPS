//
//  plasma.cpp
//  
//
//  Created by Heather Ratcliffe on 07/10/2015.
/** \brief Plasma parameters and dispersion
*
*This class will take care of solving the plasma dispersion roots, selecting a wave mode if necessary etc etc. The base calculations are modified from file mufunctions3.f90 author  Clare E. J. Watt                                   date 18/05/10. From there: Note that mu is calculated using the Appleton-Hartree relation, and the choice of sign is obtained from Albert [2005].
* @author Heather Ratcliffe  @date 07/10/2015
*/

#include <stdio.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <boost/math/special_functions/sign.hpp>
#include "support.h"
#include "plasma.h"

extern deck_constants my_const;

plasma::plasma(){

  pmass[0] = me;
  pmass[1] = mp;
  pmass[2] = mp*4.0;
  pmass[3] = me*16.0;
  
  pcharge[0] = -1.0*q0;
  pcharge[1] = 1.0*q0;
  pcharge[2] = 1.0*q0;
  pcharge[3] = 1.0*q0;
  
  calc_type ref_dens = my_const.omega_pe * my_const.omega_pe * eps0 * me / q0/q0;

  pdens[0] = 1.0 * ref_dens;
  pdens[1] = 1.00 * ref_dens;
  pdens[2] = 0.0 * ref_dens;
  pdens[3] = 0.0*ref_dens;
  
  B0 = my_const.omega_ce * me/std::abs(q0);

  this->om_ce = (pcharge[0]) * this->B0 / pmass[0];
  /**reference electron cyclotron freq \todo FIX! FAKENUMBERS */
}
plasma::~plasma(){


}

mu plasma::get_root(calc_type th, calc_type w, calc_type psi){
/** Duplicated from mufunctions by CEJ Watt
*
*@param th Polar coordinate, or latitude @param w Wave frequency @param psi Wave pitch angle k to B0 \todo Check constant types OK \todo Better root selection
*
*On notation: within this routine and plasma::get_phi we use notation as from mufunctions3.f90. In the return values as defined in support.h we match with Lyons and Albert. Thus in my_mu, we have lat, r, theta, omega for polar coordinate, r, wave normal angle and wave frequency
*/
  mu mu_ret;
  
  std::cout<<"THIS DOES NOT WORK>>>>>>>>>>>>>>>"<<std::endl;
  //REWRITE THIS against get_phi_mu_om. THIS IS WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  calc_type dndr[ncomps], dndth[ncomps];
  calc_type dB0dr, dB0dth;
  
  /** \todo get or calc these...*/
  //std::cout<<"These numbers are wrong, remember"<<std::endl;
  // FAKENUMBERS
  for(int i=0;i<ncomps; i++){
    dndr[i] = 1.0;
    dndth[i] = 1.0;
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
    wp[i] = this->get_omega_ref("pe");
    wp2[i] = wp[i]*wp[i];
    wc[i] = this->get_omega_ref("ce");
    
    X[i] = wp2[i]/(w*w);
    Y[i] = wc[i]/w;

  }
//We loop over components and trust compiler to unroll for us :) most of these will trivially vectorise anyway.

  s2psi = pow(sin(psi), 2);
  c2psi = pow(cos(psi), 2);
 /** \todo Check these are correct interpretation... */

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
  
  if( (mua2 > 0.0) && (mub2 > 0.0) ){
  
    if(D < 0.0 ){ smu = 1.0; mu2 = mua2;} //see Albert [2005]
    else{smu = -1.0; mu2 = mub2;}
    
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
      dXdw[i] = -2.0*wp2[i]/(pow(w, 3));
      dYdw[i] = -wc[i]/(pow(w, 2));

    }
    
    dAdpsi = sin(2.0*psi)*(S - P);
    dBdpsi = sin(2.0*psi)*(R*L - P*S);
    dFdpsi = 2.0*(dAdpsi - dBdpsi);
    dGdpsi = 2.0*dAdpsi - dBdpsi + (smu/J)*(B*dBdpsi - 2.0*C*dAdpsi);
    mu_ret.dmudtheta = (0.5/mu_ret.mu)*(dHdF*dFdpsi + dHdG*dGdpsi);

    dpsidth = -2.0/(1.0 + 3.0*pow(cos(th), 2));
    dmudw = 0.0;
    
    mu_ret.dmudlat = mu_ret.dmudtheta*dpsidth;
    //even if this one can be folded into above, keep it out as not vectorisable
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

//  this->my_mu = mu_ret;
//  mu_set = true;
//store away mu for furture use, along with key params used. 
  last_mu = mu_ret;
  last_th = th;
  last_w = w;
  last_psi = psi;
  
  return mu_ret;
}

calc_type plasma::get_phi(calc_type th, calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n){
/**Get's the Phi defined by Lyons 1974.
*Will be clumsy for now, because each call recalls mu, and we need to sum over n in the end. And it duplicates the STIX params calcs
*Also needs particle pitch angle alpha \todo Fix relativistic gamma...
 */
 
  calc_type wp[ncomps], wp2[ncomps], wc[ncomps];
  calc_type R, L, P, S, D, ret, term1, term2, term3, denom, tmp_bes, bessel_arg, sin2psi, D_mu2S, gamma;

  for(int i=0; i<ncomps; ++i){
    wp[i] = my_const.omega_pe;
    wp2[i] = wp[i]*wp[i];
    wc[i] = my_const.omega_ce;
  }
  for(int i=0; i<ncomps; i++){
    R = R - wp2[i]/(w*(w + wc[i]));
    L = L - wp2[i]/(w*(w - wc[i]));
    P = P - wp2[i]/(w*w);
  }

  S = 0.5*(R + L);
  D = 0.5*(R - L);

  mu my_mu;

  if((last_th != th) || (last_w != w) || (last_psi != psi)){
    std::cout<<"Regetting"<<std::endl;
    my_mu = this->get_root(th, w, psi);}
  else{my_mu = last_mu;}
  
  calc_type mu2 = pow(my_mu.mu, 2);

  sin2psi = pow(sin(psi), 2);
  D_mu2S = D / (mu2 - S);
  gamma = 1;
  omega_n = -1.0 * n * my_const.omega_ce/gamma;
  //temporaries for simplicity

  term1 = (mu2* sin2psi - P)/(mu2);
  
  denom = pow(D_mu2S*term1, 2) + pow((P*cos(psi)/mu2), 2);
  
  bessel_arg = n* tan(psi)*tan(alpha) * (w - omega_n)/omega_n;// n x tan alpha (om - om_n)/om_n
  
  tmp_bes = boost::math::cyl_bessel_j(n+1, bessel_arg);
  term2 = (1 + D_mu2S)*tmp_bes;
  
  tmp_bes = boost::math::cyl_bessel_j(n-1, bessel_arg);
  term2 += (1 - D_mu2S)*tmp_bes;

  tmp_bes = boost::math::cyl_bessel_j(n, bessel_arg);
  term3 = sin(psi)*cos(psi)*tmp_bes/tan(alpha);

  ret = pow((0.5*term1*term2 + term3), 2)/denom;
  return ret;

}

mu_dmudom plasma::get_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n){
/**Get's the Phi defined by Lyons 1974, and mu, dmu/dom i.e. the set needed for D.
*Also needs particle pitch angle alpha \todo Fix relativistic gamma... WATCH for Clares version which uses a different alpha entirely...

 */
 
  calc_type term1, term2, term3, denom, tmp_bes, bessel_arg, sin2psi, D_mu2S, gamma;

  mu_dmudom my_mu;

  calc_type dndr[ncomps], dndth[ncomps];
  calc_type dB0dr, dB0dth;
  
  /** \todo get or calc these...*/
  //std::cout<<"These numbers are wrong, remember"<<std::endl;
  // FAKENUMBERS
  for(int i=0;i<ncomps; i++){
    dndr[i] = 0.0;
    dndth[i] = 0.0;
  }
  dB0dr = 0.0;
  dB0dth = 0.0;
  
  calc_type R=1.0, L=1.0, P=1.0, J, S, D, A, B, C, s2psi, c2psi, mua2, mub2, mu2;
  calc_type F, G, smu;
  calc_type dHdF, dHdG, dmudw;
  calc_type wp[ncomps], wp2[ncomps], wc[ncomps], X[ncomps], Y[ncomps];
  
  calc_type dmudX[ncomps], dmudY[ncomps];
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
  c2psi = std::pow(cos(psi), 2);
 /** \todo Check these are correct interpretation... */

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
  my_mu.mu = 1.0;
  my_mu.dmudom = 0.0;
  my_mu.err = 1;

  if( (mua2 > 0.0) || (mub2 > 0.0) ){
    if(D < 0.0 ){ smu = 1.0; mu2 = mua2;} //see Albert [2005]
    else{smu = -1.0; mu2 = mub2;}
    
    my_mu.mu = sqrt(mu2);
    
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
      dmudX[i] = (0.5/my_mu.mu)*(dHdF*dFdX[i] + dHdG*dGdX[i]);

      dRdY[i] = X[i]/( pow(1.0 + Y[i], 2) );
      dLdY[i] = -X[i]/( pow(1.0 - Y[i], 2) );
      dSdY[i] = 0.5*(dRdY[i] + dLdY[i]);
      dDdY[i] = 0.5*(dRdY[i] - dLdY[i]);

      dAdY[i] = dSdY[i]*s2psi;
      dBdY[i] = (L*dRdY[i] + R*dLdY[i])*s2psi + P*dSdY[i]*(1.0 + c2psi);
      dCdY[i] = P*(L*dRdY[i] + R*dLdY[i]);

      dFdY[i] = 2.0*(dAdY[i] - dBdY[i] + dCdY[i]);
      dGdY[i] = 2.0*dAdY[i] - dBdY[i] + (smu/J)*(B*dBdY[i] - 2.0*(A*dCdY[i] + C*dAdY[i]));

      dmudY[i] = (0.5/my_mu.mu)*(dHdF*dFdY[i] + dHdG*dGdY[i]);


      dXdw[i] = -2.0*wp2[i]/(pow(w, 3));
      dYdw[i] = -wc[i]/(pow(w, 2));

    }
    
    dmudw = 0.0;
    for(int i=0; i<ncomps; i++) dmudw = dmudw + dmudX[i]*dXdw[i] + dmudY[i]*dYdw[i];
    
    my_mu.dmudom = dmudw;
    my_mu.err = 0;
  
    sin2psi = std::pow(sin(psi), 2);
    D_mu2S = D / (mu2 - S);
    gamma = 1;
    omega_n = -1.0 * n * my_const.omega_ce/gamma;
    //temporaries for simplicity

    term1 = (mu2* sin2psi - P)/(mu2);
    
    denom = pow(D_mu2S*term1, 2) + pow((P*cos(psi)/mu2), 2);
    
    bessel_arg = n* tan(psi)*tan(alpha) * (w - omega_n)/omega_n;// n x tan alpha (om - om_n)/om_n
    
    tmp_bes = boost::math::cyl_bessel_j(abs(n)+1, bessel_arg);
    term2 = (1 + D_mu2S)*tmp_bes;
    
    tmp_bes = boost::math::cyl_bessel_j(abs(n)-1, bessel_arg);
    term2 += (1 - D_mu2S)*tmp_bes;

    tmp_bes = boost::math::cyl_bessel_j(abs(n), bessel_arg);
    term3 = sin(psi)*cos(psi)*tmp_bes/tan(alpha);

    my_mu.phi = pow((0.5*term1*term2 + term3), 2)/denom;

  }

  
  return my_mu;
}

std::vector<calc_type> plasma::get_omega(calc_type x, calc_type v_par, calc_type n){
/**Get resonant frequency for particular x...
*
*Solve high density approx to get omega. for pure electron proton plasma....
* Calls cubic_solve Note for slowly changing v_par, suggests Newtons method might be more efficient. Although much of this could be precomputed for given grids.
Return empty vector if no valid solutions \todo Extend to general case?
*/

  std::vector<calc_type> ret_vec;

  calc_type wc = this->get_omega_ref("ce");
  calc_type omega_pe_loc = this->get_omega_ref("pe");
  calc_type om_ref_ce = my_const.omega_ce;

  if(std::abs(v_par) < tiny_calc_type){
    if( std::abs(n)-1.0 < tiny_calc_type ) ret_vec.push_back(wc* n);
    return ret_vec;
  }
  else if(std::abs(n) < tiny_calc_type){
  //also special...
  //coeff d in cubic is 0 and we reduce to quadratic assuming omega != 0
    
    return ret_vec;
    //if omega_pe > omega_ce no solution...???
  /** \todo add case */
  
  }

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
  
  return ret_vec;

}

calc_type plasma::get_omega_ref(std::string code){

  if(code == "ce") return this->om_ce;
  if(code == "pe") return my_const.omega_pe;
  else return 0.0;


}

calc_type plasma::get_dispersion(my_type k, int wave_type){
/** \brief Gets omega for given k
*
* Uses local refernce cyclotron and plasma frequencies and works with UNNORMALISED quantitites
*/
  calc_type ret = 0.0;
  calc_type om_ce_loc, om_pe_loc;
  
  om_ce_loc = this->get_omega_ref("ce");
  om_pe_loc = this->get_omega_ref("pe");
  
  switch(wave_type){

    case WAVE_WHISTLER :
      ret = v0*v0*k*k*om_ce_loc/(v0*v0 + om_pe_loc*om_pe_loc)/om_ce_loc;
      
      //here goes dispersion in suitable normed units.
      break;

    case WAVE_PLASMA :
      ret =1 ;
      break;
  }

  return ret;

}

