//Removed code doing the full mu calculation. This may be needed, so is worth keeping around, but is uneeded as now
/** \brief Full refractive index
*
*Contains refractive index mu and all derivatives, plus error flag
*/
struct mu{
  calc_type mu; /**< Refractive index */
  calc_type mug; /**< */
  calc_type dmudr; /**< d mu /dr (radial distance) */
  calc_type dmudlat; /**< d mu/ d latitude */
  calc_type dmudtheta; /**< d mu / d theta (wave normal angle) */
  calc_type dmudom; /**< d mu / d omega (wave frequency)*/
  calc_type alpha; /**< */
  int err; /**< 0 if mu found successfully, 1 else*/

};


mu plasma::get_root(calc_type th, calc_type w, calc_type psi, bool Righthand){
/** \brief Solve plasma dispersion 
*
*Solves Appleton-Hartree plasma dispersion and returns struct containing mu, its derivatives and error code. See \ref str
*Duplicated from mufunctions by CEJ Watt
*
*@param th Polar coordinate, or latitude @param w Wave frequency @param psi Wave pitch angle k to B0 @param Righthand Mode handedness
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
  //I don't think there's an upper or lower bound on w we need to enforce, only positivity. We'll take abs below and succeed, but in debug mode also warn
  if(w < 0) my_error_print("!!!!!!Error in get_root, wave frequency (w="+mk_str(w)+") is negative!!!!!!", 0);

#endif

  w = std::abs(w);
  
  for(int i=0;i<ncomps; i++){
    dndr[i] = 0.0;
    dndth[i] = 0.0;
  }
  dB0dr = 0.0;
  dB0dth = .0;
  
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

