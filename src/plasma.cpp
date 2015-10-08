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

#include "support.h"
#include "plasma.h"

extern deck_constants my_const;

plasma::plasma(){

//allocate position, density, B axes
//set other parameters

  ncomps = 1;
  pmass = 1.0 * mp;
  // H plasma
  pcharge = 1.0 * q0;

}

//std::vector<calc_type>
//Note regardless of ncomps we only return one collective value...
mu plasma::get_root(calc_type th, calc_type w, calc_type psi){
/** Duplicated from mufunctions by CEJ Watt
*
*@param th Polar coordinate @param w Wave frequency @param psi Wave pitch angle k to B0
*/
  mu mu_ret;
  calc_type R=1.0, L=1.0, P=1.0, S, D, A, B, C, J, s2psi, c2psi, mua2, mub2, mu2;
  calc_type F, G, smu;
  calc_type dHdF, dHdG, dmudw, dmudpsi,dAdpsi,dBdpsi,dFdpsi,dGdpsi,dpsidth;
  calc_type wp[ncomps], wp2[ncomps], wc[ncomps], X[ncomps], Y[ncomps];
  
  calc_type dmudX[ncomps], dXdr[ncomps], dXdth[ncomps], dmudY[ncomps], dYdr[ncomps], dYdth[ncomps];
  calc_type dPdX[ncomps], dLdX[ncomps], dRdX[ncomps], dSdX[ncomps], dDdX[ncomps];
  calc_type dAdX[ncomps], dBdX[ncomps], dCdX[ncomps], dFdX[ncomps], dGdX[ncomps];
  calc_type dLdY[ncomps], dRdY[ncomps], dSdY[ncomps], dDdY[ncomps];
  calc_type dAdY[ncomps], dBdY[ncomps], dCdY[ncomps], dFdY[ncomps], dGdY[ncomps], dXdw[ncomps], dYdw[ncomps];
  
  //These will hold suitably calc'd plasma frequency, square and cyclotron freq. If we have to derive from position, we do...
  for(int i=0; i<ncomps; ++i){
    wp[i] = my_const.omega_pe;
    wp2[i] = wp[i]*wp[i];
    wc[i] = my_const.omega_ce;
    
    X[i] = wp2[i]/(w*w);
    Y[i] = wc[i]/w;

  }
//We loop over components and trust compiler to unroll for us :) most of these will trivially vectorise anyway.

  s2psi = pow(sin(psi), 2);
  c2psi = pow(cos(psi), 2);
 /** \todo Check these are correct interpretation... */

  for(int i=0; i<ncomps; i++){
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
  mu_ret.dmudth = 0.0;
  mu_ret.alpha = 0.0;
  mu_ret.err = 1;
  
  //Now we fill...
  
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
}
    
    mu_ret.err = 0;
  
  }
  /*
    dPdX = -one
    dLdX = -one/(one - Y)
    dRdX = -one/(one + Y)
    dSdX = half*(dRdX + dLdX)
    dDdX = half*(dRdX - dLdX)
    dAdX = half*(dRdX + dLdX)*s2psi + dPdX*c2psi
    dBdX = (L*dRdX + R*dLdX)*s2psi + (P*dSdX + S*dPdX)*(one + c2psi)
    dCdX = P*R*dLdX + P*L*dRdX + R*L*dPdX
    dFdX = two*(dAdX - dBdX + dCdX)
    dGdX = two*dAdX - dBdX + (smu/J)*(B*dBdX - two*(A*dCdX + C*dAdX)) 
    dmudX = (half/mu)*(dHdF*dFdX + dHdG*dGdX)
!
    dRdY = X/((one + Y)**2)
    dLdY = -X/((one - Y)**2)
    dSdY = half*(dRdY + dLdY)
    dDdY = half*(dRdY - dLdY)
    dAdY = dSdY*s2psi
    dBdY = (L*dRdY + R*dLdY)*s2psi + P*dSdY*(one + c2psi)
    dCdY = P*(L*dRdY + R*dLdY)
    dFdY = two*(dAdY - dBdY + dCdY)
    dGdY = two*dAdY - dBdY + (smu/J)*(B*dBdY - two*(A*dCdY + C*dAdY)) 
    dmudY = (half/mu)*(dHdF*dFdY + dHdG*dGdY)
!
    dXdr = ev*ev*dndr/(w*w*e0*pmass)
    dYdr = pcharge*dB0dr/(w*pmass)
    dXdth = ev*ev*dndth/(w*w*e0*pmass)
    dYdth = pcharge*dB0dth/(w*pmass)
    dXdw = -two*wp2/(w*w*w)
    dYdw = -wc/(w*w)
!
    dAdpsi = sin(two*psi)*(S - P)
    dBdpsi = sin(two*psi)*(R*L - P*S)
    dFdpsi = two*(dAdpsi - dBdpsi)
    dGdpsi = two*dAdpsi - dBdpsi + (smu/J)*(B*dBdpsi - two*C*dAdpsi)    
    dmudpsi = (half/mu)*(dHdF*dFdpsi + dHdG*dGdpsi)
!
    dpsidth = -two/(one + 3.0_dp*cos(th)*cos(th))
    dmudr = zero
    dmudw = zero
    dmudth = dmudpsi*dpsidth
    do i=1,ncomps
       dmudr = dmudr + dmudX(i)*dXdr(i) + dmudY(i)*dYdr(i)
       dmudth = dmudth + dmudX(i)*dXdth(i) + dmudY(i)*dYdth(i)
       dmudw = dmudw + dmudX(i)*dXdw(i) + dmudY(i)*dYdw(i)
    enddo
    mug = mu + w*dmudw
    alpha = -dmudpsi/mu
!
 endif


*/
  return mu_ret;
}


/*
subroutine mufunctions(th,w,ncomps,wp,wp2,wc,wc2,pcharge,pmass,psi, &
&                      s2psi,c2psi,dndr,dndth,dB0dr,dB0dth,mu,mug, &
&                      dmudr,dmudth,alpha,merror)
!
!
 real(DP),intent(in)::th                        !< polar coordinate
 real(DP),intent(in)::w                         !< wave frequency
 integer(I4B),intent(in)::ncomps                !< number of components
 real(DP),dimension(ncomps),intent(in)::wp,wp2  !< plasma frequency, sq
 real(DP),dimension(ncomps),intent(in)::wc,wc2  !< gyrofrequency, sq
 real(DP),dimension(ncomps),intent(in)::pcharge !< particle charge
 real(DP),dimension(ncomps),intent(in)::pmass   !< particle mass
 real(DP),intent(in)::psi                       !< angle between k & B0
 real(DP),intent(in)::s2psi,c2psi               !< trig psi
 real(DP),dimension(ncomps),intent(in)::dndr    !< deriv of n wrt r
 real(DP),dimension(ncomps),intent(in)::dndth   !< deriv of n wrt th
 real(DP),intent(in)::dB0dr,dB0dth              !< deriv of B0 wrt r,th
 
 
 real(DP),intent(out)::mu,mug                   !< refractive index, group
 real(DP),intent(out)::dmudr,dmudth             !< derivs
 real(DP),intent(out)::alpha                    !< angle between k and prop
 integer(I4B),intent(out)::merror               !< testing for mu
 
 
 real(DP),dimension(:),allocatable::dmudX,dXdr,dXdth,dmudY,dYdr,dYdth
 real(DP),dimension(:),allocatable::dPdX,dLdX,dRdX,dSdX,dDdX
 real(DP),dimension(:),allocatable::dAdX,dBdX,dCdX,dFdX,dGdX
 real(DP),dimension(:),allocatable::dLdY,dRdY,dSdY,dDdY
 real(DP),dimension(:),allocatable::dAdY,dBdY,dCdY,dFdY,dGdY
 real(DP),dimension(:),allocatable::dXdw,dYdw
 real(DP)::dHdG,dHdF,dmudw
 real(DP)::dmudpsi,dAdpsi,dBdpsi,dFdpsi,dGdpsi,dpsidth
 real(DP)::mu2                                  ! mu squared
 real(DP)::mua2,mub2
 real(DP)::S,D,R,L,P                            ! Stix parameters 
 real(DP)::A,B,C,J                              ! For calculation of mu2
 real(DP)::F,G                                  ! For calculation of mu2
 real(DP),dimension(ncomps)::X,Y                ! ratio of (wp/w)^2 and wc/w
 integer(I4B)::i                                ! count variable
 integer(I4B)::smu
!
! Allocate all allocatable arrays
!
!
! Calculate mu
! 
 R = one
 L = one
 P = one
 do i=1,ncomps
    R = R - wp2(i)/(w*(w+wc(i)))
    L = L - wp2(i)/(w*(w-wc(i)))
    P = P - wp2(i)/(w*w)
 enddo 
////////////////////////

 X = wp2/(w*w)
 Y = wc/w
 S = half*(R + L)
 D = half*(R - L)
 A = S*s2psi + P*c2psi
 B = R*L*s2psi + P*S*(one+c2psi)
 C = P*R*L
 J = sqrt(B*B - 4.0_dp*A*C)
 //////////////////
 
 mua2 = one - two*(A - B + C)/(two*A - B + J)
 mub2 = one - two*(A - B + C)/(two*A - B - J)
 if ((mua2<zero).and.(mub2<zero)) then
    merror = 1
    mu = one                              ! These are just placeholders
    mug = zero
    dmudr = zero
    dmudth = zero
    alpha = zero    
 else 
    merror = 0
    if (D<0) then                              ! see Albert [2005]
       smu = +1
       mu2 = mua2
    else
       smu = -1
       mu2 = mub2
    endif
    mu = sqrt(mu2)
    F = two*(A - B + C)
    G = two*A - B + smu*J
    dHdF = -one/G
    dHdG = F/(G*G)
!
    dPdX = -one
    dLdX = -one/(one - Y)
    dRdX = -one/(one + Y)
    dSdX = half*(dRdX + dLdX)
    dDdX = half*(dRdX - dLdX)
    dAdX = half*(dRdX + dLdX)*s2psi + dPdX*c2psi
    dBdX = (L*dRdX + R*dLdX)*s2psi + (P*dSdX + S*dPdX)*(one + c2psi)
    dCdX = P*R*dLdX + P*L*dRdX + R*L*dPdX
    dFdX = two*(dAdX - dBdX + dCdX)
    dGdX = two*dAdX - dBdX + (smu/J)*(B*dBdX - two*(A*dCdX + C*dAdX)) 
    dmudX = (half/mu)*(dHdF*dFdX + dHdG*dGdX)
!
    dRdY = X/((one + Y)**2)
    dLdY = -X/((one - Y)**2)
    dSdY = half*(dRdY + dLdY)
    dDdY = half*(dRdY - dLdY)
    dAdY = dSdY*s2psi
    dBdY = (L*dRdY + R*dLdY)*s2psi + P*dSdY*(one + c2psi)
    dCdY = P*(L*dRdY + R*dLdY)
    dFdY = two*(dAdY - dBdY + dCdY)
    dGdY = two*dAdY - dBdY + (smu/J)*(B*dBdY - two*(A*dCdY + C*dAdY)) 
    dmudY = (half/mu)*(dHdF*dFdY + dHdG*dGdY)
!
    dXdr = ev*ev*dndr/(w*w*e0*pmass)
    dYdr = pcharge*dB0dr/(w*pmass)
    dXdth = ev*ev*dndth/(w*w*e0*pmass)
    dYdth = pcharge*dB0dth/(w*pmass)
    dXdw = -two*wp2/(w*w*w)
    dYdw = -wc/(w*w)
!
    dAdpsi = sin(two*psi)*(S - P)
    dBdpsi = sin(two*psi)*(R*L - P*S)
    dFdpsi = two*(dAdpsi - dBdpsi)
    dGdpsi = two*dAdpsi - dBdpsi + (smu/J)*(B*dBdpsi - two*C*dAdpsi)    
    dmudpsi = (half/mu)*(dHdF*dFdpsi + dHdG*dGdpsi)
!
    dpsidth = -two/(one + 3.0_dp*cos(th)*cos(th))
    dmudr = zero
    dmudw = zero
    dmudth = dmudpsi*dpsidth
    do i=1,ncomps
       dmudr = dmudr + dmudX(i)*dXdr(i) + dmudY(i)*dYdr(i)
       dmudth = dmudth + dmudX(i)*dXdth(i) + dmudY(i)*dYdth(i)
       dmudw = dmudw + dmudX(i)*dXdw(i) + dmudY(i)*dYdw(i)
    enddo
    mug = mu + w*dmudw
    alpha = -dmudpsi/mu
!
 endif
!
! Deallocate all the arrays
!
!
end subroutine mufunctions
*/