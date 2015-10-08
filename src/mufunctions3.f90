!***************************************************************************
!> This program calculates all the functions of mu that are needed for the 
!! RHS of the ODEs. Note that mu is calculated using the Appleton-Hartree 
!! relation, and the choice of sign is obtained from Albert [2005].
!! Replaces implementations in \ref mufunctions.f90 and \ref mufunctions2.f90
!! @author  Clare E. J. Watt                                   @date 18/05/10

subroutine mufunctions(th,w,ncomps,wp,wp2,wc,wc2,pcharge,pmass,psi, &
&                      s2psi,c2psi,dndr,dndth,dB0dr,dB0dth,mu,mug, &
&                      dmudr,dmudth,alpha,merror)
!
!  This program calculates all the functions of mu that are needed for the 
!  RHS of the ODEs. Note that mu is calculated using the Appleton-Hartree 
!  relation, and the choice of sign is obtained from Albert [2005].
!
!***************************************************************************
!
 use nrtype; use parameters
 implicit none
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
 allocate(dmudX(ncomps)); allocate(dXdr(ncomps)); allocate(dXdth(ncomps))
 allocate(dmudY(ncomps)); allocate(dYdr(ncomps)); allocate(dYdth(ncomps))
 allocate(dPdX(ncomps)); allocate(dLdX(ncomps)); allocate(dRdX(ncomps))
 allocate(dDdX(ncomps)); allocate(dSdX(ncomps))
 allocate(dAdX(ncomps)); allocate(dBdX(ncomps)); allocate(dCdX(ncomps))
 allocate(dFdX(ncomps)); allocate(dGdX(ncomps))
 allocate(dLdY(ncomps)); allocate(dRdY(ncomps))
 allocate(dDdY(ncomps)); allocate(dSdY(ncomps))
 allocate(dAdY(ncomps)); allocate(dBdY(ncomps)); allocate(dCdY(ncomps))
 allocate(dFdY(ncomps)); allocate(dGdY(ncomps))
 allocate(dXdw(ncomps)); allocate(dYdW(ncomps))
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
 X = wp2/(w*w)
 Y = wc/w
 S = half*(R + L)
 D = half*(R - L)
 A = S*s2psi + P*c2psi
 B = R*L*s2psi + P*S*(one+c2psi)
 C = P*R*L
 J = sqrt(B*B - 4.0_dp*A*C)
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
 deallocate(dmudX); deallocate(dXdr); deallocate(dXdth)
 deallocate(dmudY); deallocate(dYdr); deallocate(dYdth)
 deallocate(dPdX); deallocate(dLdX); deallocate(dRdX)
 deallocate(dSdX); deallocate(dDdX)
 deallocate(dLdY); deallocate(dRdY)
 deallocate(dSdY); deallocate(dDdY)
 deallocate(dAdX); deallocate(dBdX); deallocate(dCdX)
 deallocate(dFdX); deallocate(dGdX)
 deallocate(dAdY); deallocate(dBdY); deallocate(dCdY)
 deallocate(dFdY); deallocate(dGdY)
 deallocate(dXdw); deallocate(dYdw)
!
end subroutine mufunctions
