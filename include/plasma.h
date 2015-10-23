//
//  plasma.h
//  
//
//  Created by Heather Ratcliffe on 07/10/2015.
//
//

#ifndef _plasma_h
#define _plasma_h

#include<vector>

class plasma{

static const int ncomps=4; 

int wave_type;

calc_type * position;
calc_type * ne;
calc_type B0;
calc_type * om_pe;
calc_type om_ce;

/*calc_type * dndr[ncomps];
calc_type * dndth[ncomps];
calc_type * dB0dr;
calc_type * dB0dth;
*/
//Might not need both options. Note also the deck constants versions...

//Depending on calls, might want to store the derivs and stuff

//Needs to get a spatial axis to work with.
//Needs to hold density profile, background B field profile.

//Other Plasma params, such as species massm charge
calc_type pmass[ncomps];
calc_type pcharge[ncomps];
calc_type pdens[ncomps];

mu last_mu;
calc_type last_th, last_w, last_psi;

public:

plasma();
void get_density(){;}
void get_B0(){;}

//obtain these profiles somehow, either from deck or otherwise. ne might be constant... B can get from file 0

//If we copy mufunctions we have something like
//We'll want an averaged thing over the spatial block we're working with, because we're assuming waves from whole block are identical.
//But one call per ptich angle and frequency might turn out to be time consuming. Profiling necessary...
mu get_root(calc_type th, calc_type w, calc_type psi);
calc_type get_phi( calc_type th, calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n=0);
mu_dmudom get_phi_mu_om( calc_type th, calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n=0);

/*probably want some parameter "which_thing" is some way of specifying what we want, probably via enum or named constants PROBABLY a mask is best, i.e. bitmask with named contants and we get each thing we want. But how to return? Rturn invalid number for those we on't request. Make it optional so by default we get all? Or we return all of: real(DP),intent(out)::mu,mug                   !< refractive index, group
 real(DP),intent(out)::dmudr,dmudth             !< derivs
 real(DP),intent(out)::alpha                    !< angle between k and prop
 integer(I4B),intent(out)::merror               !< testing for mu
*/

calc_type get_omega(calc_type x, calc_type v_par, calc_type n);


};

#endif
