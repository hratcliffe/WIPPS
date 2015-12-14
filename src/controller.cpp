//
//  controller.cpp
//  
//
//  Created by Heather Ratcliffe on 19/11/2015.
//
//

#include <stdio.h>
#include <cmath>
#include "support.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"

extern deck_constants my_const;

controller::controller(){
/** \brief Setup
*
*\todo Plasma object should be setup from files or such
*/
  my_plas = new plasma(my_const.omega_ce * me/std::abs(q0));
  my_spect = nullptr;
  my_d = nullptr;

};

controller::~controller(){

  if(my_plas) delete my_plas;
  if(my_spect) delete my_spect;
  if(my_d) delete my_d;
};

void controller::add_spectrum(int nx, int n_ang){
/** \brief Create and add spectrum
*
*rectangular spectrum, for when angular distribution depends on omega
*/
  my_spect = new spectrum(nx, n_ang);
  my_spect->my_controller = this;
}

void controller::add_spectrum(int * row_lengths, int ny){
/** \brief Create and add spectrum
*
*Spectrum for when angular distribution does not depend on omega
*/

  my_spect = new spectrum(row_lengths, ny);
  my_spect->my_controller = this;
  
}

void controller::add_d(int nx, int n_angs){
/** \brief Create and addd diffusion_coefficient and create it's basic axes etc
*
*
*/

  my_d = new diffusion_coeff(nx, n_angs);
  my_d->my_controller = this;
  my_d->make_velocity_axis();
  my_d->make_pitch_axis();

}
/*
A bounce average function here would: assemble ordered stack of d's and spectrums.
Run calculate on each, matching proerply
*** If required, supply v_par axes to them so they line up as wanted
Bounce average into global D interpolating v as necessary etc
\todo Write this...
*/

