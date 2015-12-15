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
  current_spect=0;
  current_d=0;

};

controller::~controller(){

  if(my_plas) delete my_plas;
  
  for(int i=0; i<my_spect.size(); i++){
    delete my_spect[i];
    my_spect[i] = nullptr;
  }
  for(int i=0; i<my_d.size(); i++){
    delete my_d[i];
    my_d[i] = nullptr;
  }

};

void controller::add_spectrum(int nx, int n_ang){
/** \brief Create and add spectrum
*
*rectangular spectrum, for when angular distribution depends on omega
*/
  spectrum * tmp_spect;
  tmp_spect = new spectrum(nx, n_ang);
  tmp_spect->my_controller = this;
  my_spect.push_back(tmp_spect);
  current_spect = my_spect.size()-1;

}

void controller::add_spectrum(int * row_lengths, int ny){
/** \brief Create and add spectrum
*
*Spectrum for when angular distribution does not depend on omega
*/
  spectrum * tmp_spect;
  tmp_spect = new spectrum(row_lengths, ny);
  tmp_spect->my_controller = this;
  my_spect.push_back(tmp_spect);
  current_spect = my_spect.size()-1;
  
}

void controller::add_d(int nx, int n_angs){
/** \brief Create and addd diffusion_coefficient and create it's basic axes etc
*
*
*/
  diffusion_coeff * tmp_d;
  tmp_d = new diffusion_coeff(nx, n_angs);
  tmp_d->my_controller = this;
  tmp_d->make_velocity_axis();
  tmp_d->make_pitch_axis();
  my_d.push_back(tmp_d);
  current_d = my_d.size()-1;

}

spectrum * controller::get_current_spectrum(){

  if(!my_spect.empty()) return my_spect[current_spect];
  else return nullptr;

}

diffusion_coeff * controller::get_current_d(){
  
  if(!my_d.empty()) return my_d[current_d];
  else return nullptr;
}

void controller::bounce_average(){
  /** \todo finish this!!*/
  if(my_d.size() ==1) return;
  //No averaging to do!

  int dims[2];
  get_size(dims);
  add_d(dims[0], dims[1]);
  //Add new averaged D
  
  for(int i=0; i<my_d.size() -1; i++){
    //Flatten down each d onto this one
  
  
  }
  

}

void controller::get_size(int dims[2]){

  dims[0] = my_d[current_d]->get_dims(0);
  dims[1] = my_d[current_d]->get_dims(1);

}
