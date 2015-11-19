//
//  controller.cpp
//  
//
//  Created by Heather Ratcliffe on 19/11/2015.
//
//

#include <stdio.h>
#include "support.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"

controller::controller(){

  my_plas = new plasma();
  my_spect = nullptr;
  my_d = nullptr;


};


controller::~controller(){

  if(my_plas) delete my_plas;
  if(my_spect) delete my_spect;
  if(my_d) delete my_d;
};

void controller::add_spectrum(int nx, int n_ang){

  my_spect = new spectrum(nx, n_ang);
  my_spect->my_controller = this;
}

void controller::add_spectrum(int * row_lengths, int ny){

  my_spect = new spectrum(row_lengths, ny);
  my_spect->my_controller = this;
  
}

void controller::add_d(int nx, int n_angs){

  my_d = new diffusion_coeff(nx, n_angs);
  my_d->my_controller = this;

}
