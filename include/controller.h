//
//  controller.h
//  
//
//  Created by Heather Ratcliffe on 19/11/2015.
//
//

#ifndef _controller_h
#define _controller_h

#include "support.h"

class spectrum;
class plasma;
class diffusion_coeff;
/*Circular dependencies, don't include headers*/

/** \brief Controls plasma, spectrum and d_coeff objects and their connections
*
*This is the public facing class controlling plasma, spectrum and d_coeff objects. It makes sure there is a plasma to provide needed functions for the latters. Will also be responislbe for forming bounce averaged coefficients from a stack of ds?
* @author Heather Ratcliffe @date 19/11/2015.
*/

class controller{
  plasma * my_plas; /**< Plasma object*/

  std::vector<spectrum *> my_spect;/**< Spectrum object or vector thereof?*/
  std::vector<diffusion_coeff *> my_d;/**< Diffusion coefficient object or vector thereof?*/
  size_t current_spect;
  size_t current_d;
  void get_size(int dims[2]);
public:

  controller(std::string file_prefix);
  ~controller();
  void add_spectrum(int nx, int n_ang,bool separable);
  void add_d(int nx, int n_angs, int pos=-1);
  spectrum * get_current_spectrum();
  diffusion_coeff * get_current_d();
  plasma * get_plasma(){return my_plas;};
  void bounce_average();
  void handle_d_mpi();
  bool save_spectra(std::string pref);
  bool save_D(std::string pref);

};


#endif
