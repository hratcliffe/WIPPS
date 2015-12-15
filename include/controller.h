//
//  controller.h
//  
//
//  Created by Heather Ratcliffe on 19/11/2015.
//
//

#ifndef _controller_h
#define _controller_h

class spectrum;
class plasma;
class diffusion_coeff;

/** \brief Controls plasma, spectrum and d_coeff objects and their connections
*
*This is the public facing class controlling plasma, spectrum and d_coeff objects. It makes sure there is a plasma to provide needed functions for the latters. Will also be responislbe for forming bounce averaged coefficients from a stack of ds?
* @author Heather Ratcliffe @date 19/11/2015.
*/

class controller{
  plasma * my_plas; /**< Plasma object*/

  std::vector<spectrum *> my_spect;/**< Spectrum object or vector thereof?*/
  std::vector<diffusion_coeff *> my_d;/**< Diffusion coefficient object or vector thereof?*/
  int current_spect;
  int current_d;
public:

  controller();
  ~controller();
  void add_spectrum(int nx, int n_ang);
  void add_spectrum(int * row_lengths, int ny);
  void add_d(int nx, int n_angs);
  spectrum * get_current_spectrum();
  diffusion_coeff * get_current_d();
  plasma * get_plasma(){return my_plas;};
};


#endif
