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

/** Ties together plasma objects with resulting spectra and diffusion coefficients */
class controller{
/** */
  plasma * my_plas;

  //These might become vectors
  spectrum * my_spect;
  diffusion_coeff * my_d;

public:


  controller();
  ~controller();
  void add_spectrum(int nx, int n_ang);
  void add_spectrum(int * row_lengths, int ny);
  void add_d(int nx, int n_angs);
  spectrum * get_current_spectrum(){return my_spect;};
  diffusion_coeff * get_current_d(){return my_d;};
  /** For Future expansion and to prevent public access*/
  plasma * get_plasma(){return my_plas;};
};


#endif
