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

class controller{

public:
// Ties togther the sets of plasma, spectrum etc objects.
//If we only ever make and delete via this, we're golden.
  plasma * my_plas;

  //These might become vectors
  spectrum * my_spect;
  diffusion_coeff * my_d;

  controller();
  ~controller();
  void add_spectrum(int nx, int n_ang);
  void add_spectrum(int * row_lengths, int ny);
  void add_d(int nx, int n_angs);
};


#endif
