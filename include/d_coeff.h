//
//  d_coeff.h
//  
//
//  Created by Heather Ratcliffe on 23/09/2015.
//
//

#ifndef ____d_coeff__
#define ____d_coeff__

#include <stdio.h>

class data_array;
class spectrum;
class plasma;
class controller;

class diffusion_coeff: public data_array{

private:

  int n_thetas;
  int n_n;
  controller * my_controller;
  friend void controller::add_d(int nx, int n_angs);
  friend controller::~controller();

  diffusion_coeff(int nx, int n_ang);
  virtual ~diffusion_coeff();


public:

  int latitude;
  int wave_id;
  //ID for which wave cutout we're going for...

  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10]);

  bool write_to_file(std::fstream &file);

  void calculate();
  void make_velocity_axis();
  
};



#endif /* defined(____d_coeff__) */
