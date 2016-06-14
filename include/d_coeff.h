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
#include "controller.h"
#include "spectrum.h"
#include "my_array.h"
#include "plasma.h"


class data_array;
class spectrum;
class plasma;
class controller;

/** \brief Diffusion coefficient object
*
* Specialised data_array containing the calculated coefficient plus relevant ids. Can be made/destroyed only by controller object. 
@author Heather Ratcliffe @date 23/09/2015.
*/
class diffusion_coeff: public data_array{

private:

  int n_thetas; /**< Number of wave normal angles to consider */
  int n_n; /**< Max number of resonances to consider */
  controller * my_controller;/**< Owning controller which gives access to plasma and spectrum*/
  friend void controller::add_d(int nx, int n_angs, int pos);
  friend controller::~controller();
  friend void controller::handle_d_mpi();
  diffusion_coeff(int nx, int n_ang);
  virtual ~diffusion_coeff();


public:

  int latitude;/**< Latitude of calculation*/
  int wave_id;/**< ID of wave mode considered*/
  std::string tag; /**<Identifies as local, averaged etc*/
  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10]);

  bool write_to_file(std::fstream &file);

  void calculate();
  void make_velocity_axis();
  void make_pitch_axis();
  /** \todo Write these eh*/
  int get_min_n(calc_type v_par, my_type k_thresh, calc_type om_ce);
  int get_max_n(calc_type v_par, my_type k_thresh, calc_type om_ce);
  //Return min and max n worth exploring for largest k_par in range
};



#endif /* defined(____d_coeff__) */
