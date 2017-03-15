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
#include "data_array.h"
#include "plasma.h"


class spectrum;
class plasma;
class controller;

/** \brief Diffusion coefficient object
*
* Specialised data_array containing the calculated coefficient plus relevant ids. Can be made/destroyed only by controller object. In general should be a 2-D array with first axis momentum, second pitch-angle \todo Consider locking away the angle axis, and/or moving away from tan
@author Heather Ratcliffe @date 23/09/2015 \ingroup cls
*/
class diffusion_coeff: public data_array{

private:

  int n_thetas; /**< Number of wave normal angles to consider for integrals*/
  int n_n; /**< Max number of resonances to consider */
  controller * my_controller;/**< Owning controller which gives access to plasma and spectrum*/
  friend class controller;
  explicit diffusion_coeff(int n_momenta, int n_ang);
  virtual ~diffusion_coeff(){;};

  int get_min_n(calc_type v_par, my_type k_thresh, calc_type om_ce);
  int get_max_n(calc_type v_par, my_type k_thresh, calc_type om_ce);

  void make_velocity_axis();
  void make_pitch_axis();
  void copy_ids(spectrum * spect);

public:

  int latitude;/**< Latitude of calculation*/
  int wave_id;/**< ID of wave mode considered*/
  std::string tag; /**<Identifies as local, averaged etc*/
  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[ID_SIZE]);

  bool write_to_file(std::fstream &file);

  d_report calculate(bool quiet=0);

  //--- Helpers for accessing and values. Because we might want a tan or plain angle axis and we don't want that to contaminate the rest of the code
  my_type get_element_by_values(my_type p, my_type alpha);
  my_type get_axis_element_ang(size_t ind);
  inline my_type angle_to_stored_angle(my_type alpha){return(alpha);}
  inline my_type stored_angle_to_angle(my_type alpha){return(alpha);}
//  inline my_type angle_to_stored_angle(my_type alpha){return tan(alpha);}
//  inline my_type stored_angle_to_angle(my_type alpha){return atan(alpha);}
};



#endif /* defined(____d_coeff__) */
