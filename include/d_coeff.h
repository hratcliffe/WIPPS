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

enum class D_type_spec {alpha_alpha, alpha_p, p_alpha, p_p};

/** \brief Progress info structure*/
struct running_report{
  size_t last_report;
  size_t report_interval;
  bool quiet;
};
  

/** \brief Diffusion coefficient object
*
* Specialised data_array containing the calculated coefficient plus relevant ids. Can be made/destroyed only by controller object. In general should be a 2-D array with first axis momentum, second pitch-angle
@author Heather Ratcliffe @date 23/09/2015 \ingroup cls
*/
class diffusion_coeff: public data_array{
private:

  int n_thetas; /**< Number of wave normal angles to consider for integrals*/
  int n_n; /**< Max number of resonances to consider */
  bool single_n;/**< Flag to force considering only a single resonance, mostly for test purposes*/
  int n_used;/**< The resonance to consider if single_n is set, mostly for test purposes*/
  controller * my_controller;/**< Owning controller which gives access to plasma and spectrum*/
  friend class controller;
  explicit diffusion_coeff(int n_momenta, int n_angs);
  virtual ~diffusion_coeff(){;};

  int get_min_n(calc_type mod_v, calc_type cos_alpha, my_type k_thresh, calc_type om_ce);
  int get_max_n(calc_type mod_v, calc_type cos_alpha, my_type k_thresh, calc_type om_ce);

  void make_velocity_axis();
  void make_pitch_axis();
  void copy_ids(spectrum * src);
  void first_running_report(size_t total_its, running_report &rept, bool quiet);
  void do_running_report(size_t v_ind, calc_type mod_v, size_t min_n, size_t max_n, running_report &rept);

public:

  int latitude;/**< Latitude of calculation*/
  int wave_id;/**< ID of wave mode considered*/
  std::string tag; /**<Identifies as local, averaged etc*/
  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[ID_SIZE]);
  /** \brief Set single resonance to consider*/
  void set_single_n(int n){single_n = true; n_used = n;}
  /** \brief Set max/min resonant number to consider*/
  void set_max_n(int n){n_n = n;}
  bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file);

  d_report calculate(D_type_spec type_of_D = D_type_spec::alpha_alpha, bool quiet=0);

  //--- Helpers for accessing and values. Because we might want a tan or plain angle axis and we don't want that to contaminate the rest of the code
  my_type get_element_by_values(my_type p, my_type alpha);
  my_type get_axis_element_ang(size_t ind);
  inline my_type angle_to_stored_angle(my_type alpha){return(alpha);}/**<Convert actual angle in radians to stored angle axis value @param alpha Actual angle value @return Converted value in axis*/
  inline my_type stored_angle_to_angle(my_type value){return(value);}/**<Convert stored angle axis value to actual angle in radians @param value Stored angle value @return Actual angle value */
//  inline my_type angle_to_stored_angle(my_type alpha){return tan(alpha);}/**<Convert actual angle in radians to stored angle axis value @param alpha Actual angle value @return Converted value in axis*/
//  inline my_type stored_angle_to_angle(my_type alpha){return atan(alpha);}/**<Convert stored angle axis value to actual angle in radians @param value Stored angle value @return Actual angle value */
};



#endif /* defined(____d_coeff__) */
