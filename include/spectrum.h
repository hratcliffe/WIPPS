//
//  spectrum2.h
//  
//
//  Created by Heather Ratcliffe on 29/07/2016  as refactor of spectrum class
//
//

#ifndef _spectrum2_h
#define _spectrum2_h

#include "support.h"
#include "controller.h"
#include "my_array.h"


class data_array;
class plasma;
class controller;

/** \brief Specialised data_array to hold spectrum
*
*Specialises shape and adds functions to process spectrum, normalise it etc. Can be created/destroyed only by controllers. IMPORTANT: because we are working with FFT data, we assume the angle/frequency axis either covers some small cutout in +ve domain, or is symmetrical in positive and negative values. A few of the specific routines here use this to simplify things. The "angle" axis is stored as tan(theta) for theta the wave normal angle. We seperate forwards and backwards wave modes by \author Heather Ratcliffe \date 24/09/2015 
*/
class spectrum{
  friend void controller::add_spectrum(int nx, int n_ang);
  friend void controller::add_spectrum(int * row_lengths, int ny);
  friend controller::~controller();

  controller * my_controller;/**< Links this to a plasma object*/
  void construct();
  spectrum(int nx, int n_ang);/**< Private because only controllers can create/destroy*/
  spectrum(int * row_lengths, int ny);/**< Private because only controllers can create/destroy*/
  my_type normB;/**< Norm of B(w)*/
  my_type* normg;/**< Norms of g_w(x) for each w*/
  bool normaliseB();/**< Fills normB*/
  bool normaliseg(my_type omega);/**< Fills normg for omega*/
  virtual ~spectrum();/**< Private because only controllers can create/destroy*/
  my_type max_power;/**<Value of maximum in spectral power*/
//  my_type k_thresh;/**<K value of where spectrum crosses threshold (noise) value*/
  data_array * g_angle_array;/**< \todo refactor to non-pointer*/
  data_array * B_omega_array;

public:
  char block_id[ID_SIZE]; /**< The field name id from SDF file*/

  float time[2];/**< Time range over which data are taken*/
  int space[2];/**< Space range over which data are taken*/
  int wave_id; /**< ID for which wave mode cutout we're going for. See support.h*/
  bool angle_is_function;/**< Says we impose g(x) rather than have one g for each w*/
  int function_type;/**< Type code for angular function. See support.h */
  int n_angs;/**< Number of angles to use*/
  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type=FUNCTION_DELTA);

  bool generate_spectrum(data_array * parent, int om_fuzz=10, int angle_type=FUNCTION_DELTA);
  bool truncate_om(my_type om_min, my_type om_max);
  bool truncate_x(my_type x_min, my_type x_max);

  my_type get_omega(my_type k, int wave_type, bool deriv=0);
  my_type get_k(my_type omega, int wave_type, bool deriv =0);

  bool make_angle_distrib();
  int where_omega(my_type value);
  std::vector<int> all_where(my_type * ax_ptr, int len, my_type target, std::function<bool(my_type,my_type)> func = std::greater<my_type>());
  
  bool write_to_file(std::fstream &file);

  void make_test_spectrum(int time[2], int space[2],int angle_type=FUNCTION_DELTA);

  calc_type get_G1(calc_type omega);
  calc_type get_G2(calc_type omega, calc_type x);
  calc_type check_upper();
  calc_type get_peak_omega();
  
  void copy_ids(data_array * src);
  bool check_ids( data_array * src);
  bool is_good(){return B_omega_array->is_good() && g_angle_array->is_good();}
  
  my_type get_B_element(int nx);
  my_type get_ang_element(int nx, int ny=0);
  my_type get_B_axis_element(int nx);
  my_type get_ang_axis_element(int nx);
  int get_ang_dims(int i=-1);
  int get_B_dims(int i=-1);

  
};


#endif
