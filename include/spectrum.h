//
//  spectrum.h
//  
//
//  Created by Heather Ratcliffe on 29/07/2016  as refactor of spectrum class
//
//

#ifndef _spectrum_h
#define _spectrum_h

#include "support.h"
#include "controller.h"
#include "my_array.h"
#include "data_array.h"


class plasma;
class controller;

/** \brief A spectrum in omega and angle
*
*Holds data on the omega and angle distributions. The latter can depend on omega! Can be created/destroyed only by controllers. IMPORTANT: because we are working with FFT data, we assume the angle/frequency axis either covers some small cutout in +ve domain, or is symmetrical in positive and negative values. A few of the specific routines here use this to simplify things. The sign of omega is simply copied from the sign of k. The "angle" axis is stored as tan(theta) for theta the wave normal angle. Access to elements should use the wrappers at the bottom of the file as internal layout may change
  \author Heather Ratcliffe \date 24/09/2015
*/
class spectrum{
  friend void controller::add_spectrum(int nx, int n_ang, bool separable);
  friend void controller::add_spectrum(std::string file);
  friend controller::~controller();

  controller * my_controller;/**< Links this to a plasma object*/
  void construct();
  spectrum();
  spectrum(int nx, int n_ang, bool separable);/**< Private because only controllers can create/destroy*/
  spectrum & operator=(const spectrum& src);/**<For setting one spectrum equal another*/
  spectrum(const spectrum &src);/**<For copying a spectrum*/
  spectrum(spectrum && src) = default;
  my_type normB;/**< Norm of B(w)*/
  my_type* normg;/**< Norms of g_w(x) for each w*/
  bool normaliseB();/**< Fills normB*/
  bool normaliseg(my_type omega);/**< Fills normg for omega*/
  my_type max_power;/**<Value of maximum in spectral power*/
  bool angle_is_function;/**< Says we impose g(x) rather than have one g for each w*/
  int function_type;/**< Type code for angular function. See support.h */
  size_t smooth;

  data_array g_angle_array;
  data_array B_omega_array;

public:

  spectrum(std::string filename);
  ~spectrum();/**< Reloading a spectrum can be done by anyone. Creating a new is restricted*/

  char block_id[ID_SIZE]; /**< The field name id from SDF file*/

  my_type time[2];/**< Time range over which data are taken*/
  size_t space[2];/**< Space range over which data are taken*/
  int wave_id; /**< ID for which wave mode cutout we're going for. See support.h*/
  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type=FUNCTION_DELTA);
  void init();
  void smooth_B(int n_pts);
  bool generate_spectrum(data_array& parent, int om_fuzz=10, int angle_type=FUNCTION_DELTA);
  bool truncate_om(my_type om_min, my_type om_max);
  bool truncate_x(my_type x_min, my_type x_max);

  my_type get_omega(my_type k, int wave_type, bool deriv=0);
  my_type get_k(my_type omega, int wave_type, bool deriv =0);

  bool make_angle_distrib();
  int where_omega(my_type value);
  std::vector<int> all_where(my_type * ax_ptr, int len, my_type target, std::function<bool(my_type,my_type)> func = std::greater<my_type>());
  
  bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file);
  
  void make_test_spectrum(int angle_type=FUNCTION_DELTA, bool two_sided=false, my_type om_ce=17000.0);

  calc_type get_G1(calc_type omega);
  calc_type get_G2(calc_type omega, calc_type x);
  calc_type check_upper();
  calc_type get_peak_omega();
  
  void copy_ids(data_array & src);
  bool check_ids( data_array & src);
  inline bool is_good(){return (B_omega_array.is_good() && g_angle_array.is_good() && normg);}
  
  
  //The following are all wrappers around the getter/setters for the g and B arrays. USE THEM!
  inline my_type get_B_element(size_t n_om)const{return B_omega_array.get_element(n_om);}
  inline my_type get_g_element(size_t n_ang)const{return g_angle_array.get_element((size_t) 0, n_ang);}
  inline my_type get_g_element(size_t n_om, size_t n_ang)const{return g_angle_array.get_element(n_om, n_ang);}
  inline void set_B_element(size_t n_om, my_type val){B_omega_array.set_element(n_om, val);}
  inline void set_g_element(size_t n_ang, my_type val){g_angle_array.set_element(0, n_ang, val);}
  inline void set_g_element(size_t n_om, size_t n_ang, my_type val){g_angle_array.set_element(n_om, n_ang, val);}

  inline my_type get_om_axis_element(size_t nx)const{return B_omega_array.get_axis_element(0, nx);}
  inline my_type get_ang_axis_element(size_t nx)const{return g_angle_array.get_axis_element(1, nx);}

  inline void set_om_axis_element(size_t nx, my_type val){B_omega_array.set_axis_element(0, nx, val);g_angle_array.set_axis_element(0, nx, val);}
  inline void set_ang_axis_element(size_t nx, my_type val){g_angle_array.set_axis_element(1, nx, val);}

  inline size_t get_g_dims()const{return this->g_angle_array.get_dims();}
  inline size_t get_g_dims(size_t i)const{return this->g_angle_array.get_dims(i);}
  inline size_t get_B_dims()const{return this->B_omega_array.get_dims();}
  inline size_t get_B_dims(size_t i)const{return this->B_omega_array.get_dims(i);}
  
  data_array  copy_out_B();
  data_array  copy_out_g();
  
};



#endif
