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
#ifdef RUN_TESTS_AND_EXIT
#include "tests.h"
#endif

/** \brief A spectrum in omega and angle
*
*Holds data on the omega and angle distributions. The latter can depend on omega! Can be created/destroyed only by controllers, so has not public constructor/destructors. IMPORTANT: because we are working with FFT data, we assume the angle/frequency axis either covers some small cutout in +ve domain, or is symmetrical in positive and negative values. A few of the specific routines here use this to simplify things. The sign of omega is simply copied from the sign of k. The "angle" axis is stored as tan(theta) for theta the wave normal angle. Access to elements should use the wrappers at the bottom of spectrum.h, described in \ref spectAcc because internal layout could change in future
  \author Heather Ratcliffe \date 24/09/2015 \ingroup cls
*/
class spectrum{
  friend class controller;/**<Controllers can create/destroy spectra and access their internals*/
  controller * my_controller;/**< Links this to a plasma object*/

/******** The data ****/
  data_array g_angle_array;/**< Angular component*/
  data_array B_omega_array;/**< Magnitude component*/

/********Tags and info ****/
  my_type max_power;/**<Value of maximum in spectral power*/
  bool angle_is_function;/**< Whether we have g(omega, x) (false) or just g(x) (true) */
  int function_type;/**< Type code for angular function. See support.h */
  size_t smooth;/**<Smoothing applied to B_omega, if any*/
  my_type normB;/**< Norm of B(w)*/
  my_type * normg;/**< Norms of g_w(x) for each w*/

/********Basic setup and allocation functions ****/
  void construct();
  void init();
  explicit spectrum();
  explicit spectrum(int nx, int n_ang, bool separable);
  explicit spectrum(std::string filename);
  ~spectrum();

/********Technical stuff making my_array a proper "object" ****/
  spectrum & operator=(const spectrum& src);
  spectrum(const spectrum &src);
  spectrum(spectrum && src) = default;/**<Move a spectrum object*/
  bool operator==(const spectrum &rhs)const;
  bool operator!=(const spectrum &rhs)const{return !(*this == rhs);}/**< See spectrum::operator==()*/

/********Setup helper functions ****/
  void make_angle_axis();
  bool make_angle_distrib();

/******** Access helper functions ****/
  int where_omega(my_type value);

/********Spectrum operation helpers ****/
  bool normaliseB();
  bool normaliseg(my_type omega);

/********Access wrappers ****/
/** \ingroup spectAcc 
*@{ */
  my_type * get_omega_axis(size_t &len){return B_omega_array.get_axis(0, len);}
  my_type * get_angle_axis(size_t &len){return g_angle_array.get_axis(1, len);}
/** @} */

public:

  char block_id[ID_SIZE]; /**< The field name id from SDF file*/
  my_type time[2];/**< Time range over which data are taken*/
  size_t space[2];/**< Space range over which data are taken*/
  int wave_id; /**< ID for which wave mode cutout we're going for. See support.h*/

/********Setup helper functions ****/
  inline bool is_good()const{return (B_omega_array.is_good() && g_angle_array.is_good() && normg);}/**<Check if a spectrum is complete and useable*/
#ifdef RUN_TESTS_AND_EXIT
  void make_test_spectrum(int angle_type=FUNCTION_DELTA, bool two_sided=false, my_type om_ce=17000.0);
#endif
  bool generate_spectrum(data_array& parent, int om_fuzz=10, int angle_type=FUNCTION_DELTA, data_array * mask=nullptr);

  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type=FUNCTION_DELTA);
  void copy_ids(const data_array & src);
  bool check_ids(const data_array & src)const;
  void copy_tags(const spectrum & src);
  bool check_tags(const spectrum & src)const;

/******** Access helper functions ****/
  my_type get_omega(my_type k, int wave_type, bool deriv=0,my_type theta=0.0);
  my_type get_k(my_type omega, int wave_type, bool deriv =0,my_type theta=0.0);

/********Spectrum operation helpers ****/
  void smooth_B(int n_pts);
  bool truncate_om(my_type om_min, my_type om_max);
  bool truncate_x(my_type x_min, my_type x_max);
  calc_type check_upper();
  calc_type get_peak_omega();

/********File IO ****/
  bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file);
  
/********Main spectral calculations ****/
  calc_type get_G1(calc_type omega);
  calc_type get_G2(calc_type omega, calc_type x);
  
/********Data release (for testing) ****/
  data_array  copy_out_B();
  data_array  copy_out_g();
  
/********Access wrappers ****/
/** \defgroup spectAcc Spectrum access wrappers
*\brief Accessors for the two parts of spectrum, B and g
*
*Spectrum does not guarantee the internal representation of the B and g parts, so these should be used to get/set the data and axes for B and g parts
*@{ */

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
/** @} */

#ifdef RUN_TESTS_AND_EXIT
  //Allow deep testing
  friend class test_entity_spectrum;
#endif
  
};



#endif
