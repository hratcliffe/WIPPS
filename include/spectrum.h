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
*Holds data on the omega and angle distributions. If fed an FFTd data array this will be X^2(omega, theta) where X is E or B. The latter can depend on omega! Can be created/destroyed only by controllers, so has no public constructor/destructors. IMPORTANT: because we are working with FFT data, we assume the angle/frequency axis either covers some small cutout in +ve domain, or is symmetrical in positive and negative values. A few of the specific routines here use this to simplify things. The sign of omega is simply copied from the sign of k. The "angle" axis is stored as tan(theta) for theta the wave normal angle. Access to elements should use the wrappers at the bottom of spectrum.h, described in \ref spectAcc because internal layout could change in future
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
  bool g_is_angle_only;/**< Whether we have g(omega, x) (false) or just g(x) (true) */
  int function_type;/**< Type code for angular function. See support.h */
  size_t smooth;/**<Smoothing applied to B_omega, if any*/
  my_type norm_B;/**< Norm of B(w)*/
  my_type * norm_g;/**< Norms of g_w(x) for each w. Note even if g itself indep of omega, this is NOT*/

/********Basic setup and allocation functions ****/
  void construct();
  void init();
  explicit spectrum();
  explicit spectrum(int n_om, int n_ang, bool separable);
  explicit spectrum(std::string filename);
  ~spectrum();

/********Technical stuff making my_array a proper "object" ****/
  spectrum & operator=(const spectrum& src);
  spectrum(const spectrum &src);
  spectrum(spectrum && src) = default;/**<Move a spectrum object*/

/********Setup helper functions ****/
  void make_angle_axis();
  bool make_angle_distrib(my_type std_dev);

/******** Access helper functions ****/
  int where_omega(my_type omega);

/********Access wrappers ****/
/** \ingroup spectAcc 
*@{ */
  /** \brief Get frequency axis 
  *
  *NB do not muck with this pointer. It's provided for ease of using where but is a const my_type * for a reason
  @param[out] len Length of axis
  @return Pointer to omega axis */
  const my_type * get_omega_axis(size_t &len){return B_omega_array.get_axis(0, len);}
  /** \brief Get angle axis 
  *
  *NB do not muck with this pointer. It's provided for ease of using where but is a const my_type * for a reason
  @param[out] len Length of axis
  @return Pointer to angle axis */
  const my_type * get_angle_axis(size_t &len){return g_angle_array.get_axis(1, len);}
/** @} */

public:

  char block_id[ID_SIZE]; /**< The field name id from SDF file*/
  my_type time[2];/**< Time range over which data are taken*/
  size_t space[2];/**< Space range over which data are taken*/
  int wave_id; /**< ID for which wave mode cutout we're going for. See support.h*/
  bool get_g_is_angle_only(){return g_is_angle_only;}/**< Get flag showing if spectrum is separable @return True if g is a function of only angle (and not frequency) false else*/

/********Technical stuff making my_array a proper "object" ****/
  bool operator==(const spectrum &rhs)const;
  bool operator!=(const spectrum &rhs)const{return !(*this == rhs);}/**< See spectrum::operator==()*/

/********Setup helper functions ****/
  inline bool is_good()const{return (B_omega_array.is_good() && g_angle_array.is_good() && norm_g);}/**<Check if a spectrum is complete and useable @return Boolean true if good, false else*/
#ifdef RUN_TESTS_AND_EXIT
  void make_test_spectrum(int angle_type=FUNCTION_DELTA, bool two_sided=false, my_type om_ce=17000.0, my_type std_dev = DEFAULT_SPECTRUM_ANG_STDDEV);
#endif
  bool generate_spectrum(data_array& parent, int om_fuzz=10, int angle_type=FUNCTION_DELTA, my_type std_dev = DEFAULT_SPECTRUM_ANG_STDDEV, data_array * mask=nullptr);

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

/********Spectrum operation helpers ****/
  bool calc_norm_B();
  bool calc_norm_g(size_t om_ind);
  /** Get the normalising constant for B part of spectrum 
  @return Current value of norm_B*/
  my_type get_norm_B(){return norm_B;}
  /** Get the normalising constant for g part of spectrum 
  @param om_ind Frequency index to get from, 0 for separable spectra
  @return Current value of norm_g at specified location */
  my_type get_norm_g(size_t om_ind){return (om_ind < get_omega_length())? norm_g[om_ind]:0.0;}
  
/********File IO ****/
  bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file);
  
/********Data release (for testing) ****/
  data_array  copy_out_B();
  data_array  copy_out_g();
  
/********Access wrappers ****/
/** \ingroup halp
*\defgroup spectAcc Spectrum access wrappers
*\brief Accessors for the two parts of spectrum, B and g
*
*Spectrum does not guarantee the internal representation of the B and g parts, so these should be used to get/set the data and axes for B and g parts
*@{ */

  inline my_type get_B_element(size_t n_om)const{return B_omega_array.get_element(n_om);}/**< Get B element frequency index n_om*/
  inline my_type get_g_element(size_t n_ang)const{return g_angle_array.get_element((size_t) 0, n_ang);}/**< Get g element at angle n_ang (separable spectra)*/
  inline my_type get_g_element(size_t n_om, size_t n_ang)const{return g_angle_array.get_element(n_om, n_ang);}/**<Get g element at frequency n_om, angle n_ang (nonseparable spectra)*/
  
  inline void set_B_element(size_t n_om, my_type val){B_omega_array.set_element(n_om, val);}/**< Set B element frequency index n_om*/
  inline void set_g_element(size_t n_ang, my_type val){g_angle_array.set_element(0, n_ang, val);}/**< Set g element at angle n_ang (separable spectra)*/
  inline void set_g_element(size_t n_om, size_t n_ang, my_type val){g_angle_array.set_element(n_om, n_ang, val);}/**<Set g element at frequency n_om, angle n_ang (nonseparable spectra)*/

  inline my_type get_om_axis_element(size_t nx)const{return B_omega_array.get_axis_element(0, nx);}/**< Get frequency axis element at index nx*/
  inline my_type get_ang_axis_element(size_t nx)const{return g_angle_array.get_axis_element(1, nx);}/**< Get angle axis element at index nx*/

  inline long get_om_axis_index_from_value(my_type omega)const{return B_omega_array.get_axis_index_from_value(0, omega);}/**< Get frequency axis index for value omega*/
  inline long get_ang_axis_index_from_value(my_type ang)const{return g_angle_array.get_axis_index_from_value(1, ang);}/**< Get angle axis index for value ang*/

  inline void set_om_axis_element(size_t nx, my_type val){B_omega_array.set_axis_element(0, nx, val);g_angle_array.set_axis_element(0, nx, val);}/**< Set frequency axis index for value omega*/
  inline void set_ang_axis_element(size_t nx, my_type val){g_angle_array.set_axis_element(1, nx, val);}/**< Set angle axis index for value ang*/

  inline size_t get_g_dims()const{return this->g_angle_array.get_dims();}/**< Get rank of g array*/
  inline size_t get_g_dims(size_t i)const{return this->g_angle_array.get_dims(i);}/**< Get dimensions of g array*/
  inline size_t get_B_dims()const{return this->B_omega_array.get_dims();}/**<Get rank of B array */
  inline size_t get_B_dims(size_t i)const{return this->B_omega_array.get_dims(i);}/**<Get dimensions of B array*/
  inline size_t get_angle_length()const{return this->g_angle_array.get_dims(1);}/**<Get length of angle axis*/
  inline size_t get_omega_length()const{return this->B_omega_array.get_dims(0);}/**<Get length of frequency axis*/
  
/** @} */
/** @} */

#ifdef RUN_TESTS_AND_EXIT
  //Allow deep testing
  friend class test_entity_spectrum;
#endif
  
};

/********Main spectral calculations ****/
/** \ingroup halp
*\defgroup spectGs Spectrum calculations
*\brief Calculate normalised spectra
*
*To calculate diffusion coefficients we need the spectrum at a point, with proper normalisation. These calculate the factors called G_1 and G_2 in Albert \cite Albert2005
*@{ */

calc_type get_G1(spectrum * my_spect, calc_type omega);
calc_type get_G2(spectrum * my_spect, calc_type omega, calc_type x);
/** @} */

#endif
