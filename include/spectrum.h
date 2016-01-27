//
//  spectrum.h
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#ifndef _spectrum_h
#define _spectrum_h

class data_array;
class plasma;
class controller;

/** \brief Specialised data_array to hold spectrum
*
*Specialises shape and adds functions to process spectrum, normalise it etc. Can be created/destroyed only by controllers \author Heather Ratcliffe \date 24/09/2015
*/
class spectrum : public data_array{
  friend void controller::add_spectrum(int nx, int n_ang);
  friend void controller::add_spectrum(int * row_lengths, int ny);
  friend controller::~controller();

  controller * my_controller;/**< Links this to a plasma object*/
  bool ax_omega;/**< Flag whether we derived in k or omega*/
  void construct();
  spectrum(int nx, int n_ang);/**< Private because only controllers can create/destroy*/
  spectrum(int * row_lengths, int ny);/**< Private because only controllers can create/destroy*/
  my_type normB;/**< Norm of B(w)*/
  my_type* normg;/**< Norms of g_w(x) for each w*/
  bool normaliseB();/**< Fills normB*/
  bool normaliseg(my_type omega);/**< Fills normg for omega*/
  virtual ~spectrum();/**< Private because only controllers can create/destroy*/

public:

  int wave_id; /**< ID for which wave mode cutout we're going for. See support.h*/
  bool angle_is_function;/**< Says we impose g(x) rather than have one g for each w*/
  int function_type;/**< Type code for angular function. See support.h */
  int n_angs;/**< Number of angles to use*/

  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type=FUNCTION_NULL);

  bool generate_spectrum(data_array * parent, int om_fuzz=10);

  my_type get_omega(my_type k, int wave_type, bool deriv=0);
  my_type get_k(my_type omega, int wave_type, bool deriv =0);

  my_type * get_angle_distrib(int &len, my_type omega=0.0);

  std::vector<int> all_where(my_type * ax_ptr, int len, my_type target, std::function<bool(my_type,my_type)> func = std::greater<my_type>());
  
  bool write_to_file(std::fstream &file);

  void make_test_spectrum(int time[2], int space[2]);

  calc_type get_G1(calc_type omega);
  calc_type get_G2(calc_type omega, calc_type x);

};


#endif
