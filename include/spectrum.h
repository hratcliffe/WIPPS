//
//  spectrum.h
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#ifndef _spectrum_h
#define _spectrum_h


//Lets have a spectrum class then
//Contains data, axis, sizes, ids (field, time range, space range)
//For now, let's stick to the 1-d, then there's only delta say in angle...
//We can impose say a narrow Gaussian later.
//NB Note in 1-d also omega and k_pllel and k are freely interchangeable
//Alright. Lets just add an extra row for the B^2, and fill the rest of the array with g_omega. It's never going to be a massive amount of data
//But if the angle profile is say a function, we don't need to...

class data_array;
class plasma;
class controller;

class spectrum : public data_array{
//specialised data array with extra ID fields to hold additional data

  friend void controller::add_spectrum(int nx, int n_ang);
  friend void controller::add_spectrum(int * row_lengths, int ny);
  friend controller::~controller();

  controller * my_controller;
  bool ax_omega;/** Flag whether we derived in k or omega*/
  void construct();
  spectrum(int nx, int n_ang);
  spectrum(int * row_lengths, int ny);
  my_type normB;
  my_type* normg;
  bool normaliseB();
  bool normaliseg(my_type omega);
  virtual ~spectrum();

public:

  int wave_id;
  //ID for which wave cutout we're going for...

  bool angle_is_function;
  //Says either we're in 1-D, and we impose functional form for anguar profile, e.g delta function or Gaussian, or we're in 2-D but still wish to impose.

  int function_type;

  int n_angs;

  void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type=FUNCTION_NULL);

  bool generate_spectrum(data_array * parent);

  my_type get_omega(my_type k, int wave_type, bool deriv=0);
  my_type get_k(my_type omega, int wave_type, bool deriv =0);

  my_type * get_angle_distrib(int &len, my_type omega=0.0);

//  int where(my_type * ax_ptr, int len, my_type target, std::function<bool(my_type,my_type)> func = std::greater<my_type>());
  int where(my_type * ax_ptr, int len, my_type target);
  std::vector<int> all_where(my_type * ax_ptr, int len, my_type target, std::function<bool(my_type,my_type)> func = std::greater<my_type>());
  
  bool write_to_file(std::fstream &file);

  void make_test_spectrum();

  calc_type get_G1(calc_type omega);
  calc_type get_G2(calc_type omega, calc_type x);

};


#endif
