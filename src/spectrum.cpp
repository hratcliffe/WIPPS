//
//  spectrum.cpp
//  
//
//  Created by Heather Ratcliffe on 29/07/2016 as refactor of spectrum class
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "plasma.h"
#include "spectrum.h"

extern deck_constants my_const;
extern const mpi_info_struc mpi_info;

void spectrum::construct(){
/** \brief Generic spectrum contruction actions
*
*These are performed regardless of dimensions
*/
  function_type = 0;
  wave_id = 0;
  my_controller = nullptr;
  normB = 0;
  normg = nullptr;
  max_power=0.0;
  g_angle_array = nullptr;
  B_omega_array = nullptr;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));


}

spectrum::spectrum(int n_om, int n_ang, bool separable){
/** \brief Construct spectrum
*
*Constructs a spectrum. If the B and angle dependences are seperable (separable = true) then this will be B^2(omega) g(theta). Else it will be B^2(omega) g(omega, theta). Both have a B array of size n_omega, the former has g of 1 x n_angles, the latter n_omega x n_angles
*/

  construct();
  this->B_omega_array = new data_array(n_om);

  if(separable){
    this->g_angle_array = new data_array(1, n_ang);
    angle_is_function = true;
    function_type = FUNCTION_DELTA;
    normg = (my_type *) calloc(1, sizeof(my_type));
    //Single row so only one norm
  }else{
    this->g_angle_array = new data_array(n_om, n_ang);
    angle_is_function = false;
    normg = (my_type *) calloc(n_ang, sizeof(my_type));
    //Norm each row
  }
}

spectrum::spectrum(std::string filename){
/** \brief Setup a spectrum from file dump \todo test!!!*/

//First we grab the position of close block. Then we attempt to read two arrays. If we reach footer after first we error, or do not after second we warn.

  std::fstream file;
  file.open(filename, std::ios::in|std::ios::binary);
  if(!file.is_open()) return;

  bool err = 0, cont = 1;
  size_t end_block=0, next_block=0;
  file.seekg(-1*sizeof(size_t), file.end);
  file.read((char*) &end_block, sizeof(size_t));
  file.seekg(0, std::ios::beg);

  std::vector<size_t> dims = B_omega_array->read_dims_from_file(file);

  if(dims.size() !=1){
    my_print("Invalid dimensions for B", mpi_info.rank);
    cont = 0;
  }
  if(cont){
    //Now set up B correct size
    file.seekg(0, std::ios::beg);
    //return to start
    B_omega_array = new data_array(dims[0]);
    err = B_omega_array->read_from_file(file);
  }
  if(cont && err){
    my_print("File read failed", mpi_info.rank);
    cont = 0;
  }
  if(cont){
    file.seekg(-1*sizeof(size_t), std::ios::cur);
    file.read((char*) &next_block, sizeof(size_t));
    if(next_block == end_block){
      //Early termination of file...
      my_print("Insufficent arrays in file", mpi_info.rank);
      cont =0;
    }
  }
  if(cont){
    dims = g_angle_array->read_dims_from_file(file);
    if(dims.size() !=2){
      my_print("Wrong dimensions for g", mpi_info.rank);
      cont = 0;
    }
  }
  if(cont){
    //Now set up g correct size
    file.seekg(0, std::ios::beg);
    //return to start

    if(dims[0] == 1){
      this->g_angle_array = new data_array(1, dims[1]);
      angle_is_function = true;
      function_type = FUNCTION_DELTA;
      normg = (my_type *) calloc(1, sizeof(my_type));
      //Single row so only one norm
    }else{
      this->g_angle_array = new data_array(dims[0], dims[1]);
      angle_is_function = false;
      normg = (my_type *) calloc(dims[1], sizeof(my_type));
      //Norm each row
    }
  
    err = g_angle_array->read_from_file(file);
  
    if(err){
      my_print("File read failed", mpi_info.rank);
      cont =0;
    }
  }
  if(cont){
    file.seekg(-1*sizeof(size_t), std::ios::cur);
    file.read((char*) &next_block, sizeof(size_t));

    if(next_block != end_block){
      //File is not done!
      my_print("Excess arrays in file", mpi_info.rank);
      cont = 0;
    }
  }
  if(cont){
    char id_in[ID_SIZE];
    file.seekg(end_block+sizeof(size_t));
    if(file) file.read(id_in, sizeof(char)*ID_SIZE);
    strcpy(this->block_id, id_in);
  }

  if(!cont){
    //Something went wrong. Clean up and quit
    if(B_omega_array) delete B_omega_array;
    if(g_angle_array) delete g_angle_array;
    return;
  }
}


void spectrum::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[ID_SIZE], int function_type){
/**\brief Set parameters
*
*Sets the time and space ranges, wave type etc attached to the spectrum. Times should be in terms of file output time. Space in terms of grid points.
*/

  this->time[0] = time1;
  this->time[1] = time2;
  this->space[0] = space1;
  this->space[1] = space2;
  strcpy(this->block_id, block_id);
  this->wave_id = wave_id;
  this->function_type = function_type;
}

spectrum::~spectrum(){

  free(normg);
  if(g_angle_array) delete g_angle_array;
  if(B_omega_array) delete B_omega_array;

}

bool spectrum::generate_spectrum(data_array * parent, int om_fuzz, int angle_type){
/**\brief Generate spectrum from data
*
*Takes a parent data array and uses the specified ids to generate a spectrum. Windows using the specified wave dispersion and integrates over frequency. Also adopts axes from parent. \todo Ensure !angle_is_function forces all other rows to be equal length... \todo Fill in the rest of logic etc @param parent Data array to read from @param om_fuzz Band width around dispersion curve in percent of central frequency \todo omega vs k, is there some normalising to do? \todo Make parent const ref
*/

  if(!g_angle_array->is_good() || !B_omega_array->is_good()){
    my_print("Angle or omega arrays invalid. Returning", mpi_info.rank);
    return 1;
  }

  if(parent && angle_is_function){
    //First we read axes from parent
    this->copy_ids(parent);
    this->function_type = angle_type;
    size_t len;

    my_type *ax_ptr = parent->get_axis(1, len);
    //y-axis to work with
//    my_type *k_ax_ptr = parent->get_axis(0, lenk);

    //Now we loop across x, calculate the wave cutout bounds, and total, putting result into data

    int j;
    int low_bnd, high_bnd;
    my_type om_disp, max=0.0;
    my_type tolerance = om_fuzz/100.0;
    my_type total;
/*    my_type max_om = parent->get_axis_element(1, len-1);
    max_om /= this->get_length(0);*/
    //for(int i=0; i<this->get_length(0); ++i) this->set_axis_element(0, i, )

    for(size_t i=0; i<get_B_dims(0); ++i){
      om_disp = get_omega(parent->get_axis_element(0,i), WAVE_WHISTLER);
      
      set_om_axis_element(i, om_disp);
      
      low_bnd = where(ax_ptr, len, om_disp *(1.0-tolerance));
      high_bnd = where(ax_ptr, len, om_disp *(1.0+tolerance));
      if(low_bnd < 0 || high_bnd< 0){
        set_B_element(i,0.0);
        continue;
      }
      //now total the part of the array between these bnds
      total=0.0;
      for(j=low_bnd; j<high_bnd; j++) total += parent->get_element(i,j);
      set_B_element(i,total);
      if(total > max) max = total;
    }

    make_angle_distrib();
    this->max_power = total;

  }else if(parent){
 /** \todo general spectrum extracttion routine */
  //TODO in this case we have to extract spectrim and angle data somehow......

    //First we read axes from parent
    size_t len;

    my_type * ax_ptr = parent->get_axis(0, len);
    //memcpy ((void *)this->axes, (void *)ax_ptr, len*sizeof(my_type));
    ax_ptr = parent->get_axis(1, len);
    //y-axis to work with

    //Generate angle axis to work with
  {int res = 1;
  //set axis resolution somehow... TODO this
  g_angle_array->make_linear_axis(0, res, 0);
  }

  //and now we extract the data at each angle...
  //tan theta = k_x/k_y



  }else{
    return 1;

  }

  return 0;

}

bool spectrum::make_angle_distrib(){
/** \brief Generate angle axis and distribution
*
*Generates an angle axis linear spaced in tan theta between MIN_ANGLE and MAX_ANGLE. Then generates and fills the angular spectrum according to function specified by function_type member variable. Options are FUNCTION_DELTA: delta function with peak at 0 angle and integral 1. FUNCTION_GAUSS: Gaussian with std-dev SPECTRUM_ANG_STDDEV, centre at 0 angle and integral 1. FUNCTION_ISO: Isotropic distribution over range considered, with integral 1. FUNCTION \todo How to handle parallel/anti. \todo Integrals sum to 1 or 0.5???
*/


  if(!angle_is_function){
    my_print("Angular distrib is not a function. Returning", mpi_info.rank);
    return 1;
  }
  if(!g_angle_array->is_good()){
    my_print("Angular array invalid. Returning", mpi_info.rank);
    return 1;
  }
  
  calc_type res = (ANG_MAX - ANG_MIN)/get_g_dims(1);
  size_t len;
  int offset = -ANG_MIN/res;
  g_angle_array->make_linear_axis(1, res, offset);
  len = get_g_dims(1);
  my_type val;

  if(function_type == FUNCTION_DELTA){
  
    for(size_t i=0; i<len; ++i) set_g_element(i,0.0);
    val = 1.0/res;
    int zero = where(g_angle_array->get_axis(1, len), len, 0.0);
    //Set_element checks bnds automagically
    set_g_element(zero, val);

  }else if(function_type == FUNCTION_GAUSS){
    my_type ax_el;
    my_type norm;
    norm = 1.0/ (std::sqrt(2.0*pi) * SPECTRUM_ANG_STDDEV);
    for(size_t i=0; i<len; ++i){
      ax_el = get_ang_axis_element(i);
      val = std::exp( -0.5 * std::pow(ax_el/SPECTRUM_ANG_STDDEV, 2)) * norm;
      set_g_element(i,val);
    }
  }else if(function_type ==FUNCTION_ISO){

    val = 1.0/ (ANG_MAX - ANG_MIN)*get_g_dims(1)/(get_g_dims(1)-1);
    for(size_t i=0; i<len; ++i) set_g_element(i,val);

  }else{
  
    my_print("Invalid function type. Returning", mpi_info.rank);
    return 1;

  }
  
  return 0;

}

my_type spectrum::get_omega(my_type k, int wave_type, bool deriv){
/** \brief Gets omega for given k
*
* Calls to plasma because approximations for density etc etc should be made there. @param k Wavenumber @param wave_type Wave species @param deriv Return v_g instead
*/
  if(my_controller && (my_controller->get_plasma())) return (my_type) my_controller->get_plasma()->get_dispersion(k, wave_type, 0, deriv);
  else return 0.0;
}

my_type spectrum::get_k(my_type omega, int wave_type, bool deriv){
/** \brief Gets omega for given k
*
* Calls to plasma because approximations for density etc etc should be made there. @param k Wavenumber @param wave_type Wave species @param deriv Return v_g instead
*/
  if(my_controller && (my_controller->get_plasma())) return (my_type) my_controller->get_plasma()->get_dispersion(omega, wave_type, 1, deriv);
  else return 0.0;
}

std::vector<int> spectrum::all_where(my_type * ax_ptr, int len, my_type target,std::function<bool(my_type,my_type)> func){
/**\brief Find all points where
*
*Finds all indices where ax_ptr[i] vs target value satisfies the function given (e.g. std::greater)
*
*/
  std::vector<int> ret;

  for(int j=0;j<len; j++) if(func(ax_ptr[j],  target)) ret.push_back(j);

  return ret;

}

bool spectrum::write_to_file(std::fstream &file){
/** \brief Write to file \todo Swap to sentinel*/

  if(!file.is_open() || (!this->B_omega_array->is_good())|| (!this->g_angle_array->is_good())) return 1;

  bool err, write_err=0;
  err = B_omega_array->write_to_file(file, false);
  
  err |= g_angle_array->write_to_file(file, false);
  //Don't close these files. Instead the below writes a single footer

  size_t ftr_start = (size_t) file.tellg();
  //Start of ftr means where to start reading block, i.e. location of the next_location tag
  size_t next_location = ftr_start+ sizeof(char)*ID_SIZE;

  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  file.write(block_id, sizeof(char)*ID_SIZE);

  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  file.write((char*) & ftr_start, sizeof(size_t));

  return err;

}

bool spectrum::read_from_file(std::fstream &file){
/** \todo test*/

//First we grab the position of close block. Then we attempt to read two arrays. If we reach footer after first we error, or do not after second we warn. If the arrays in file are the wrong size, we have to error out
  bool err = 0;
  size_t end_block=0, next_block=0;
  file.seekg(-1*sizeof(size_t), file.end);
  file.read((char*) &end_block, sizeof(size_t));
  file.seekg(0, std::ios::beg);

  err = B_omega_array->read_from_file(file);
  if(err){
    my_print("File read failed", mpi_info.rank);
    return err;
  }
  file.seekg(-1*sizeof(size_t), std::ios::cur);
  file.read((char*) &next_block, sizeof(size_t));
  if(next_block == end_block){
    //Early termination of file...
    my_print("Insufficent arrays in file", mpi_info.rank);
    if(B_omega_array) delete B_omega_array;
    return 1;
  }

  err = g_angle_array->read_from_file(file);
  if(err){
    my_print("File read failed", mpi_info.rank);
    return err;
  }
  file.seekg(-1*sizeof(size_t), std::ios::cur);
  file.read((char*) &next_block, sizeof(size_t));

  if(next_block != end_block){
    //File is not done!
    my_print("Excess arrays in file", mpi_info.rank);
  }

  char id_in[ID_SIZE];
  file.seekg(end_block+sizeof(size_t));
  if(file) file.read(id_in, sizeof(char)*ID_SIZE);
  strcpy(this->block_id, id_in);

  return err;


}

void spectrum::make_test_spectrum(int time[2], int space[2],int angle_type){
/** \brief Generate dummy spectrum
*
*Makes a basic spectrum object with suitable number of points, and Gaussian(s) centred at fixed k/freq and x value \todo Should we use negative freqs?? If so use the offset as parameter @param time Time range (number of points) @param space Space range (number of points) @param angle_type Function to use for angular distrib
*/

  if(!this->B_omega_array->is_good() || !this->g_angle_array->is_good()){
    my_print("Array allocations failed, cannot make spectrum", mpi_info.rank);
    return;
  }
  char id[ID_SIZE] = "ex";

  this->set_ids(time[0], time[1], space[0], space[1], WAVE_WHISTLER, id, angle_type);
  
  //setup axes
  size_t len0;
  my_type * ax_ptr;

  ax_ptr = B_omega_array->get_axis(0, len0);
  my_type res;
  bool offset = true;
  //whether to have even ±pm axes or start from 0;
  res = 17000.0*1.0/(my_type)len0;
  offset= false;
  //res to cover range from offset to max in len0 steps

  //Rough value for length of
  for(size_t i=0; i<len0; i++) *(ax_ptr+i) = res*((my_type)i - (my_type)offset *(my_type)len0/2.);

  make_angle_distrib();

  //Generate the negative k data
  
  my_type centre, width, background = 0.0;
  centre = 14000, width=0.1*centre;
  
  my_type * data_tmp, *ax_tmp;
  my_type * data_ptr;
  data_ptr = (my_type *) malloc(len0*sizeof(my_type));
  data_tmp = data_ptr;
  ax_tmp = ax_ptr;

  if(offset){
    for(size_t i=0; i<=len0/2; i++, ax_tmp++, data_tmp++){
      *(data_tmp) = exp(-pow((*(ax_tmp) + centre), 2)/width/width) + background;
    }
    data_tmp--;
    //we've gone one past our termination condition...
    for(size_t i=1; i<len0/2; i++) *(data_tmp + i) = *(data_tmp - i);
  }else{
  
    for(size_t i=0; i<len0; i++, ax_tmp++, data_tmp++){
      *(data_tmp) = exp(-pow((*(ax_tmp) - centre), 2)/width/width) + background;
    }
  }
  this->B_omega_array->populate_data(data_ptr, len0);
  free(data_ptr);
  //reflect onto +ve axis
//+ 0.25*exp(-pow((*(ax_tmp) + centre), 2)/width*50.0)
  max_power = 1.0;
  //Store value of maximum
}

bool spectrum::truncate_om(my_type om_min, my_type om_max){
/** \brief Truncate omega distribution at om_min and om_max.
*
*Zeros all elements outside the range [abs(om_min), abs(om_max)]. Zeros are ignored. om_min must be < om_max.
*/

  if(om_min < 0) om_min = std::abs(om_min);
  if(om_max < 0) om_max = std::abs(om_max);
  if(om_min >= om_max){
    my_print("Invalid omega range, aborting truncate", mpi_info.rank);
    return 1;
  
  }

  int index = -1;
  int len = get_B_dims(0);

  if(om_min != 0.0){
    index=where_omega(om_min);
    if(index != -1) for(int i=0; i< index; i++) set_B_element(i, 0.0);
    //Zero up to om_min
  }

  if(om_max != 0.0 && om_max < get_om_axis_element(len-1)){
    index=where_omega(om_max);
    if(index != -1) for(int i = index; i< len; i++) set_B_element(i, 0.0);
    //Zero after to om_max
  }

  normaliseB();
  //Re-do normalisation
  return 0;

}

bool spectrum::truncate_x(my_type x_min, my_type x_max){
/** \brief Truncate angle distribution at x_min and x_max.
*
*Zeros all elements outside the range [abs(x_min), abs(x_max)]. Zeros are ignored. x_min must be < x_max.
*/

  if(x_min >= x_max){
    my_print("Invalid x range, aborting truncate", mpi_info.rank);
    return 1;
  }

  int index = -1;
  size_t len = get_g_dims(1), om_len = get_g_dims(0);

  if(x_min > ANG_MIN){
    index = where(g_angle_array->get_axis(1, len), len, x_min);
    if(index != -1){
      for(size_t i=0; i< (size_t)index; i++){
        for(size_t j=0; j< om_len; j++){
          set_g_element(i, j, 0.0);
        }
      }
    }
  
  }
  if(x_max < ANG_MAX){
    index = where(g_angle_array->get_axis(0, len), len, x_max);
    if(index != -1){
      for(size_t i=index; i< len; i++){
        for(size_t j=0; j< om_len; j++){
          set_g_element(i, j, 0.0);
        }
      }
    }
  }

  return 0;

}

int spectrum::where_omega(my_type omega){
/** \brief Gets where omega axis exceeds passed value
*
*Finds where frequency or wavenumber axis exceeds passed omega, using dispersion relation to transform k to omega if necessary.
*/
  int index;
  size_t len;
  B_omega_array->get_axis(0, len);

  index = where(B_omega_array->get_axis(0, len), len, omega);
  return index;

}

bool spectrum::normaliseB(){
/** \brief Normalise B(w)
*
* Calculate the total square integral of values over range \todo Is this data bare or squared?
*/

  int len = get_B_dims(0);
  my_type * d_axis = (my_type *) calloc(len, sizeof(my_type));
  my_type * data = (my_type *) malloc(len*sizeof(my_type));

  for(int i=0; i<len-1; i++) d_axis[i] = get_om_axis_element(i+1) - get_om_axis_element(i);
  for(int i=0; i<len; i++) data[i] = get_B_element(i);
  
  normB = integrator(data, len, d_axis);
  free(d_axis);
  free(data);
  return 0;
}

bool spectrum::normaliseg(my_type omega){
/** \brief Normalise g_w(x)
*
*Calculate the norm of g used in e.g. denom of Albert eq 3 or calc'd in derivations.tex. We assume omega, x are off the axes already so no interpolation  \todo Catch zero norms \todo CHECK
*/

  size_t len = get_g_dims(0);
  plasma * plas =my_controller->get_plasma();

  my_type * d_axis = (my_type *) calloc(len, sizeof(my_type));
  my_type * integrand = (my_type *) calloc(len, sizeof(my_type));

  for(size_t i=0; i<len-1; i++) d_axis[i] = get_om_axis_element(i+1) - get_om_axis_element(i);
  //Construct dx axis for integration /** \todo omega or ang???*/

  mu my_mu;

  int om_ind = 1;
  //skip over B data
  
  size_t lena = get_B_dims(0);
  if(!angle_is_function){
    om_ind = where(B_omega_array->get_axis(0, lena), lena, omega);
  }
  if(om_ind<0) return 1;
  //break if Omega is out of range

  my_type x, psi;
  size_t inda, indb;

  //Addressing changes if we have g(x) or g(w, x)
  angle_is_function ? indb=om_ind: inda=om_ind;
  /** \todo Not anymore!!! */
  for(size_t i=0; i<len; i++){
    x = get_ang_axis_element(i);
    psi = atan(x);
    my_mu = plas->get_root(0.0, omega, psi);

    angle_is_function ? inda=i: indb=i;

    if(!my_mu.err){
      integrand[i] = get_g_element(inda, indb) * x * std::pow((std::pow(x, 2)+1.0), -1.5)*std::pow(my_mu.mu, 2) * std::abs( my_mu.mu + omega*my_mu.dmudom);
    }
    //product of g(theta) * x (x^2+1)^-(3/2) * mu^2 |mu+omega dmu/domega|
  }
  
  //integrate
  my_type normg_tmp = integrator(integrand, len, d_axis);

  if(angle_is_function) om_ind -=1;
  normg[om_ind] = normg_tmp;
  //store into right place
  
  //clean up
  free(d_axis);
  free(integrand);
  return 0;
}

calc_type spectrum::get_G1(calc_type omega){
/** \brief G1 from Albert 2005.
*
*Gets the value of B(w) (interpolated if necessary) and the normalising constant from normB
\todo Does it matter that our k is limited? Do waves really go to low intensity in bit we see \todo Do we need the vg conversion factor?
*/

  calc_type B2;
  my_type tmpB2;
  if(normB ==0.0) normaliseB();
  size_t len, offset;
  my_type data_bit[2];
  my_type ax_val;
  my_type * axis = B_omega_array->get_axis(0, len);

  //If we have B(k) we need to change to B(w)
  my_type change_of_vars = 1.0;

  ax_val = (my_type) omega;
  
  offset = where(axis, len, ax_val);
  //Interpolate if possible, else use the end
  if(offset > 0 && offset < len){
    data_bit[0] = get_B_element(offset-1);
    data_bit[1] = get_B_element(offset);
    tmpB2 = interpolate(B_omega_array->get_axis(0, len) + offset-1, data_bit, ax_val, 2);
    //tmpB2 = data_bit[0];
  }else if(offset == 0){
    //we're right at end, can't meaningfully interpolate, use raw
    tmpB2 = get_B_element(0);
  }else{
    //offset <0 or > len, value not found
    tmpB2 = 0.0;
  }

  //Add change_of_vars constant in case we have k axis
  B2 = (calc_type) tmpB2 * change_of_vars;

  //Add norm. constant
  return B2/normB;

}

calc_type spectrum::get_G2(calc_type omega, calc_type x){
/** \brief Get G2 from Albert 2005
*
* Gets the value of g(w, x) and the normalising constant from normg \todo IS THIS OMEGA OR do we calc omega according to conditions on integral??? \todo interpolate on omega? or angle or both. Or fix angle axis as matched to D. In some sense we want to minimise work here...
*/


  int om_ind, offset;
  size_t len;
  my_type tmpg;
  my_type data_bit[2];

  len = get_B_dims(0);

  if(!angle_is_function){
    om_ind = where(B_omega_array->get_axis(0, len), len, omega);
  }
  else om_ind = 0;
  if(om_ind>=0 && normg[om_ind] == 0.0){
    normaliseg(omega);
  }
 // std::cout<< normg[0]<<" "<<std::endl;
  
  //Bump up to miss B row
  my_type * axis = g_angle_array->get_axis(0, len);
  offset = where(axis, len, x);
  
  //Interpolate if possible, else use the end
  if(offset > 0 && offset < len){
    data_bit[0] = get_g_element(offset-1, om_ind);
    data_bit[1] = get_g_element(offset, om_ind);
    tmpg = interpolate(axis + offset-1, data_bit, (my_type)x, 2);

  }else if(offset==0){
    //we're right at end, can't meaningfully interpolate, use raw
    tmpg = get_g_element(offset, om_ind);

  }else{
    //offset <0 or > len, value not found
    tmpg = 0.0;
  }

  
  if(offset >=0 && om_ind >=0)return tmpg/normg[om_ind];
  else return 0.0;

}

calc_type spectrum::check_upper(){
/** \brief Check upper k limit of spectral power
*
* Checks the upper bound of region of significant spectral power, i.e. above SPECTRUM_THRESHOLD*peak_power \todo This wont work with +ve only spectrum
*/

//First move up from bottom in strides and find where _first_ rises above threshold. Specifcally upwards.

  size_t stride = 1;
  size_t ax_len = get_B_dims(0);
  size_t len = ax_len/2 - stride;
  my_type threshold = SPECTRUM_THRESHOLD*this->max_power, k_thresh;
  size_t index=0;

  //Naive check for symmetric or +ve only. Allow up to 2 d ax mismatch for odd/even lengths
  my_type d_axis = std::abs(get_om_axis_element(1) -get_om_axis_element(2));

  if(std::abs(get_om_axis_element(1) + get_om_axis_element(ax_len-1)) > 2.0 * d_axis){

  //Axis presumed to be +ve only, move in from top
    for(size_t i=0; i< ax_len; i+=stride){
      if((get_B_element(ax_len - i) < threshold && get_B_element(ax_len - i - stride) > threshold)){
        index = i;
        break;
      }
    }
  
  
  }else{
  //Axis presumed to be symmetrical

    for(size_t i=0; i< len; i+=stride){
      if((get_B_element(i) < threshold && get_B_element(i+stride) > threshold) ||(get_B_element(ax_len - i) < threshold && get_B_element(ax_len - i - stride) > threshold)){
        index = i;
        break;
      }
    }
  }
  //get omega value and convert to k...
  if(index != 0){
    my_type tmp = std::abs(get_om_axis_element(index));
    k_thresh = get_k(tmp, WAVE_WHISTLER);
  }else{
  //there either isn't a peak, or isn't waves or whatever. So we pick something arbitrary...
    my_type tmp = my_const.omega_ce * 0.999;
    k_thresh = get_k(tmp, WAVE_WHISTLER);
  
  }
  return k_thresh;
}

calc_type spectrum::get_peak_omega(){
/** Find position of spectral peak
*
*Finds location of highest peak in spectrum.
*/

  calc_type value = -1.0, tmp;
  int index;
  for(size_t i=0; i<get_B_dims(0); ++i){
    tmp = get_B_element(i);
    if(tmp > value){
      index = i;
      value = tmp;
    }

  }

  return get_om_axis_element(index);

}

void spectrum::copy_ids( data_array * src){
/** Copies ID fields from src array to this */

  strcpy(this->block_id, src->block_id);
  
  std::copy(src->time, src->time + 2, this->time);
  for(int i=0; i < 2; ++i) this->space[i] = src->space[i];
  //if(angle)
}

bool spectrum::check_ids( data_array * src){
/** Checks ID fields match src */

  bool err=false;
  if(strcmp(this->block_id, src->block_id) != 0) err =true;
  for(int i=0; i< 3; i++) if(src->time[i] != this->time[i]) err=true;
  for(int i=0; i < 2; ++i) if(this->space[i] != src->space[i]) err=true;

  return err;
}


