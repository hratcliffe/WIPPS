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
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));
  smooth=0;
}

spectrum::spectrum(){
  construct();
}
spectrum::spectrum(int n_om, int n_ang, bool separable){
/** \brief Construct spectrum
*
*Constructs a spectrum. If the B and angle dependences are seperable (separable = true) then this will be B^2(omega) g(theta). Else it will be B^2(omega) g(omega, theta). Both have a B array of size n_omega, the former has g of 1 x n_angles, the latter n_omega x n_angles
*/

  construct();
  this->B_omega_array = data_array(n_om);

  if(separable){
    this->g_angle_array = data_array(1, n_ang);
    angle_is_function = true;
    function_type = FUNCTION_DELTA;
    normg = (my_type *) calloc(1, sizeof(my_type));
    //Single row so only one norm
  }else{
    this->g_angle_array = data_array(n_om, n_ang);
    angle_is_function = false;
    normg = (my_type *) calloc(n_ang, sizeof(my_type));
    //Norm each row
  }
}

spectrum::spectrum(std::string filename){
/** \brief Setup a spectrum from file dump*/

//First we grab the position of close block. Then we attempt to read two arrays. If we reach footer after first we error, or do not after second we warn.

  construct();
  std::fstream file;
  file.open(filename, std::ios::in|std::ios::binary);
  if(!file.is_open()) return;

  bool err = 0, cont = 1;
  size_t end_block=0, next_block=0;
  size_t jump_pos=0;
  file.seekg(-1*sizeof(size_t), file.end);
  file.read((char*) &end_block, sizeof(size_t));
  file.seekg(0, std::ios::beg);

  std::vector<size_t> dims = B_omega_array.read_dims_from_file(file);

  if(dims.size() !=1){
    my_print("Invalid dimensions for B", mpi_info.rank);
    cont = 0;
  }else{
    my_print("B size "+mk_str(dims[0]), mpi_info.rank);
  
  }
  if(cont){
    //Now set up B correct size
    file.seekg(0, std::ios::beg);
    //return to start
    B_omega_array = data_array(dims[0]);
    err = B_omega_array.read_from_file(file);
  }
  if(cont && err){
    my_print("File read failed", mpi_info.rank);
    cont = 0;
  }
  if(cont){
    file.read((char*) &next_block, sizeof(size_t));
    if(next_block == end_block){
      //Early termination of file...
      my_print("Insufficent arrays in file", mpi_info.rank);
      cont =0;
    }
  }
  if(cont){
    file.seekg(-1*sizeof(size_t), std::ios::cur);
    //Last read was invalid, so we seek back

    jump_pos = (size_t)file.tellg();
    dims = g_angle_array.read_dims_from_file(file);

    if(dims.size() !=2){
      my_print("Wrong dimensions for g", mpi_info.rank);
      cont = 0;
    }else{
      my_print("g size "+mk_str(dims[0])+'x'+mk_str(dims[1]), mpi_info.rank);
    }
  
  }
  if(cont){
    file.seekg(jump_pos);
    //Now set up g correct size
    if(dims[0] == 1){
      this->g_angle_array = data_array(1, dims[1]);
      angle_is_function = true;
      function_type = FUNCTION_DELTA;
      normg = (my_type *) calloc(1, sizeof(my_type));
      //Single row so only one norm
    }else{
      this->g_angle_array = data_array(dims[0], dims[1]);
      angle_is_function = false;
      normg = (my_type *) calloc(dims[1], sizeof(my_type));
      //Norm each row
    }
  
    err = g_angle_array.read_from_file(file);
  
    if(err){
      my_print("File read failed", mpi_info.rank);
      cont =0;
    }
  }
  if(cont){
    file.read((char*) &next_block, sizeof(size_t));
    next_block= (size_t)file.tellg() - sizeof(size_t);
    //If we're done, this should be the
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
    g_angle_array = data_array();
    B_omega_array = data_array();
    //Return to size-less state
    
  }
  return;

}

spectrum & spectrum::operator=(const spectrum& src){
  
  if(this->normg) free(normg);
  construct();
  if(!is_good()) return *this;
  //Stop if src not good

  size_t g_sz = src.get_g_dims(0);
  normg = (my_type *) calloc(g_sz, sizeof(my_type));
  std::copy(src.normg, src.normg + g_sz, this->normg);
  return *this;
}

spectrum::spectrum(const spectrum &src){
/** \brief Copy constructor
*
*Copy src to a new instance. Requires we copy the normg block
*/
  construct();
  if(!is_good()) return;
  //Stop if src not good

  size_t g_sz = src.get_g_dims(0);
  normg = (my_type *) calloc(g_sz, sizeof(my_type));
  std::copy(src.normg, src.normg + g_sz, this->normg);
}

spectrum::~spectrum(){
  if(normg) free(normg);
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

void spectrum::init(){
/** \brief Cache spectrum values
*
*Fill cached values for max power, norms etc*/

  max_power = B_omega_array.maxval();
  normaliseB();
  size_t n = get_B_dims();
  my_type omega;
  for(size_t i=0; i< n; i++){
    omega = get_om_axis_element(i);
    normaliseg(omega);
  }
}

bool spectrum::generate_spectrum(data_array &parent, int om_fuzz, int angle_type, data_array * mask){
/**\brief Generate spectrum from data
*
*Takes a parent data array and uses the specified ids to generate a spectrum. Windows using the specified wave dispersion and integrates over frequency. Also adopts axes from parent. IMPORTANT: when using real angular data we roughly fuzz around the correct k values, but this is not uniform! Non-smooth or rapidly varying spectra may give odd results @param parent Data array to read from @param om_fuzz Band width around dispersion curve in percent of central frequency \todo omega vs k, is there some normalising to do? 
*/

  if(!this->is_good()){
    my_print("Spectrum object invalid. Returning", mpi_info.rank);
    return 1;
  }

  if(parent.is_good() && angle_is_function){
    //First we read axes from parent
    this->copy_ids(parent);
    this->function_type = angle_type;
    size_t len;

    my_type *ax_ptr = parent.get_axis(1, len);
    //y-axis to work with

    //Now we loop across x, calculate the wave cutout bounds, and total, putting result into data

    int j;
    int low_bnd, high_bnd;
    my_type om_disp, max=0.0;
    my_type tolerance = om_fuzz/100.0;
    my_type total = 0.0;

    calc_type sgn;
    bool one_sided = (std::abs(ax_ptr[0]) < std::abs(ax_ptr[len/2]));

    for(size_t i=0; i<get_B_dims(0); ++i){
      om_disp = get_omega(parent.get_axis_element(0,i), WAVE_WHISTLER);
      if(!one_sided && i<= get_B_dims(0)/2) sgn = -1.0;
      else sgn = 1.0;
      
      om_disp *= sgn;
      set_om_axis_element(i, om_disp);
      
      low_bnd = where(ax_ptr, len, om_disp *(1.0-sgn*tolerance));
      high_bnd = where(ax_ptr, len, om_disp *(1.0+sgn*tolerance));
      
      //Adjust for low off bottom or high off top
      //Low off top or high off bottom will still be caught and skipped
      if(low_bnd < 0 && om_disp *(1.0-sgn*tolerance) < ax_ptr[0]) low_bnd = 0;
      if(high_bnd < 0 && om_disp *(1.0+sgn*tolerance) > ax_ptr[len-1]) high_bnd = len-1;

      if(low_bnd < 0 || high_bnd< 0){
        set_B_element(i,0.0);
        continue;
      }
      //now total the part of the array between these bnds
      total=0.0;
      for(j=low_bnd; j<high_bnd; j++) total += parent.get_element(i,j);
      if(mask) for(j=low_bnd; j<high_bnd; j++) mask->set_element(i, j, mask->get_element(i, j)+1);
      
      set_B_element(i,total);
      if(total > max) max = total;
    }

    make_angle_distrib();
    this->max_power = max;

  }else if(parent.is_good() && parent.get_dims()==3){
  //This will do a rough integral around the correct omega. Note it is not a neat sampling of the domain, so may give odd results for non-smooth cases
  
    this->copy_ids(parent);

    my_type tolerance = om_fuzz/100.0;

    //Get the k_x and k_y axes for lookup
    size_t len_x, len_y, len_om;
    my_type * kx_ax = parent.get_axis(0, len_x);
    my_type * ky_ax = parent.get_axis(1, len_y);
    my_type * om_ax = parent.get_axis(2, len_om);
    
    bool one_sided = (std::abs(om_ax[0]) < std::abs(om_ax[len_om/2]));
    my_type om_disp;

    make_angle_axis();
    
    my_type k, k_y, tantheta, tmp, tmp_sum, decrement = 1.0+GEN_PRECISION, max_power=0.0, theta;
    int kx_low, kx_high, ky_low, ky_high, i_sgn=1, om_ind, om_low, om_high;
    //Now we do a double loop
    
    for(size_t i = 0; i< len_x; i++){
      //Signs so we always integrate from k_low to k_high, if one sided is always 1
      if(!one_sided) i_sgn = (i<len_x/2 ? -1: 1);
      if(i_sgn > 0) decrement = 1.0 - GEN_PRECISION;

      //use the k_x axis to generate an omega one. Thus if fuzz=0 we match 1-1
      om_disp = get_omega(parent.get_axis_element(0,i), WAVE_WHISTLER);
      set_om_axis_element(i, om_disp*i_sgn);

      //We can assume either +ve and -ve k (include i_sgn in kx_high, or omega (include in om_low and high)
      //Omega in parent data matching current value
      om_ind = where(om_ax, len_om, om_disp);
      
      tmp_sum = 0.0;
      for(size_t j = 0; j< get_g_dims(1); j++){
        tantheta = get_ang_axis_element(j);
        theta = atan(tantheta);
//----------Fuzz omega--------------------
        om_low = where(om_ax, len_om, om_disp*(1.0-i_sgn*tolerance)*i_sgn);
        om_high = where(om_ax, len_om, om_disp*(1.0+i_sgn*tolerance)*i_sgn);
//        om_low = where(om_ax, len_om, om_disp*(1.0-tolerance));
//        om_high = where(om_ax, len_om, om_disp*(1.0+tolerance));
        if(om_low < 0) om_low = 0;
        if(om_disp*(1.0+tolerance) < *om_ax) om_high = 0;
        else if(om_high < 0) om_high = len_om - 1;


//-------No fuzz in k_x -----------------------
        k = get_k(om_disp, WAVE_WHISTLER, 0, theta);//*i_sgn;//*decrement;
        
        kx_high = where(kx_ax, len_x, k*cos(theta));
        kx_low = kx_high;
        
        if(k*cos(theta)< *(kx_ax) || std::abs(k)<GEN_PRECISION || k*cos(theta)> *(kx_ax+len_x -1)) kx_high = -1;
        else if(kx_high <0) kx_high = len_x-1;
        if(kx_low < 0) kx_low = 0;
        
//-------- k_y is fuzzy----------
        k = get_k(om_disp*(1.0+i_sgn*tolerance), WAVE_WHISTLER, 0, theta)*i_sgn;//*decrement;
        ky_high = where(ky_ax, len_y, k*sin(theta));
        if(std::abs(k)<GEN_PRECISION || std::abs(k*sin(theta))> *(ky_ax+len_y -1)) ky_high = -1;
        else if(ky_high <0) ky_high = len_y-1;


        k = get_k(om_disp*(1.0-i_sgn*tolerance), WAVE_WHISTLER, 0, theta)*i_sgn;//*decrement;
        ky_low = where(ky_ax, len_y, k*sin(theta));
        if(ky_low < 0) ky_low = 0;
        
//-------Total on small range and create mask if requested-----
        tmp = 0.0;
        if(kx_high != -1 && ky_high != -1){
          for(int ii=om_low; ii<=om_high; ii++){
            for(int jj=ky_low; jj<=ky_high; jj++){
              tmp += parent.get_element(kx_high, jj, ii);
            }
          }
        }
        if(mask){
          if(kx_high != -1 && ky_high != -1){
            for(int ii=om_low; ii<=om_high; ii++){
              for(int jj=ky_low; jj<=ky_high; jj++){
                mask->set_element(kx_high, jj, ii, mask->get_element(kx_high, jj, ii)+1);
              }
            }
          }
       }
//------------------
        //Tmp is now the convolved B(omega)g(omega, theta). We dump it into g as is for later norming and add to the sum which will end up as B.
        /** \todo Is this the correct norm?*/
        set_g_element(i, j, tmp);
        tmp_sum += tmp;
      }
      //Now the angles are done, we norm. g and set B
      //Use some small value to prevent a divide by zero
      if(std::abs(tmp_sum) > tiny_my_type){
        for(size_t j = 0; j< get_g_dims(1); j++) set_g_element(i, j, get_g_element(i, j)/tmp_sum);
      }
      set_B_element(i, tmp_sum);
      //Store the max power for later
      if(tmp_sum > max_power) max_power = tmp_sum;

    }
    this->max_power = max_power;
  }else{
    return 1;
  }

  return 0;

}

void spectrum::make_angle_axis(){
/** \brief Generate angle axis
*
* Uniform angle axis from ANG_MIN to ANG_MAX
*/
  calc_type res = (ANG_MAX - ANG_MIN)/get_g_dims(1);
  size_t len;
  int offset = -ANG_MIN/res;
  g_angle_array.make_linear_axis(1, res, offset);

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
  if(!g_angle_array.is_good()){
    my_print("Angular array invalid. Returning", mpi_info.rank);
    return 1;
  }
  make_angle_axis();
  size_t len = get_g_dims(1);
  my_type val;

  if(function_type == FUNCTION_DELTA){
  
    for(size_t i=0; i<len; ++i) set_g_element(i,0.0);
    val = 1.0/(get_ang_axis_element(1)-get_ang_axis_element(0));
    int zero = where(g_angle_array.get_axis(1, len), len, 0.0);
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

void spectrum::smooth_B(int n_pts){
/** \brief Smooth B
*
*Apply a box-car smoothing to B with specified number of pts. Store n_pts
*/

  this->smooth = n_pts;
  B_omega_array.smooth_1d(n_pts);
}

my_type spectrum::get_omega(my_type k, int wave_type, bool deriv, my_type theta){
/** \brief Gets omega for given k
*
* Calls to plasma because approximations for density etc etc should be made there. @param k Wavenumber @param wave_type Wave species @param deriv Return v_g instead
*/
  if(my_controller && (my_controller->is_good())) return (my_type) my_controller->get_plasma().get_dispersion(k, wave_type, 0, deriv, theta);
  else return 0.0;
}

my_type spectrum::get_k(my_type omega, int wave_type, bool deriv, my_type theta) {
/** \brief Gets omega for given k
*
* Calls to plasma because approximations for density etc etc should be made there. @param k Wavenumber @param wave_type Wave species @param deriv Return v_g instead
*/
  if(my_controller && my_controller->is_good()) return (my_type) my_controller->get_plasma().get_dispersion(omega, wave_type, 1, deriv, theta);
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
/** \brief Write to file*/

  if(!file.is_open() || (!this->B_omega_array.is_good())|| (!this->g_angle_array.is_good())) return 1;

  bool err, write_err=0;
  err = B_omega_array.write_to_file(file, false);
  
  err |= g_angle_array.write_to_file(file, false);
  //Don't close these files. Instead the below writes a single footer

  size_t ftr_start = (size_t) file.tellg();
  //Start of ftr means where to start reading block, i.e. location of the next_location tag
  size_t next_location = ftr_start+ sizeof(char)*ID_SIZE +sizeof(size_t)*2;

  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section, i.e. of file end val
  file.write(block_id, sizeof(char)*ID_SIZE);
  
  file.write((char*) & smooth, sizeof(size_t));

  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  file.write((char*) & ftr_start, sizeof(size_t));

  return err;

}

bool spectrum::read_from_file(std::fstream &file){
/** \brief Initialise spectrum from dump
*
*Reads a dump file which is expected to contain two arrays, first B then g, as written by spectrum->write_to_file and constructs spectrum from data
*/

//First we grab the position of close block. Then we attempt to read two arrays. If we reach footer after first we error, or do not after second we warn. If the arrays in file are the wrong size, we have to error out
  bool err = 0;
  size_t end_block=0, next_block=0;
  file.seekg(-1*sizeof(size_t), file.end);
  file.read((char*) &end_block, sizeof(size_t));
  file.seekg(0, std::ios::beg);

  err = B_omega_array.read_from_file(file);
  if(err){
    my_print("File read failed", mpi_info.rank);
    return err;
  }
  file.seekg(-1*sizeof(size_t), std::ios::cur);
  file.read((char*) &next_block, sizeof(size_t));
  if(next_block == end_block){
    //Early termination of file...
    my_print("Insufficent arrays in file", mpi_info.rank);
    B_omega_array = data_array();
    return 1;
  }

  err = g_angle_array.read_from_file(file);
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

void spectrum::make_test_spectrum(int angle_type, bool two_sided, my_type om_ce){
/** \brief Generate dummy spectrum
*
*Makes a basic spectrum object with suitable number of points, and Gaussian(s) centred at fixed k/freq and x value @param time Time range (number of points) @param space Space range (number of points) @param angle_type Function to use for angular distrib @param two_sided Whether to generate symmetric spectrum or one-sided
*/

  if(!this->is_good()){
    my_print("Spectrum object invalid, returning", mpi_info.rank);
    return;
  }
  char id[ID_SIZE] = "ex";

  this->set_ids(0, 1, 0, 1, WAVE_WHISTLER, id, angle_type);
  
  //setup axes
  size_t len0;
  my_type * ax_ptr;

  ax_ptr = B_omega_array.get_axis(0, len0);
  my_type res = om_ce*(1.0+(my_type)two_sided)/(my_type)len0;
  //res to cover range from either 0 or -max to max in len0 steps

  //Generate axes for one or two-sided
  for(size_t i=0; i<len0; i++) *(ax_ptr+i) = res*((my_type)i - (my_type)two_sided *(my_type)len0/2.);

  make_angle_distrib();
  
  my_type centre, width, background = 0.0;
  centre = om_ce*14.0/17.0, width=0.1*centre;
  //Matches our test spectrum data
  
  my_type * data_tmp, *ax_tmp;
  my_type * data_ptr;
  data_ptr = (my_type *) malloc(len0*sizeof(my_type));

  //Initial values
  data_tmp = data_ptr;
  ax_tmp = ax_ptr;

  if(two_sided){
    //Fill a two-sided spectrum with peaks where specified
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
  this->B_omega_array.populate_data(data_ptr, len0);
  free(data_ptr);

  max_power = 1.0;
  //Store value of maximum
}

bool spectrum::truncate_om(my_type om_min, my_type om_max){
/** \brief Truncate omega distribution at om_min and om_max.
*
*Zeros all elements outside the range [om_min, om_max]. Zeros are ignored. om_min must be < om_max. Om_min or max out of axis range does nothing on that end. NB B is renormalised after the truncation
*/

  if(om_min >= om_max){
    my_print("Invalid omega range, aborting truncate", mpi_info.rank);
    return 1;
  }

  int index = -1;
  int len = get_B_dims(0);

  if(om_min > get_om_axis_element(1)){
    index=where_omega(om_min);
    if(index != -1) for(int i=0; i< index; i++) set_B_element(i, 0.0);
    //Zero up to om_min
  }

  if(om_max < get_om_axis_element(len-1)){
    index=where_omega(om_max);
    if(index != -1) for(int i = index; i< len; i++) set_B_element(i, 0.0);
    //Zero after om_max
  }

  normaliseB();
  //Re-do normalisation
  return 0;
}

bool spectrum::truncate_x(my_type x_min, my_type x_max){
/** \brief Truncate angle distribution at x_min and x_max.
*
*Zeros all elements outside the range [x_min, x_max]. x_min must be < x_max. If x_min or max are out of range, nothing is done at that end. \todo Should we renorm. g now?
*/

  if(x_min >= x_max){
    my_print("Invalid x range, aborting truncate", mpi_info.rank);
    return 1;
  }

  int index = -1;
  size_t len = get_g_dims(1), om_len = get_g_dims(0);

  if(x_min > ANG_MIN){
    index = where(g_angle_array.get_axis(1, len), len, x_min);
    if(index != -1){
      for(size_t i=0; i< (size_t)index; i++){
        for(size_t j=0; j< om_len; j++){
          set_g_element(i, j, 0.0);
        }
      }
    }
  
  }
  if(x_max < ANG_MAX){
    index = where(g_angle_array.get_axis(0, len), len, x_max);
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
  B_omega_array.get_axis(0, len);

  index = where(B_omega_array.get_axis(0, len), len, omega);
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
*Calculate the norm of g used in e.g. denom of Albert eq 3 or calc'd in derivations.tex. We assume omega, x are off the axes already so no interpolation.  \todo Catch zero norms \todo test \todo Taking omega wont work how we want
*/

  size_t len = get_g_dims(1);
  if(!my_controller) return 1;
  plasma plas = my_controller->get_plasma();

  my_type * d_axis = (my_type *) calloc(len, sizeof(my_type));
  my_type * integrand = (my_type *) calloc(len, sizeof(my_type));

  for(size_t i=0; i<len-1; i++) d_axis[i] = get_ang_axis_element(i+1) - get_ang_axis_element(i);
  //Construct dx axis for integration

  mu my_mu;

  long om_ind = 0;
  
  size_t lena = get_B_dims(0);
  if(!angle_is_function){
    om_ind = where(B_omega_array.get_axis(0, lena), lena, omega);
  }
  if(om_ind<0) return 1;
  //break if Omega is out of range

  my_type x, psi;

  for(size_t i=0; i<len; i++){
    x = get_ang_axis_element(i);
    psi = atan(x);
    my_mu = plas.get_root(0.0, omega, psi);

    if(!my_mu.err){
      integrand[i] = get_g_element(om_ind, i) * x * std::pow((std::pow(x, 2)+1.0), -1.5)*std::pow(my_mu.mu, 2) * std::abs( my_mu.mu + omega*my_mu.dmudom);
    }
    //product of g(theta) * x (x^2+1)^-(3/2) * mu^2 |mu+omega dmu/domega|
  }
  
  //integrate
  my_type normg_tmp = integrator(integrand, len, d_axis);

  normg[om_ind] = normg_tmp;
  
  //clean up
  free(d_axis);
  free(integrand);
  return 0;
}

calc_type spectrum::get_G1(calc_type omega){
/** \brief G1 from Albert 2005.
*
*Gets the value of B(w) (interpolated if necessary) and the normalising constant from normB
\todo Does it matter that our k is limited? Do waves really go to low intensity in bit we see \todo Do we need the vg conversion factor? \todo CHECK and FIXXXX and test
*/

  calc_type B2;
  my_type tmpB2;
  if(normB ==0.0) normaliseB();
  size_t len, offset;
  my_type data_bit[2];
  my_type ax_val;
  my_type * axis = B_omega_array.get_axis(0, len);

  //If we have B(k) we need to change to B(w)
  my_type change_of_vars = 1.0;

  ax_val = (my_type) omega;
  
  offset = where(axis, len, ax_val);
  //Interpolate if possible, else use the end
  if(offset > 0 && offset < len){
    data_bit[0] = get_B_element(offset-1);
    data_bit[1] = get_B_element(offset);
    tmpB2 = interpolate(B_omega_array.get_axis(0, len) + offset-1, data_bit, ax_val, 2);
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
* Gets the value of g(w, x) and the normalising constant from normg \todo IS THIS OMEGA OR do we calc omega according to conditions on integral??? \todo interpolate on omega? or angle or both. Or fix angle axis as matched to D. In some sense we want to minimise work here... \todo CHECK and FIXXXX and test
*/


  int om_ind, offset;
  size_t len;
  my_type tmpg;
  my_type data_bit[2];

  len = get_B_dims(0);

  if(!angle_is_function){
    om_ind = where(B_omega_array.get_axis(0, len), len, omega);
  }
  else om_ind = 0;
  
  if(om_ind>=0 && normg[om_ind] == 0.0){
    normaliseg(omega);
  }
 // std::cout<< normg[0]<<" "<<std::endl;
  
  //Bump up to miss B row
  my_type * axis = g_angle_array.get_axis(1, len);
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
* Checks the upper bound of region of significant spectral power, i.e. above SPECTRUM_THRESHOLD*peak_power
*/

  size_t stride = 1;
  size_t ax_len = get_B_dims(0);
  size_t len = ax_len/2 - stride;
  my_type threshold = SPECTRUM_THRESHOLD*this->max_power;
  size_t index = 0;

  if(get_om_axis_element(1) < get_om_axis_element(2)){
  //Axis presumed to be +ve only, move in from top to find first rise above threshold
    for(size_t i=0; i< ax_len; i+=stride){
      if((get_B_element(ax_len - i) < threshold && get_B_element(ax_len - i - stride) > threshold)){
        index = i;
        break;
      }
    }
  }else{
  //Axis presumed to be symmetrical. Move in from both ends
    for(size_t i=0; i< len; i+=stride){
      if((get_B_element(i) < threshold && get_B_element(i+stride) > threshold) ||(get_B_element(ax_len - i) < threshold && get_B_element(ax_len - i - stride) > threshold)){
        index = i;
        break;
      }
    }
  }
  //get omega value and convert to k.
  my_type omega_max;
  
  if(index != 0){
    omega_max = std::abs(get_om_axis_element(index));
  }else{
  //there either isn't a peak, or isn't waves or whatever. So we return axis end
    omega_max = get_om_axis_element(ax_len-1);
  }
  return get_k(omega_max, WAVE_WHISTLER);
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

void spectrum::copy_ids( data_array & src){
/** Copies ID fields from src array to this */

  strncpy(this->block_id, src.block_id, ID_SIZE);

  std::copy(src.time, src.time + 2, this->time);
  for(int i=0; i < 2; ++i) this->space[i] = src.space[i];

  g_angle_array.copy_ids(src);
  B_omega_array.copy_ids(src);

}

bool spectrum::check_ids( data_array & src){
/** Checks ID fields match src */

  bool err=false;
  if(strcmp(this->block_id, src.block_id) != 0) err =true;
  for(int i=0; i< 3; i++) if(src.time[i] != this->time[i]) err=true;
  for(int i=0; i < 2; ++i) if(this->space[i] != src.space[i]) err=true;
  err |= g_angle_array.check_ids(src);
  err |= B_omega_array.check_ids(src);
  return err;
}

data_array spectrum::copy_out_B(){
  
  data_array B_copy = this->B_omega_array;
  return B_copy;

}
data_array spectrum::copy_out_g(){
  
  data_array g_copy = this->g_angle_array;
  return g_copy;

}
