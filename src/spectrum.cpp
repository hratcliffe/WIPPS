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
#include <cstring>
#include <limits>
#include "plasma.h"
#include "spectrum.h"

/********Basic setup and allocation functions ****/
void spectrum::construct(){
/** \brief Generic spectrum contruction actions
*
*These are performed regardless of dimensions
*/
  function_type = 0;
  wave_id = 0;
  my_controller = nullptr;
  norm_B = 0;
  norm_g = nullptr;
  max_power=0.0;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));
  smooth=0;
}

void spectrum::init(){
/** \brief Cache spectrum values
*
*Fill cached values for max power, norms etc*/

  max_power = B_omega_array.maxval();
  calc_norm_B();
  size_t len = get_omega_length();
  for(size_t i = 0; i < len; i++){
    calc_norm_g(i);
  }
}

spectrum::spectrum(){
/** \brief Default constructor
*
*Sets up empty but valid spectrum
*/
  construct();
}

spectrum::spectrum(int n_om, int n_ang, bool separable){
/** \brief Construct spectrum
*
*Constructs a spectrum of the given dimensions, n_om, n_ang. If the B and angle dependences are seperable (separable = true) then this will be of the form B^2(omega) g(theta). Else it will be B^2(omega) g(omega, theta). Both have a B array of size n_om, the former has g of 1 x n_angles, the latter n_om x n_angles
@param n_om Number of frequency points in spectrum
@param n_ang Number of angle points
@param separable Whether g is function of omega
*/

  //Set up spectrum and the B data array
  construct();
  this->B_omega_array = data_array(n_om);

  //Setup correct shape/size g. Allocates memory to store the normalising values so that int_x g = 1 for each omega
  if(separable){
    this->g_angle_array = data_array(1, n_ang);
    g_is_angle_only = true;
    function_type = FUNCTION_DELTA;
    norm_g = (my_type *) calloc(n_om, sizeof(my_type));
  }else{
    this->g_angle_array = data_array(n_om, n_ang);
    g_is_angle_only = false;
    norm_g = (my_type *) calloc(n_om, sizeof(my_type));
  }
  //Cache the norms and maxes
  this->init();
}

spectrum::spectrum(std::string filename){
/** \brief Setup a spectrum from a file
*
*Reads spectrum from file, written using spectrum::write_to_file() or IDL write routines
@param filename Full filename to read from
*/

//First we grab the position of close block. Then we attempt to read two arrays. If we reach footer after first we error, or do not after second we warn.

  construct();
  std::fstream file;
  file.open(filename, std::ios::in|std::ios::binary);
  if(!file.is_open()) return;

  bool err = 0, cont = 1;
  size_t end_block=0, next_block=0, om_sz = 0;
  size_t jump_pos=0, g_first_dim = 0;

  end_block = get_file_endpos(file);
  file.seekg(0, std::ios::beg);

  std::vector<size_t> dims = B_omega_array.read_dims_from_file(file);

  if(dims.size() !=1){
    my_error_print("Invalid dimensions for B", mpi_info.rank);
    cont = 0;
  }else{
    my_print("B size "+mk_str(dims[0]), mpi_info.rank);
  
  }
  if(cont){
    //Now set up B correct size
    file.seekg(0, std::ios::beg);
    //return to start
    B_omega_array = data_array(dims[0]);
    om_sz = dims[0];
    err = B_omega_array.read_from_file(file);
  }
  if(cont && err){
    my_error_print("File read failed", mpi_info.rank);
    cont = 0;
  }
  if(cont){
    file.read((char*) &next_block, sizeof(size_t));
    if(next_block == end_block){
      //Early termination of file...
      my_error_print("Insufficent arrays in file", mpi_info.rank);
      cont =0;
    }
  }
  if(cont){
    file.seekg(-1*sizeof(size_t), std::ios::cur);
    //Last read was invalid, so we seek back

    jump_pos = (size_t)file.tellg();
    dims = g_angle_array.read_dims_from_file(file);

    if(dims.size() !=2){
      my_error_print("Wrong dimensions for g", mpi_info.rank);
      cont = 0;
    }else{
      my_print("g size "+mk_str(dims[0])+'x'+mk_str(dims[1]), mpi_info.rank);
    }
  
  }
  if(dims.size() > 1) g_first_dim = dims[0];//We need this later
  if(g_first_dim != 1 && g_first_dim != om_sz) cont = false;//g size doesn't match to B size
  if(cont){
    file.seekg(jump_pos);
    //Now set up g correct size
    if(g_first_dim){
      this->g_angle_array = data_array(1, dims[1]);
      g_is_angle_only = true;
      function_type = FUNCTION_DELTA;
      //Single row so only one norm
    }else{
      this->g_angle_array = data_array(dims[0], dims[1]);
      g_is_angle_only = false;
      //Norm each row
    }
    norm_g = (my_type *) calloc(om_sz, sizeof(my_type));
    err = g_angle_array.read_from_file(file);
  
    if(err){
      my_error_print("File read failed", mpi_info.rank);
      cont = 0;
    }
  }
  if(cont){
    file.read((char*) &next_block, sizeof(size_t));
    next_block= (size_t)file.tellg() - sizeof(size_t);
    //If we're done, this should be the
     if(next_block != end_block){
      //File is not done!
      my_error_print("Excess data in file", mpi_info.rank);
      cont = 0;
    }
  }
  if(cont){
    char id_in[ID_SIZE];
    file.seekg(end_block+sizeof(size_t));
    if(file){
      file.read(id_in, sizeof(char)*ID_SIZE);
      strcpy(this->block_id, id_in);
    }

    if(file) file.read((char*)this->time, sizeof(my_type)*2);
    if(file) file.read((char*)this->space, sizeof(size_t)*2);
    size_t wave_id_tmp = 0;
    if(file) file.read((char*)&wave_id_tmp, sizeof(size_t));
    wave_id = (int) wave_id_tmp;
    size_t function_tmp = 0;
    if(file) file.read((char*)&function_tmp, sizeof(size_t));
    function_type = (int) function_tmp;
    //If g_first_dim = 1 we probably have a test file. But this may not be specifying all the extra flags, so we keep but effectively ignore the function_type
    this->g_is_angle_only = true;
    if(g_first_dim > 1) this->g_is_angle_only = false;
    if(file) file.read((char*) &this->smooth, sizeof(size_t));
  }

  if(!cont){
    //Something went wrong. Clean up and quit
    g_angle_array = data_array();
    B_omega_array = data_array();
    //Return to size-less state
  }
  
  if(cont){
    //Before version v1.1 spectra were so that the total was the total power. For D etc we want rather the integral to be the power, so we convert files on read
    if(compare_as_version_string(read_wipps_version_string(filename), "v1.1", true) < 0){
      this->convert_FFT_to_integral();
    }
    this->init();
  }
  return;

}

spectrum::~spectrum(){
/** \brief Delete spectrum
*
*Free any allocated memory
*/
  if(norm_g) free(norm_g);
}

/********Technical stuff making my_array a proper "object" ****/
spectrum & spectrum::operator=(const spectrum& src){
/** \brief Copy assignment
*
*Sets this equal to a (deep) copy of source, i.e duplicates the B and g arrays and all other fields
@param src Spectrum to copy from
@return Copy of input spectrum
*/
  
  //Trap self-assign or bad copy before destructing
  if(&src == this || !src.is_good()) return *this;
  if(this->norm_g){
    free(norm_g);
    norm_g =nullptr;
  }
  construct();

  //(Deep) copy all fields
  size_t g_sz = src.get_omega_length();
  norm_g = (my_type *) calloc(g_sz, sizeof(my_type));
  std::copy(src.norm_g, src.norm_g + g_sz, this->norm_g);
  this->g_angle_array = src.g_angle_array;
  this->B_omega_array = src.B_omega_array;
  this->copy_ids(src.B_omega_array);
  this->copy_tags(src);

  return *this;
}

spectrum::spectrum(const spectrum &src){
/** \brief Copy constructor
*
*(Deep) copy src to a new instance.
@param src Spectrum to copy from
*/

  //Initialise fields
  construct();

  //If src is not good can't set anything
  if(!src.is_good()) return;

  //Copy all fields
  size_t g_sz = src.get_omega_length();
  norm_g = (my_type *) calloc(g_sz, sizeof(my_type));
  std::copy(src.norm_g, src.norm_g + g_sz, this->norm_g);
  this->g_angle_array = src.g_angle_array;
  this->B_omega_array = src.B_omega_array;
  this->copy_ids(src.B_omega_array);
  this->copy_tags(src);
}

bool spectrum::operator==(const spectrum &rhs)const{
/** \brief Equality operator
*
* Check this is equal to rhs. Since copies are always deep, we check values, not data pointers. We ignore the derived things such as smooth and norm_g
@param rhs Spectrum to compare to
@return True if equal, false else
*/

  if(this->my_controller != rhs.my_controller) return false;
  if(this->B_omega_array != rhs.B_omega_array) return false;
  if(this->g_angle_array != rhs.g_angle_array) return false;
  if(this->check_ids(rhs.B_omega_array)) return false;
  if(this->check_tags(rhs)) return false;

  return true;
}

/********Setup helper functions ****/
void spectrum::make_angle_axis(){
/** \brief Generate angle axis
*
* Create uniform angle axis from TAN_MIN to TAN_MAX (see support.h). In general these are interpreted as tan(theta) values.
*/
  calc_type res = (TAN_MAX - TAN_MIN)/get_angle_length();
  int offset = -TAN_MIN/res;
  g_angle_array.make_linear_axis(1, res, offset);

}

bool spectrum::make_angle_distrib(my_type std_dev){
/** \brief Generate angle axis and distribution
*
*Generates an angle axis linear spaced in tan theta between TAN_MIN and TAN_MAX. Then generates and fills the angular spectrum according to function specified by function_type member variable. Options are in support.h and are FUNCTION_DELTA: delta function with peak at 0 angle. FUNCTION_GAUSS: Gaussian with centre at 0 angle. FUNCTION_ISO: Isotropic distribution over range considered. The integrals are either 0.5 for 0 to TAN_MAX or 1 for -TAN_MAX to TAN_MAX. Other combinations are not guaranteed
@param std_dev Standard deviation of angular spectrum (Gauss only)
@return 0 (success), 1 (error)
*/


  if(!g_is_angle_only){
    my_error_print("Angular distrib is not a function. Returning", mpi_info.rank);
    return 1;
  }
  if(!g_angle_array.is_good()){
    my_error_print("Angular array invalid. Returning", mpi_info.rank);
    return 1;
  }
  make_angle_axis();
  size_t len = get_angle_length();
  my_type val;

  if(function_type == FUNCTION_DELTA){
  
    for(size_t i=0; i<len; ++i) set_g_element(i,0.0);
    val = 1.0/(get_ang_axis_element(1)-get_ang_axis_element(0));
    int zero = where_angle(0.0);
    //Set_element checks bnds automatically
    set_g_element(zero, val);

  }else if(function_type == FUNCTION_GAUSS){
    my_type ax_el;
    my_type norm;
    norm = 1.0/ (std::sqrt(2.0*pi) * std_dev);
    for(size_t i=0; i<len; ++i){
      ax_el = get_ang_axis_element(i);
      val = std::exp( -0.5 * std::pow(ax_el/std_dev, 2)) * norm;
      set_g_element(i,val);
    }
  }else if(function_type == FUNCTION_ISO){

    val = 1.0/ (TAN_MAX - TAN_MIN)*get_angle_length()/(get_angle_length()-1);
    if(std::abs(TAN_MIN) < GEN_PRECISION) val = val / 2.0;
    for(size_t i=0; i<len; ++i) set_g_element(i,val);

  }else{
  
    my_error_print("Invalid function type. Returning", mpi_info.rank);
    return 1;

  }
  
  return 0;

}

#ifdef RUN_TESTS_AND_EXIT
void spectrum::make_test_spectrum(int angle_type, bool two_sided, my_type om_ce, my_type std_dev){
/** \brief Generate dummy spectrum
*
*Makes a basic spectrum object with suitable number of points, and Gaussian(s) centred at fixed k/freq and x value 
@param angle_type Function to use for angular distrib 
@param two_sided Whether to generate symmetric spectrum or one-sided 
@param om_ce Plasma frequency to use
@param std_dev Standard deviation of angular functional form (where applicable)
*/

  if(!this->is_good()){
    my_error_print("Spectrum object invalid, returning", mpi_info.rank);
    return;
  }
  if(!this->g_is_angle_only){
    my_error_print("Test spectrum must be separable", mpi_info.rank);
    return;
  }

  //Setup the ID fields
  char id[ID_SIZE] = "ex";
  this->set_ids(0, 1, 0, 1, WAVE_WHISTLER, id, angle_type);
  
  //Setup axes and required angular distrib
  size_t len0 = get_omega_length();
  //res to cover range from either 0 or -max to max in len0 steps
  my_type res = om_ce*(1.0+(my_type)two_sided)/(my_type)len0;
  //Generate axes for one or two-sided
  for(size_t i = 0; i < len0; i++) set_om_axis_element(i, res*((my_type)i - (my_type)two_sided *(my_type)len0/2.));
  make_angle_distrib(std_dev);

  //Set values to match our test spectrum data
  my_type centre, width, background = 0.0;
  centre = om_ce*14.0/17.0, width=0.1*centre;
  
  my_type * data_tmp, * data_ptr;
  data_ptr = (my_type *) malloc(len0*sizeof(my_type));
  data_tmp = data_ptr;
  my_type ax_val = 0.0;
  if(two_sided){
    //Fill a two-sided spectrum with peaks where specified
    for(size_t i=0; i<=len0/2; i++, data_tmp++){
      ax_val = get_om_axis_element(i);
      *(data_tmp) = exp(-pow((ax_val + centre), 2)/width/width) + background;
    }
    data_tmp--;
    //we've gone one past our termination condition...
    for(size_t i = 1; i < len0/2; i++) *(data_tmp + i) = *(data_tmp - i);
  }else{
  
    for(size_t i = 0; i < len0; i++, data_tmp++){
      ax_val = get_om_axis_element(i);
      *(data_tmp) = exp(-pow((ax_val - centre), 2)/width/width) + background;
    }
  }
  //Copy the data into B
  this->B_omega_array.populate_data(data_ptr, len0);
  free(data_ptr);

  max_power = 1.0;
  //Store value of maximum
  this->init();

}
#endif

bool spectrum::generate_spectrum(data_array &parent, int om_fuzz, int angle_type, my_type std_dev, data_array * mask){
/**\brief Generate spectrum from data
*
*Takes a parent data array and generates the corresponding spectrum. Windows using the specified wave dispersion and integrates over frequency using om_fuzz percent band. Axes are copied from the parent. If the spectrum is of separable type (g_is_angle_only = true), the angular distribution is generated with functional form specified by angle_type.  IMPORTANT: when using real angular data we roughly fuzz around the correct k values, but this is not uniform! Non-smooth or rapidly varying data may give odd results 
@param parent Data array to read from. Spectrum will have the same units as this 
@param om_fuzz Band width around dispersion curve in percent of central frequency 
@param angle_type Angular distribution functional form 
@param std_dev Standard deviation of angular functional form (where applicable)
@param mask (optional) data_array matching sizes of parent, will be filled with the masking array used for spectrum generation. If nullptr or nothing is supplied, no mask is output. 
@return 0 for successful calculation, 1 for error
\todo 2-d and 3-d extractions don't quite agree at k=0. factor ~10 and variations near 0
*/

  if(!this->is_good()){
    my_error_print("Spectrum object invalid. Returning", mpi_info.rank);
    return 1;
  }
  if(parent.get_dims() == 0 || parent.get_dims(0) != this->get_omega_length()){
    my_error_print("0-dim of parent must match omega_length", mpi_info.rank);
    return 1;
  }


  if(parent.is_good() && g_is_angle_only){
    //Require parent to be 2-D
    if(parent.get_dims() != 2){
      my_error_print("Functional angle form requires 2-D parent", mpi_info.rank);
      return 1;
    }

    //First we read ids from parent
    this->copy_ids(parent);
    this->set_extra_ids(WAVE_WHISTLER, angle_type);
    //...and omega axis
    size_t len;
    const my_type * om_ax = parent.get_axis(1, len);

    //Now we loop across x, calculate the wave cutout bounds, and total, putting result into data

    int j;
    int low_bnd, high_bnd;
    my_type om_disp, max=0.0;
    my_type tolerance = om_fuzz/100.0;
    my_type total = 0.0;

    calc_type sgn;
    //One sided axis (runs 0 to max) rather than -max to max. These are the only options
    bool one_sided = (std::abs(parent.get_axis_element(1,0)) < std::abs(parent.get_axis_element(1,len/2)));

    for(size_t i = 0; i < get_omega_length(); ++i){
      //Convert parent k-axis value to omega, preserving sign
      om_disp = get_omega(parent.get_axis_element(0,i), WAVE_WHISTLER);
      if(!one_sided && i<= get_omega_length()/2) sgn = -1.0;
      else sgn = 1.0;
      om_disp *= sgn;
      //Set spectrum axis to this value
      set_om_axis_element(i, om_disp);
      
      //Lookup the parent omega axis values for both ends of band
      low_bnd = where(om_ax, len, om_disp *(1.0-sgn*tolerance));
      high_bnd = where(om_ax, len, om_disp *(1.0+sgn*tolerance));
      
      //Adjust for low off bottom or high off top
      //Low off top or high off bottom will still be caught and skipped
      if(low_bnd < 0 && om_disp *(1.0-sgn*tolerance) < parent.get_axis_element(1,0)) low_bnd = 0;
      if(high_bnd < 0 && om_disp *(1.0+sgn*tolerance) > parent.get_axis_element(1,len-1)) high_bnd = len-1;

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
    //Generate the angular distrib function
    make_angle_distrib(std_dev);
    this->max_power = max;

  }else if(parent.is_good() && parent.get_dims()==3){
  //This will do a rough integral around the correct omega. Note it is not a neat sampling of the domain, so may give odd results for non-smooth cases
  
    this->copy_ids(parent);
    this->set_extra_ids(WAVE_WHISTLER, FUNCTION_NULL);

    //Convert om_fuzz from percentage to fraction
    my_type tolerance = om_fuzz/100.0;

    //Get the k_x and k_y axes for lookup
    size_t len_x, len_y, len_om;
    const my_type * kx_ax = parent.get_axis(0, len_x);
    const my_type * ky_ax = parent.get_axis(1, len_y);
    const my_type * om_ax = parent.get_axis(2, len_om);
    
    //One sided axis (runs 0 to max) rather than -max to max. These are the only options
    bool one_sided = (std::abs(om_ax[0]) < std::abs(om_ax[len_om/2]));

    //Create angular axis for spectrum
    make_angle_axis();
    
    my_type k, tantheta, tmp, tmp_sum, max_power=0.0, theta, om_disp;
    int kx_high, ky_low, ky_high, i_sgn=1, om_low, om_high;

    //Now we do a double loop over k_x and angle
    for(size_t i = 0; i< len_x; i++){
      //Signs so we always integrate from k_low to k_high, if one sided is always 1
      if(!one_sided) i_sgn = (i<len_x/2 ? -1: 1);

      //Use the k_x axis to generate an omega one. Thus if fuzz=0 we match 1-1
      om_disp = get_omega(parent.get_axis_element(0,i), WAVE_WHISTLER);
      set_om_axis_element(i, om_disp*i_sgn);

      //We can assume either +ve and -ve k (include i_sgn in kx_high, or omega (include in om_low and high)
      
      tmp_sum = 0.0;
      for(size_t j = 0; j< get_angle_length(); j++){
        tantheta = get_ang_axis_element(j);
        theta = atan(tantheta);
//----------Fuzz omega--------------------
        om_low = where(om_ax, len_om, om_disp*(1.0-i_sgn*tolerance)*i_sgn);
        om_high = where(om_ax, len_om, om_disp*(1.0+i_sgn*tolerance)*i_sgn);
        if(om_low < 0) om_low = 0;
        if(om_disp*(1.0+tolerance) < *om_ax) om_high = 0;
        else if(om_high < 0) om_high = len_om - 1;

//-------No fuzz in k_x -----------------------
        k = get_k(om_disp, WAVE_WHISTLER, 0, theta);
        
        kx_high = where(kx_ax, len_x, k*cos(theta));
        
        if(k*cos(theta)< *(kx_ax) || std::abs(k)<GEN_PRECISION || k*cos(theta)> *(kx_ax+len_x -1)) kx_high = -1;
        else if(kx_high <0) kx_high = len_x-1;
       
//-------- k_y is fuzzy----------
        k = get_k(om_disp*(1.0+tolerance), WAVE_WHISTLER, 0, theta);
        ky_high = where(ky_ax, len_y, k*sin(theta));
        if(std::abs(k)<GEN_PRECISION || std::abs(k*sin(theta))> *(ky_ax+len_y -1)) ky_high = -1;
        else if(ky_high <0) ky_high = len_y-1;

        k = get_k(om_disp*(1.0-tolerance), WAVE_WHISTLER, 0, theta);
        ky_low = where(ky_ax, len_y, k*sin(theta));
        if(ky_low < 0) ky_low = 0;
        
//-------Total on small range and create mask if requested-----
        tmp = 0.0;
        if(kx_high != -1 && ky_high != -1){
          for(int ii=om_low; ii<=om_high; ii++){
            for(int jj=ky_low; jj<=ky_high; jj++){
              tmp += parent.get_element(kx_high, jj, ii);
            }
            if(!one_sided){
            //Also include negative k_y
              for(int jj=len_y-ky_high; jj<(long)len_y-ky_low; jj++){
                tmp += parent.get_element(kx_high, jj, ii);
              }
            }
          }
        }
        if(mask){
          if(kx_high != -1 && ky_high != -1){
            for(int ii=om_low; ii<=om_high; ii++){
              for(int jj=ky_low; jj<=ky_high; jj++){
                if(!one_sided) mask->set_element(kx_high, jj, ii, mask->get_element(kx_high, jj, ii)+0.5);
                else mask->set_element(kx_high, jj, ii, mask->get_element(kx_high, jj, ii)+1.0);
              }
              if(!one_sided){
              //Also include negative k_y
                for(int jj=len_y-ky_high; jj<(long)len_y-ky_low; jj++){
                  mask->set_element(kx_high, jj, ii, mask->get_element(kx_high, jj, ii)+0.5);
                }
              }
            }
          }
       }
//------------------
        //Tmp is now the convolved B(omega)g(omega, theta). We dump it into g as is for later norming and add to the sum which will end up as B.
        if(one_sided) tmp *= 2.0;//If one_sided axis, assume there's the same power at -ve omega too
        set_g_element(i, j, tmp);
        tmp_sum += tmp;
      }
      //Now the angles are done, we norm. g and set B
      //Use some small value to prevent a divide by zero
      if(std::abs(tmp_sum) > tiny_my_type){
        for(size_t j = 0; j< get_angle_length(); j++) set_g_element(i, j, get_g_element(i, j)/tmp_sum);
      }
      set_B_element(i, tmp_sum);
      //Store the max power for later
      if(tmp_sum > max_power) max_power = tmp_sum;

    }
    this->max_power = max_power;
  }else{
    return 1;
  }
  
  //Correct the spectrum to be properly X^2(omega) dOmega which is what we ultimately want
  this->convert_FFT_to_integral();

  //Calc and store normalising info etc
  this->init();

  return 0;

}

void spectrum::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[ID_SIZE], int function_type){
/**\brief Set id fields
*
*Sets the time and space ranges, wave type etc attached to this spectrum. Times should be in seconds. Space in terms of grid points.
@param time1 Initial time of data used
@param time2 End time of data used
@param space1 Start index of space range of data used
@param space2 End index of space range of data used
@param wave_id Wave type (see support.h)
@param block_id Name of block used to derive this spectrum
@param function_type Functional form of spectrum in angle
*/

  this->time[0] = time1;
  this->time[1] = time2;
  this->space[0] = space1;
  this->space[1] = space2;
  strncpy(this->block_id, block_id, ID_SIZE);
  this->wave_id = wave_id;
  this->function_type = function_type;
}

void spectrum::set_extra_ids(int wave_id, int function_type){
/**\brief Set id fields
*
*Sets the wave type and angle type attached to spectra (i.e. the ids that can't be inherited from parent data)
@param wave_id Wave type (see support.h)
@param function_type Functional form of spectrum in angle
*/

  this->wave_id = wave_id;
  this->function_type = function_type;
}

void spectrum::copy_ids( const data_array & src){
/**\brief Copy id fields
*
*Copies ID fields from src array to this. ID fields are the block_id, and the space and time values. These are attached to the spectrum AND to the B and g arrays it holds. See also spectrum::copy_tags()
@param src Array to copy ids from
*/

  strncpy(this->block_id, src.block_id, ID_SIZE);
  std::copy(src.time, src.time + 2, this->time);
  for(int i=0; i < 2; ++i) this->space[i] = src.space[i];
  g_angle_array.copy_ids(src);
  B_omega_array.copy_ids(src);

}

bool spectrum::check_ids(const data_array & src)const{
/** \brief Check ids match
*
*Checks ID fields match src. ID fields are the block_id, and the space and time values. See also spectrum::check_tags()
@param src Array to check ids against
@return True if ids are equal, false else
*/

  bool err = false;
  if(strcmp(this->block_id, src.block_id) != 0) err = true;
  for(int i=0; i< 2; i++) if(src.time[i] != this->time[i]) err = true;
  for(int i=0; i < 2; ++i) if(this->space[i] != src.space[i]) err = true;
  err |= g_angle_array.check_ids(src);
  err |= B_omega_array.check_ids(src);
  return err;
}

void spectrum::copy_tags(const spectrum & src){
/** \brief Copy tags
*
*Copies tag fields from src array to this. Tags are the g_is_angle_only, the function_type and the wave_id. See also spectrum::copy_ids()
@param src Spectrum to copy tags from
*/

  this->g_is_angle_only = src.g_is_angle_only;
  this->function_type = src.function_type;
  this->wave_id = src.wave_id;

}

bool spectrum::check_tags(const spectrum & src)const{
/** \brief Check tags
*
*Check tag fields match src. Tags are the g_is_angle_only, the function_type and the wave_id. See also spectrum::check_ids()
@param src Spectrum to compare tags against
@return True if equal, false else
*/

  bool err = false;
  if(this->g_is_angle_only != src.g_is_angle_only) err = true;
  if(this->function_type != src.function_type) err = true;
  if(this->wave_id != src.wave_id) err = true;

  return err;
}

/********Access helper functions ****/
int spectrum::where_omega(my_type omega){
/** \brief Find where omega axis exceeds omega
*
*Finds where spectrum frequency axis exceeds the value of omega.
@param omega Value of omega to find
@return Index where omega axis exceeds given value
*/
  int index;
  size_t len;
  const my_type * om_axis = get_omega_axis(len);
  index = where(om_axis, len, omega);
  return index;

}

int spectrum::where_angle(my_type ang){
/** \brief Find where angle axis exceeds ang
*
*Finds where spectrum angle exceeds the value of ang
@param ang Value of ang to find
@return Index where angle axis exceeds given value
*/
  int index;
  size_t len;
  const my_type * ang_axis = get_angle_axis(len);
  index = where(ang_axis, len, ang);
  return index;
}

my_type spectrum::get_omega(my_type k, int wave_type, bool deriv, my_type theta){
/** \brief Gets omega for given k
*
*Uses dispersion relation for given wave_type to convert k to omega. Calls to plasma because approximations for density etc etc should be made there. 
@param k Wavenumber 
@param wave_type Wave species 
@param deriv Return v_g instead 
@param theta Wave normal angle
@return Value of omega for given k and wave type etc
*/
  if(my_controller && (my_controller->is_good())) return (my_type) my_controller->get_plasma().get_dispersion(k, wave_type, 0, deriv, theta);
  else return 0.0;
}

my_type spectrum::get_k(my_type omega, int wave_type, bool deriv, my_type theta) {
/** \brief Gets k for given omega
*
* Uses dispersion relation for given wave_type to convert omega to k. Calls to plasma because approximations for density etc etc should be made there.
@param omega Frequency 
@param wave_type Wave species 
@param deriv Return v_g instead 
@param theta Wave normal angle
@return Value of k for given omega and wave type etc
*/
  if(my_controller && my_controller->is_good()) return (my_type) my_controller->get_plasma().get_dispersion(omega, wave_type, 1, deriv, theta);
  else return 0.0;
}

/********Spectrum operation helpers ****/

void spectrum::convert_FFT_to_integral(){
/** \brief Convert B data to integral form
*
* Converts the data in B^2 so that its integral is power, assuming the source is an FFT. FFT's give the power at a frequency, so we have to divide by d_omega
*/
  my_type d_om;
  size_t len = get_omega_length();
  if(len >= 2){
    d_om = std::abs(get_om_axis_element(1) - get_om_axis_element(0));
    set_B_element(0, get_B_element(0)/d_om/d_om);
  }
  for(size_t i = 1; i < len; i++){
    d_om = std::abs(get_om_axis_element(i) - get_om_axis_element(i - 1));
    set_B_element(i, get_B_element(i)/d_om/d_om);
  }

}

bool spectrum::calc_norm_B(){
/** \brief Calculate norming of B(w)
*
* Calculate the total square integral of values over range, int_{om_min}^{om_max} B^2(omega) d omega.
@return 0 (success), 1 (error)
*/

  int len = get_omega_length();
  my_type * d_axis = (my_type *) calloc(len, sizeof(my_type));
  my_type * data = (my_type *) malloc(len*sizeof(my_type));

  //Assemble d_om and data in contiguous memory
  for(int i=0; i<len-1; i++) d_axis[i] = get_om_axis_element(i+1) - get_om_axis_element(i);
  for(int i=0; i<len; i++) data[i] = get_B_element(i);
  
  norm_B = integrator(data, len, d_axis);
  free(d_axis);
  free(data);
  return 0;
}

bool spectrum::calc_norm_g(size_t om_ind){
/** \brief Normalise g_w(x)
*
*Calculate the norm of g used in e.g. denom of Albert \cite Albert2005 Eq 3 or calc'd in Derivations.tex \cite Derivations. Contains one value for each omega entry.
@param om_ind Omega index to calculate norm at 
@return 0 (success), 1 (error e.g. out of range)
*/

  size_t len = get_angle_length();
  if(!my_controller) return 1;
  plasma plas = my_controller->get_plasma();

  my_type * d_axis = (my_type *) calloc(len, sizeof(my_type));
  my_type * integrand = (my_type *) calloc(len, sizeof(my_type));

  for(size_t i = 0; i < len-1; i++) d_axis[i] = get_ang_axis_element(i+1) - get_ang_axis_element(i);
  //Construct dx axis for integration. Note x = tan theta

  mu_dmudom my_mu;
  my_type x, psi, omega = std::abs(get_om_axis_element(om_ind)), g_el;
  
  size_t lena = get_omega_length();
  //Check omega index is in range
  if(om_ind >= lena) return 1;

  for(size_t i = 0; i < len; i++){
    x = get_ang_axis_element(i);
    psi = atan(x);
    my_mu = plas.get_mu(omega, psi);
    if(!my_mu.err){
      if(g_is_angle_only){
        g_el = get_g_element(i);
      }else{
        g_el = get_g_element(om_ind, i);
      }
      integrand[i] = g_el * x / std::pow((std::pow(x, 2)+1.0), 1.5)*std::pow(my_mu.mu, 2) * std::abs( my_mu.mu + omega*my_mu.dmudom);

    }else{
      integrand[i] = 0.0;
    }
    //product of g(theta) * x (x^2+1)^-(3/2) * mu^2 |mu+omega dmu/domega|
    //See Derivations#Evaluation_of_G2 for details
  }
  
  //integrate
  my_type norm_g_tmp = integrator(integrand, len, d_axis);

  //Soften so can never be zero. Should only approach this if g_om(theta) is everywhere zero, so set simply to something tiny
  if(norm_g_tmp < std::numeric_limits<my_type>::min()){
    norm_g_tmp = std::numeric_limits<my_type>::min();
  }
  
  norm_g[om_ind] = norm_g_tmp;
  //clean up
  free(d_axis);
  free(integrand);
  return 0;
}

void spectrum::smooth_B(int n_pts){
/** \brief Smooth B
*
*Apply a box-car smoothing to B with specified number of pts. Store n_pts in smooth field
@param n_pts Boxcar smoothing width to apply
*/

  this->smooth = n_pts;
  B_omega_array.smooth_1d(n_pts);
}

bool spectrum::truncate_om(my_type om_min, my_type om_max){
/** \brief Truncate omega distribution at om_min and om_max.
*
*Zeros all elements outside the range [om_min, om_max]. Zeros are ignored. om_min must be < om_max. Om_min or max out of axis range does nothing on that end. NB B is renormalised after the truncation
@param om_min Minimum physical omega to truncate at
@param om_max Maximum physical omega to truncate at
@return 0 (success), 1 (range error)
*/

  if(om_min >= om_max){
    my_error_print("Invalid omega range, aborting truncate", mpi_info.rank);
    return 1;
  }

  int index = -1;
  size_t len = get_omega_length();
  size_t len_ang = get_angle_length();

  if(om_min > get_om_axis_element(1)){
    index = where_omega(om_min);
    if(index != -1){
      for(int i=0; i< index; i++){
        set_B_element(i, 0.0);
        if(!g_is_angle_only){
          for(size_t j = 0; j< len_ang; j++){
            set_g_element(i, j, 0.0);
          }
        }
      }
    //Zero up to om_min
    }
  }

  if(om_max < get_om_axis_element(len-1)){
    index = where_omega(om_max);
    if(index != -1){
      for(size_t i = index; i< len; i++){
        set_B_element(i, 0.0);
        if(!g_is_angle_only){
          for(size_t j = 0; j< len_ang; j++){
            set_g_element(i, j, 0.0);
          }
        }
      }
    //Zero after om_max
    }
  }

  calc_norm_B();
  for(size_t i = 0; i< get_omega_length(); i++){
    calc_norm_g(i);
  }

  //Re-do normalisation
  return 0;
}

bool spectrum::truncate_x(my_type x_min, my_type x_max){
/** \brief Truncate angle distribution at x_min and x_max.
*
*Zeros all elements outside the range [x_min, x_max]. x_min must be < x_max. If x_min or max are out of range, nothing is done at that end
@param x_min Minimum tan theta to truncate at
@param x_max Maximum tan theta to truncate at
@return 0 (success), 1 (range error)
*/

  if(x_min >= x_max){
    my_error_print("Invalid x range, aborting truncate", mpi_info.rank);
    return 1;
  }

  int index = -1;
  size_t len = get_angle_length(), om_len;

  om_len =  (g_is_angle_only? 1 : get_omega_length());
  
  if(x_min > TAN_MIN){
    index = where_angle(x_min);
    if(index != -1){
      for(size_t i=0; i< (size_t)index; i++){
        for(size_t j=0; j< om_len; j++){
          set_g_element(j, i, 0.0);
        }
      }
    }
  
  }
  if(x_max < TAN_MAX){
    index = where_angle(x_max);
    if(index != -1){
      for(size_t i=index; i< len; i++){
        for(size_t j=0; j< om_len; j++){
          set_g_element(j, i, 0.0);
        }
      }
    }
  }
  for(size_t i = 0; i< get_omega_length(); i++){
    calc_norm_g(i);
  }
  return 0;

}

calc_type spectrum::check_upper(){
/** \brief Check upper k limit of spectral power
*
* Checks the upper bound of region of significant spectral power, i.e. above SPECTRUM_THRESHOLD*peak_power
@return The physical wavenumber where spectrum drops below threshold
*/

  size_t stride = 1;
  size_t ax_len = get_omega_length();
  size_t len = ax_len/2 - stride;
  my_type threshold = SPECTRUM_THRESHOLD*this->max_power;
  size_t index = 0;

  if(get_om_axis_element(1) < get_om_axis_element(2)){
  //Axis presumed to be +ve only, move in from top to find first rise above threshold
    index = ax_len;
    for(size_t i=1; i< ax_len; i+=stride){
      if((get_B_element(ax_len - i) < threshold && get_B_element(ax_len - i - stride) > threshold)){
        index = ax_len - i;
        break;
      }
    }
  }else{
  //Axis presumed to be symmetrical. Move in from both ends
    for(size_t i=1; i< len; i+=stride){
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
@return The physical axis value where the spectral peak occurs
*/

  calc_type value = -1.0, tmp;
  int index;
  for(size_t i = 0; i < get_omega_length(); ++i){
    tmp = get_B_element(i);
    if(tmp > value){
      index = i;
      value = tmp;
    }
  }

  return get_om_axis_element(index);

}

/********Spectrum mutation helpers ****/

void spectrum::apply(spectrum::part subarr, std::function<my_type(my_type arg)> func){
/** Apply a function to each element of selected spectrum part. func must take and return a my_type or type convertible to this
@param subarr Which part to apply to, see spectrum::part
@param func Function to apply 
*/

  if(subarr == spectrum::part::B) B_omega_array.apply(func);
  else if(subarr == spectrum::part::ang) g_angle_array.apply(func);
  this->renormalise();
}

void spectrum::apply(spectrum::part subarr, std::function<my_type(my_type arg, my_type arg2)> func, my_type arg){
/** Apply a function to each element of selected spectrum part. func must take two my_type and return one my_type or type convertible to this
@param subarr Which part to apply to, see spectrum::part
@param func Function to apply
@param arg Second argument to function 
*/

  if(subarr == spectrum::part::B) B_omega_array.apply(func, arg);
  else if(subarr == spectrum::part::ang) g_angle_array.apply(func, arg);
  this->renormalise();
}

/********File IO ****/
bool spectrum::write_to_file(std::fstream &file){
/** \brief Write to file
*
* Spectra are written by writing out the B array, the g array, and then writing a single closing footer containing the id values again.
@param file Filestream to write to
@return 0 (success), 1 (error)
*/

  if(!file.is_open() || !this->is_good()) return 1;

  //Write B and g without closing files.
  bool err, write_err=0;
  err = B_omega_array.write_to_file(file, false);
  err |= g_angle_array.write_to_file(file, false);

  //Grab the start position of the final "footer" block
  size_t ftr_start = (size_t) file.tellg();
  size_t next_location = ftr_start+ sizeof(char)*ID_SIZE +sizeof(size_t)*6 + sizeof(my_type)*2;

  //Write the footer
  file.write((char*) & next_location, sizeof(size_t));
  file.write(block_id, sizeof(char)*ID_SIZE);
  file.write((char*) &time, sizeof(my_type)*2);
  file.write((char*) &space, sizeof(size_t)*2);
  //Convert wave id to size_t before writing
  size_t wave_id_tmp = this->wave_id;
  file.write((char*) &wave_id_tmp, sizeof(size_t));
  size_t function_tmp = this->function_type;
  file.write((char*) &function_tmp, sizeof(size_t));

  file.write((char*) & smooth, sizeof(size_t));

  //Check for errors
  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_error_print("Error writing offset positions", mpi_info.rank);
  //Write final closing value
  file.write((char*) & ftr_start, sizeof(size_t));

  return err;
}

bool spectrum::read_from_file(std::fstream &file){
/** \brief Initialise spectrum from file
*
*Reads a dump file which is expected to contain two arrays, first B then g, as written by spectrum->write_to_file, and constructs spectrum from data
@param file Filestream to read from
@return 0 (success), 1 (error)
*/

//First we grab the position of close block. Then we attempt to read two arrays. If we reach the closing block position after reading the first we error; if we do not after the second we warn. If the arrays in file are the wrong size, we have to quit with an error
  bool err = 0;
  size_t end_block=0, next_block=0;
  file.seekg(-1*sizeof(size_t), file.end);
  file.read((char*) &end_block, sizeof(size_t));

  //Rewind file and read B array
  file.seekg(0, std::ios::beg);

  err = B_omega_array.read_from_file(file);
  if(err){
    my_error_print("File read failed", mpi_info.rank);
    return err;
  }
  //Rewind and grab the next_position field to check if we're at the footer yet
  file.seekg(-1*sizeof(size_t), std::ios::cur);
  file.read((char*) &next_block, sizeof(size_t));
  if(next_block == end_block){
    //Early termination of file...
    my_error_print("Insufficent arrays in file", mpi_info.rank);
    B_omega_array = data_array();
    return 1;
  }

  //Read the g array
  err = g_angle_array.read_from_file(file);
  if(err){
    my_error_print("File read failed", mpi_info.rank);
    return err;
  }
  //Rewind and grab the next_position field to check if we're at the footer yet
  file.seekg(-1*sizeof(size_t), std::ios::cur);
  file.read((char*) &next_block, sizeof(size_t));
  if(next_block != end_block){
    //File is not done!
    my_error_print("Excess arrays in file", mpi_info.rank);
  }
  //Read the footer data, until we run out
  char id_in[ID_SIZE];
  file.seekg(end_block+sizeof(size_t));
  if(file) file.read(id_in, sizeof(char)*ID_SIZE);
  strcpy(this->block_id, id_in);
  if(file){
    char id_in[ID_SIZE];
    file.seekg(end_block+sizeof(size_t));
    if(file) file.read(id_in, sizeof(char)*ID_SIZE);
    strcpy(this->block_id, id_in);

    if(file) file.read((char*)this->time, sizeof(my_type)*2);
    if(file) file.read((char*)this->space, sizeof(size_t)*2);
    size_t wave_id_tmp = 0;
    if(file) file.read((char*) &wave_id_tmp, sizeof(size_t));
    wave_id = (int) wave_id_tmp;
    size_t function_tmp = 0;
    if(file) file.read((char*) &function_tmp, sizeof(size_t));
    function_type = (int) function_tmp;
    this->g_is_angle_only = true;
    if(function_type == FUNCTION_NULL) this->g_is_angle_only = false;
    if(file) file.read((char*) &this->smooth, sizeof(size_t));
  }

  return err;
}

/********Data release (for testing) ****/
data_array spectrum::copy_out_B(){
/** \brief Return a copy of B array
*
* Make a copy of the B part of data
@return A data array containing a copy of the B data
*/
  data_array B_copy = this->B_omega_array;
  return B_copy;

}
data_array spectrum::copy_out_g(){
/** \brief Return a copy of g array
* Make a copy of the g part of data
@return A data array containing a copy of the g data
*/
  
  data_array g_copy = this->g_angle_array;
  return g_copy;

}

/********Main spectral calculations ****/

calc_type get_G1(spectrum * my_spect, calc_type omega){
/** \brief G1 from \cite Albert2005
*
*Gets the value of B^2(w) (interpolated if necessary) and the normalising constant from norm_B. NB NB we omit the \f$ \Delta\omega \f$ factor since it cancels in the D calc
@param my_spect Input spectrum to work on
@param omega Frequency to eval. at
@return Normalised abs-square wave power
\caveat We assume omega symmetry here and just take abs(omega). If this is changed, best to add a omega_symm flag to spectra and select based on that. Note also that this routine uses the given spectrums omega range and assumes that beyond this there is no wave power
*/

  omega = std::abs(omega);//Assume symmetric in omega

  my_type tmpB2;
  if(my_spect->get_norm_B() == 0.0) my_spect->calc_norm_B();
  size_t len = my_spect->get_omega_length(), offset;
  my_type data_bit[2], ax_bit[2], ax_val;
  
  offset = my_spect->get_om_axis_index_from_value(omega);
  //Interpolate if possible, else use the end
  if(offset > 0 && offset < len){
    data_bit[0] = my_spect->get_B_element(offset-1);
    data_bit[1] = my_spect->get_B_element(offset);
    ax_bit[0] = my_spect->get_om_axis_element(offset-1);
    ax_bit[1] = my_spect->get_om_axis_element(offset);
    ax_val = (my_type) omega;
    tmpB2 = interpolate_linear(ax_bit, data_bit, ax_val);
  }else if(offset == 0){
    //we're right at end, can't meaningfully interpolate, use raw
    tmpB2 = my_spect->get_B_element(0);
  }else{
    //offset <0 or >= len, value not found
    tmpB2 = 0.0;
  }

  //Add norm. constant
  return tmpB2/my_spect->get_norm_B();

}

calc_type get_G2(spectrum * my_spect, calc_type omega, calc_type x){
/** \brief Get G2 from \cite Albert2005
*
* Gets the value of g(w, x) and the normalising constant from norm_g 
@param my_spect Input spectrum to work on
@param omega Frequency to eval. at
@param x Tan(theta) to eval. at
@return Normalised angular contribution
\caveat We assume omega symmetry here and just take abs(omega). See get_G1 for more

\todo Currently interpolates angle only. Perhaps interpolate on omega too?
*/

  omega = std::abs(omega);//Assume symmetric in omega

  long om_ind, norm_ind, offset;
  size_t len = my_spect->get_omega_length();
  my_type tmpg;
  my_type data_bit[2], ax_bit[2];

  len = my_spect->get_omega_length();

  //If seperable, we have only one (omega) entry in g and norm_g in position 0, so that's where we lookup
  norm_ind = my_spect->get_om_axis_index_from_value(omega);
  if(omega < 0) norm_ind --;//Fix the discrepancy between +ve and -ve omega because where takes first greater than.

  if(norm_ind < 0 ) return 0.0;
  om_ind = (my_spect->get_g_is_angle_only()? 0: norm_ind);

  //Calc norm if hasn't been
  if(norm_ind >= 0 && my_spect->get_norm_g(norm_ind) == 0.0){
    my_spect->calc_norm_g(norm_ind);
    //Already checked norm_ind is non-negative
    if((size_t)norm_ind < len) my_spect->calc_norm_g(norm_ind + 1);
  }
  
  len = my_spect->get_angle_length();
  offset = my_spect->get_ang_axis_index_from_value(x);
  //Interpolate if possible, else use the end
  if(offset > 0 && offset < (long)len){
    data_bit[0] = my_spect->get_g_element(om_ind, offset-1);
    data_bit[1] = my_spect->get_g_element(om_ind, offset);
    ax_bit[0] = my_spect->get_ang_axis_element(offset-1);
    ax_bit[1] = my_spect->get_ang_axis_element(offset);
    tmpg = interpolate_linear(ax_bit, data_bit, (my_type)x);

  }else if(offset == 0){
    //we're right at end, can't meaningfully interpolate, use raw
    tmpg = my_spect->get_g_element(om_ind, offset);

  }else{
    //offset <0 or > len, value not found
    tmpg = 0.0;
  }

  if(offset >= 0 && norm_ind >= 0){
    return tmpg/my_spect->get_norm_g(norm_ind);
  }
  else return 0.0;

}

