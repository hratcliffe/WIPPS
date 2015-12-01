//
//  spectrum.cpp
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "support.h"
#include "my_array.h"
#include "controller.h"
#include "plasma.h"
#include "spectrum.h"

extern deck_constants my_const;

//OK now it makes fine sense for spectrum to know about parent class as it's already descended from it. So we can happily have it take one of those and self generate from it so to speak
//or vice versa?


void spectrum::construct(){
/** \brief Generic specturm contruction actions
*
*These are performed regardless of dimensions
*/
  function_type = 0;
  wave_id = 0;
  my_controller = nullptr;
  ax_omega = true;
  normB = 0;
}

spectrum::spectrum(int * row_lengths, int ny):data_array(row_lengths, ny){
/** \brief Construct ragged spectrum
*
*Constructs a spectrum as a ragged array, for e.g. two independent functions of omega and angle. Spectra are always in the form B^2(omega) g(theta, omega). Here g(theta, omega) = g(theta). Normally the first row is then the number of omega points, the second the number of angles.
*/
  construct();

  angle_is_function = true;
  n_angs = row_lengths[1];
  function_type = FUNCTION_DELTA;
}

spectrum::spectrum(int nx, int n_ang):data_array(nx, n_ang+1){
/** \brief Construct rectangular spectrum
*
*Constructs a spectrum as a rectangle when for instance the angle dependence varies with k. Spectra are always in the form B^2(omega) g(theta, omega). Here we require to hold both an n_omega x 1 array for B plus the n_omega*n_angs array for g.
*/

  construct();
  angle_is_function = false;
  n_angs = n_ang;


}

void spectrum::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type){
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


}

bool spectrum::generate_spectrum(data_array * parent){
/**\brief Generate spectrum from data
*
*Takes a parent data array and uses the specified ids to generate a spectrum. Windows using the specified wave dispersion and integrates over frequency. Also adopts axes from parent.
*/

  if(parent && angle_is_function){
    //First we read axes from parent
    int len;
    ax_omega = false;

    my_type * ax_ptr = parent->get_axis(0, len);
    memcpy ((void *)this->axes, (void *)ax_ptr, len*sizeof(my_type));
    ax_ptr = parent->get_axis(1, len);
    //y-axis to work with

    //Now we loop across x, calculate the wave cutout bounds, and total, putting result into data

    int j;
    int low_bnd, high_bnd;
    float om_disp;
    float tolerance = 0.05;
    float total;
    for(int i=0; i<this->dims[0]; ++i){
      
      om_disp = get_omega(this->axes[i], WAVE_WHISTLER);

      low_bnd = where(ax_ptr, len, om_disp *(1.0-tolerance));
      high_bnd = where(ax_ptr, len, om_disp *(1.0+tolerance));
      
      //now total the part of the array between these bnds
      total=0;
      for(j=low_bnd; j<high_bnd; j++) total += parent->get_element(i,j);
      this->set_element(i,0,total);
    }

    //Now we generate evenly spaced angle axis, and generate required function ...
    //NB What to work in? tan theta, but from 0 to infty. Need a cut off.... and that only covers paralllel, not antiparallel. We'll need to extend to that sooner or later...
    {
      int res = 1;
      //set axis resolution somehow... TODO this
      make_linear_axis(1, res, 0);

      //Now generate the function data.
      if(function_type == FUNCTION_DELTA){
      //Approx delta function, round k_ll. I.e. one cell only. And size is 1/d theta
      
        for(int i=1; i<this->dims[0]; ++i) this->set_element(i,1,0);
        //zero all other elements
        float val;
        val = 1.0/res;
        //TODO this is wrong value. Wants to make integral 1...
        this->set_element(0, 1, val);
      }else if(function_type == FUNCTION_GAUSS){


      }else if(function_type == FUNCTION_DELTA){


      }else{

      }

    }

  }else if(parent){
 /** \todo general spectrum extracttion routine */
  //TODO in this case we have to extract spectrim and angle data somehow......

    //First we read axes from parent
    int len;
    ax_omega = false;

    my_type * ax_ptr = parent->get_axis(0, len);
    memcpy ((void *)this->axes, (void *)ax_ptr, len*sizeof(my_type));
    ax_ptr = parent->get_axis(1, len);
    //y-axis to work with

    //Generate angle axis to work with
  {int res = 1;
  //set axis resolution somehow... TODO this
  make_linear_axis(1, res, 0);
  }

  //and now we extract the data at each angle...
  //tan theta = k_x/k_y



  }else{
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

my_type * spectrum::get_angle_distrib(int &len, my_type omega){
//Now if it's a single function we return just a row, and dont need the omega param
//If it's not, we select row by omega. The length will always be n_angs, but we return it for clarity

  my_type * ret = NULL;

  if(angle_is_function){

    ret = data + dims[0];

  }else if(omega !=0.0){
  //select row by omega...
    int offset = where(axes+ dims[0], n_angs, omega);
    ret = data + offset*dims[0];

  }

  len = n_angs;
  return ret;

}

int spectrum::where(my_type * ax_ptr, int len, my_type target,std::function<bool(my_type,my_type)> func){
/**\brief Finds first index where ax_ptr[i] vs target value satisfies the function given (e.g. std::greater)
*
*
*/
  int j;
//  for(j=0;j<len; j++) if(ax_ptr[j]> target) break;
  for(j=0;j<len; j++) if(func(ax_ptr[j],  target)) break;

  return j;

}

std::vector<int> spectrum::all_where(my_type * ax_ptr, int len, my_type target,std::function<bool(my_type,my_type)> func){
/**\brief Finds all indices where ax_ptr[i] vs target value satisfies the function given (e.g. std::greater)
*
*
*/
  std::vector<int> ret;

  for(int j=0;j<len; j++) if(func(ax_ptr[j],  target)) ret.push_back(j);

  return ret;

}

bool spectrum::write_to_file(std::fstream &file){

/**IMPORTANT: the VERSION specifier links output files to code. If modifying output or order commit and clean build before using.
*/
  if(!file.is_open()) return 1;
  data_array::write_to_file(file);
  return 0;

}

void spectrum::make_test_spectrum(){
/**Makes a basic spectrum object with suitable number of points, and twin, symmetric Gaussians centred at fixed x.
*/

  char id[10] = "ex";

  this->set_ids(0, 100, 0, dims[0], WAVE_WHISTLER, id);
  
  ax_omega = false;

  //setup axes
  int len0, len1;
  my_type * ax_ptr;

  ax_ptr = get_axis(1, len1);
  for(int i=0; i<len1; i++) *(ax_ptr+i) = (float)i;

  ax_ptr = get_axis(0, len0);
  float res = 1.0/(float)len0;

  for(int i=0; i<len0; i++) *(ax_ptr+i) = res*((float)i - (float)len0/2.0);

  
  //Generate the angle function data.
  if(function_type == FUNCTION_DELTA){
  //Approx delta function, round k_ll. I.e. one cell only. And size is 1/d theta
    float res = 1.0;
    for(int i=1; i<this->row_lengths[1]; ++i) this->set_element(i,1,10);
    //zero all other elements
    float val = 1.0/res;
    /** \todo this is wrong value. Wants to make integral 1...*/
    this->set_element(0, 1, val);
  }else if(function_type == FUNCTION_GAUSS){


  }else{

  }
  //Generate the negative k data
  
  float centre = 0.2, width=0.005, background = 0.5;
  my_type * data_ptr = get_ptr(0, 0);
  my_type * data_tmp, *ax_tmp;
  data_tmp = data_ptr;
  ax_tmp = ax_ptr;

  for(int i=0; i<=len0/2; i++, ax_tmp++, data_tmp++) *(data_tmp) = exp(-pow((*(ax_tmp) + centre), 2)/width) + background;
  data_tmp--;
  //we've gone one past our termination condition...
  for(int i=1; i<len0/2; i++) *(data_tmp + i) = *(data_tmp - i);
  //reflect onto +ve k
//+ 0.25*exp(-pow((*(ax_tmp) + centre), 2)/width*50.0)

}

bool spectrum::normaliseB(){
/** Calculate the total square integral of values over range \todo Is this data bare or squared?*/
//calc_type integrator(calc_type * start, int len, calc_type * increment){

  my_type * d_axis = (my_type *) calloc(row_lengths[0], sizeof(my_type));
  for(int i=0; i<row_lengths[0]-1; i++) d_axis[i] = get_axis_element(0, i+1) - get_axis_element(0, i);

  normB = integrator(get_ptr(0, 0), row_lengths[0], d_axis);

  return 0;
}

calc_type spectrum::get_G1(calc_type omega){
/**returns G1 calculated as in Albert 2005.
\todo Does it matter that our k is limited? Do waves really go to low intensity in bit we see
*/

  calc_type B2;
  my_type tmpB2;
  if(normB ==0.0) normaliseB();

  if(ax_omega){
    //look up omega and interpolate
    int len;
    my_type data_bit[2];
    my_type omega_tmp = (my_type) omega;
    my_type * axis = this->get_axis(0, len);
    //get the last place where axis is not larger than omega
    int offset = where(axis, len, omega_tmp);
    if(offset > 0 && offset < len){
      data_bit[0] = get_element(0, offset-1);
      data_bit[1] = get_element(0, offset);
      //Get interpolated value of B for this k
      tmpB2 = interpolate(axis + offset-1, data_bit, omega_tmp, 2);
    }else{
      //we're right at end, can't meaningfully interpolate, use raw
      tmpB2 = get_element(0, offset);
    }
    //Cast to type and multiply vg to finish change vars
    B2 = (calc_type) tmpB2;
    
  }else{
    //We have k, need to translate via dispersion relation to get the required index and add the v_g factor
    my_type k = get_k((my_type)omega, WAVE_WHISTLER);
    //find this k in axis
    int len;
    my_type data_bit[2];
    my_type * axis = this->get_axis(0, len);
    //get the last place where axis is not larger than k
    int offset = where(axis, len, k);
    if(offset > 0 && offset < len){
      data_bit[0] = get_element(0, offset-1);
      data_bit[1] = get_element(0, offset);
      //Get interpolated value of B for this k
      tmpB2 = interpolate(axis + offset-1, data_bit, k, 2);
    }else{
      //we're right at end, can't meaningfully interpolate, use raw
      tmpB2 = get_element(0, offset);
    }
    //Cast to type and multiply vg to finish change vars
    B2 = (calc_type) tmpB2 * get_omega(k, WAVE_WHISTLER, 1);
  }

  return B2/normB;

}

calc_type spectrum::get_G2(calc_type omega, mu_dmudom my_mu){
/**returns G2 calculated as in Albert 2005. */

  calc_type a;
  a = my_mu.mu;
  //to silence unused warning temporarily
  return 0.0;
}
