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
  normg = nullptr;

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
  normg = (my_type *) calloc(1, sizeof(my_type));
  //Single row so only one norm
}

spectrum::spectrum(int nx, int n_ang):data_array(nx, n_ang+1){
/** \brief Construct rectangular spectrum
*
*Constructs a spectrum as a rectangle when for instance the angle dependence varies with k. Spectra are always in the form B^2(omega) g(theta, omega). Here we require to hold both an n_omega x 1 array for B plus the n_omega*n_angs array for g.
*/

  construct();
  angle_is_function = false;
  n_angs = n_ang;
  normg = (my_type *) calloc(n_ang, sizeof(my_type));
  //Norm each row


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

  free(normg);

}

bool spectrum::generate_spectrum(data_array * parent){
/**\brief Generate spectrum from data
*
*Takes a parent data array and uses the specified ids to generate a spectrum. Windows using the specified wave dispersion and integrates over frequency. Also adopts axes from parent. \todo Ensure !angle_is_function forces all other rows to be equal length... \todo Fill in the rest of logic etc
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
    my_type om_disp;
    my_type tolerance = 0.05;
    my_type total;
    for(int i=0; i<this->dims[0]; ++i){
      
      om_disp = get_omega(this->axes[i], WAVE_WHISTLER);

      low_bnd = where(ax_ptr, len, om_disp *(1.0-tolerance));
      high_bnd = where(ax_ptr, len, om_disp *(1.0+tolerance));
      if(low_bnd < 0 || high_bnd< 0){
        this->set_element(i,0,0.0);
        break;
      }
      //now total the part of the array between these bnds
      total=0.0;
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
/** \brief Return g_w(x) for given omega
*
*Returns the row corresponding to omega. If the angle distribution is a single function, the omega param may be omitted. Len is set to the axis length, which is always = n_angs
*/

  my_type * ret = NULL;

  if(angle_is_function){

    ret = data + dims[0];

  }else if(omega !=0.0){
  //select row by omega...
    int offset = where(axes+ dims[0], n_angs, omega);
    if(offset>0) ret = data + offset*dims[0];

  }

  len = n_angs;
  return ret;

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
  if(!file.is_open()) return 1;
  data_array::write_to_file(file);
  return 0;

}

void spectrum::make_test_spectrum(){
/** \brief Generate dummy spectrum
*
*Makes a basic spectrum object with suitable number of points, and twin, symmetric Gaussians centred at fixed x. \todo Finish cases!!!
*/

  char id[10] = "ex";

  this->set_ids(0, 100, 0, get_length(0), WAVE_WHISTLER, id);
  
  ax_omega = false;

  //setup axes
  int len0, len1;
  my_type * ax_ptr;

  ax_ptr = get_axis(1, len1);
  my_type res_x = 4.0/len1;
  for(int i=0; i<len1; i++) *(ax_ptr+i) = (my_type)i*res_x;

  ax_ptr = get_axis(0, len0);
  my_type res_k = 1.0/(my_type)len0;
  for(int i=0; i<len0; i++) *(ax_ptr+i) = res_k*((my_type)i - (my_type)len0/2.0);
  
  //Generate the angle function data.
  if(function_type == FUNCTION_DELTA){
  //Approx delta function, round k_ll. I.e. one cell only. And size is 1/d theta
    for(int i=1; i<len1; ++i) this->set_element(i,1,1);
    //zero all other elements
    my_type val = 1.0/res_x;
    /** \todo this is wrong value. Wants to make integral 1...*/
    this->set_element(0, 1, val);
  }else if(function_type == FUNCTION_GAUSS){
    for(int i=1; i<len1; ++i) this->set_element(i,1,1);


  }else{
    for(int i=1; i<len1; ++i) this->set_element(i,1,1);

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
/** \brief Normalise B(w)
*
* Calculate the total square integral of values over range \todo Is this data bare or squared?
*/

  int len = get_length(0);
  my_type * d_axis = (my_type *) calloc(len, sizeof(my_type));
  for(int i=0; i<len-1; i++) d_axis[i] = get_axis_element(0, i+1) - get_axis_element(0, i);
  normB = integrator(get_ptr(0, 0), len, d_axis);
  
  free(d_axis);
  return 0;
}

bool spectrum::normaliseg(my_type omega){
/** \brief Normalise g_w(x)
*
*Calculate the norm of g used in e.g. denom of Albert eq 3 or calc'd in derivations.tex. We assume omega, x are off the axes already so no interpolation  \todo Catch zero norms
*/

  int len=get_length(1);
  plasma * plas =my_controller->get_plasma();

  my_type * d_axis = (my_type *) calloc(len, sizeof(my_type));
  my_type * integrand = (my_type *) calloc(len, sizeof(my_type));

  for(int i=0; i<len-1; i++) d_axis[i] = get_axis_element(1, i+1) - get_axis_element(1, i);
  //Construct dx axis for integration

  mu my_mu;

  int om_ind = 1;
  int lena=get_length(0);
  if(!angle_is_function){
    if(ax_omega){
      om_ind = where(get_axis(0, lena), lena, omega);
    }else{
      my_type k = get_k(omega, WAVE_WHISTLER);
      om_ind = where(get_axis(0, lena), lena, k);
      
    }
  }
  if(om_ind<0) return 1;
  //break if Omega is out of range

  my_type x, psi;
  int inda, indb;

  //Addressing changes if we have g(x) or g(w, x)
  angle_is_function ? indb=om_ind: inda=om_ind;

  for(int i=0; i<len; i++){
    x = get_axis_element(1, i);
    psi = atan(x);
    my_mu = plas->get_root(0, omega, psi);

    angle_is_function ? inda=i: indb=i;

    if(!my_mu.err){
      integrand[i] = get_element(inda, indb) * x * std::pow((std::pow(x, 2)+1), -1.5)*std::pow(my_mu.mu, 2) * std::abs( my_mu.mu + omega*my_mu.dmudom);
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
\todo Does it matter that our k is limited? Do waves really go to low intensity in bit we see
*/

  calc_type B2;
  my_type tmpB2;
  if(normB ==0.0) normaliseB();

  int len, offset;
  my_type data_bit[2];
  my_type ax_val;
  my_type * axis = this->get_axis(0, len);

  //If we have B(k) we need to change to B(w)
  my_type change_of_vars = 1.0;

  if(ax_omega){
    ax_val = (my_type) omega;
  }else{
    //We have k, need to translate via dispersion relation to get the required index and add the v_g factor
    ax_val = get_k((my_type)omega, WAVE_WHISTLER);
    change_of_vars = get_omega(ax_val, WAVE_WHISTLER, 1);
  }

  offset = where(axis, len, ax_val);

  //Interpolate if possible, else use the end
  if(offset > 0 && offset < len){
    data_bit[0] = get_element(0, offset-1);
    data_bit[1] = get_element(0, offset);
    tmpB2 = interpolate(axis + offset-1, data_bit, ax_val, 2);
  }else if(offset==0){
    //we're right at end, can't meaningfully interpolate, use raw
    tmpB2 = get_element(0, offset);
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
* Gets the value of g(w, x) and the normalising constant from normg \todo Interpolate! \todo IS THIS OMEGA OR do we calc omega according to conditions on integral???
*/


  int om_ind, x_ind, len;
  len=get_length(0);
  if(!angle_is_function){
    if(ax_omega) om_ind = where(get_axis(0, len), len, omega);
    else om_ind = where(get_axis(0, len), len, this->get_k(omega, WAVE_WHISTLER));
  }
  else om_ind = 0;
  if(om_ind>=0 && normg[om_ind] == 0.0){
    normaliseg(omega);
//    std::cout<<"norming"<<om_ind<<std::endl;
  }

  len=get_length(1);
  x_ind = where(get_axis(1, len), len, x);
  
  
  if(x_ind >=0 && om_ind >=0)return get_element(om_ind, x_ind)/normg[om_ind];
  else return 0.0;

}
