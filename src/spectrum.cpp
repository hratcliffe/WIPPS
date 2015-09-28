//
//  spectrum.cpp
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#include <stdio.h>
#include "../include/support.h"
#include "../include/my_array.h"
#include "../include/spectrum.h"


extern deck_constants my_const;

//OK now it makes fine sense for spectrum to know about parent class as it's already descended from it. So we can happily have it take one of those and self generate from it so to speak
//or vice versa?


spectrum::spectrum(int nx):data_array(nx+DEFAULT_N_ANG, 1){
  //Assume by default we''l integrate over second dim...

angle_is_function = true;
//angle profile will be specified as a function and any val is product spect(omega, theta) = B^2(omega) * g(theta). So we only need 2 one-d arrays
n_angs = DEFAULT_N_ANG;
function_type = FUNCTION_DELTA;

}

spectrum::spectrum(int nx, int n_ang):data_array(nx, n_ang+1){
  //Assume by default we''l integrate over second dim...

angle_is_function = false;

//angle profile is array, spect(omega, theta) = B^2(omega) * g(theta, omega). So we need a full nx by n_ang array for g, and an extra row for B^2
n_angs = n_ang;


}

spectrum::spectrum(int nx, data_array* parent):data_array(nx, 1){
  //Assume by default we''l integrate over second dim...
  //parent can be NULL. Caller knows if it isn't we assume. Then call can use first argument as parent->nx if valid. Etc
//if parent not null we borrow the block id and first axis from it. We copy, not refer to same memory as we don't know when parent is destroyed.
//And we have no need to keep all the parent data around, so we don't make this part of that

if(parent){
  strcpy(this->block_id, block_id);
  
  int len;
  my_type * ax_ptr = parent->get_axis(0, len);
  memcpy ((void *)this->axes, (void *)ax_ptr, len*sizeof(my_type));

}

}

void spectrum::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type){
//set id params for later...
//times can be normed however we want. Space is relative to grid...
//wave id is defined in header and defines the cutout we use

this->time[0] = time1;
this->time[1] = time2;
this->space[0] = space1;
this->space[1] = space2;
strcpy(this->block_id, block_id);
this->wave_id = wave_id;
this->function_type = function_type;
}

bool spectrum::generate_spectrum(data_array * parent){
//windows and integrates over frequency by default
//windowing is controlled by the wave_id

if(parent && angle_is_function){
  //First we read axes from parent
  int len;
  my_type * ax_ptr = parent->get_axis(0, len);
  memcpy ((void *)this->axes, (void *)ax_ptr, len*sizeof(my_type));
  ax_ptr = parent->get_axis(1, len);
  //y-axis to work with

  //Now we loop across x, calculate the wave cutout bounds, and total, putting result into data

  int j;
  int low_bnd, high_bnd;
  float om_disp, om_target;
  float tolerance = 0.05;
  float total;
  for(int i=0; i<this->dims[0]; ++i){
    
    om_disp = get_dispersion(this->axes[i], WAVE_WHISTLER);
//    om_target = om_disp *(1.0-tolerance);
    //This is a temporary, slow version. Might change to a binary bisect if needed etc...
 //   for(j=0;j<len; j++) if(ax_ptr[j]> om_target) break;
//    low_bnd = j;
    low_bnd = where(ax_ptr, len, om_disp *(1.0-tolerance));
    high_bnd = where(ax_ptr, len, om_disp *(1.0+tolerance));
    
  //  om_target = om_disp *(1.0+tolerance);
    //for(j=0;j<len; j++) if(ax_ptr[j]> om_target) break;
    //high_bnd = j;
    
    //now total the part of the array between these bnds
    total=0;
    for(j=low_bnd; j<high_bnd; j++) total += parent->get_element(i,j);
    this->set_element(i,0,total);
  }

//Now we generate evenly spaced angle axis, and generate required function ...
//NB What to work in? tan theta, but from 0 to infty. Need a cut off.... and that only covers paralllel, not antiparallel. We'll need to extend to that sooner or later...


//void data_array::make_linear_axis(int dim, float res, int offset){
{int res = 1;
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

//TODO in this case we have to extract spectrim and angle data somehow......

  //First we read axes from parent
  int len;
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

float spectrum::get_dispersion(my_type k, int wave_type){

float ret=0;
float tmp;
switch(wave_type){

  case WAVE_WHISTLER :
//    om_whis = v0^2*k_ax^2*om_ce/(v0^2*k_ax^2+om_pe^2)/om_ce
    ret = v0*v0*k*k*my_const.omega_ce/(v0*v0 + my_const.omega_pe*my_const.omega_pe)/my_const.omega_ce;
    
    //here goes dispersion in suitable normed units.
    break;

  case WAVE_PLASMA :
    ret =1 ;
    break;
}

return ret;

}

my_type * spectrum::get_angle_distrib(my_type ang, my_type omega){
//Now if it's a single function we return just a row, and dont need the omega param
//If it's not, we select row by omega

my_type * ret = NULL;

if(angle_is_function){
  ret = data + dims[0];

}else{

//select row by omega...
  int offset = where(axes+ dims[0], n_angs, omega);
  ret = data + offset*dims[0];

}

return ret;

}

int spectrum::where(my_type * ax_ptr, int len, my_type target){
//method so we can upgrade easily...

  int j;
  for(j=0;j<len; j++) if(ax_ptr[j]> target) break;
  
  return j;

}
