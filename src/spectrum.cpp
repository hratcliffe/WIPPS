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


spectrum::spectrum(int nx):data_array(nx, 1){
  //Assume by default we''l integrate over second dim...

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

void spectrum::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10]){
//set id params for later...
//times can be normed however we want. Space is relative to grid...
//wave id is defined in header and defines the cutout we use

this->time[0] = time1;
this->time[1] = time2;
this->space[0] = space1;
this->space[1] = space2;
strcpy(this->block_id, block_id);
this->wave_id = wave_id;

}

bool spectrum::generate_spectrum(data_array * parent){
//windows and integrates over frequency by default
//windowing is controlled by the wave_id

if(parent){
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
    om_target = om_disp *(1.0-tolerance);
    //This is a temporary, slow version. Might change to a binary bisect if needed etc...
    for(j=0;j<len; j++) if(ax_ptr[j]> om_target) break;
    low_bnd = j;
    om_target = om_disp *(1.0+tolerance);
    for(j=0;j<len; j++) if(ax_ptr[j]> om_target) break;
    high_bnd = j;
    
    //now total the part of the array between these bnds
    total=0;
    for(j=low_bnd; j<high_bnd; j++) total += parent->get_element(i,j);
    this->set_element(i,j,total);
  }

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



