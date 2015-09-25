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



//OK now it makes fine sense for spectrum to know about parent class as it's already descended from it. So we can happily have it take one of those and self generate from it so to speak
//or vice versa?


spectrum::spectrum(int nx, data_array* parent):data_array(nx, 1){
  //Assume by default we''l integrate over second dim...
  //parent can be NULL. Caller knows if it isn't we assume. Then call can use first argument as parent->nx if valid. Etc
//if parent not null we borrow the block id and first axis from it. We copy, not refer to same memory as we don't know when parent is destroyed
//And we have no need to keep all the parent data around, so we don't make this part of that
if(parent){
  
  this->block_id = parent->block_id;
  int len;
  my_type * ax_ptr = parent->get_axis(0, len);
  memcpy ((void *)this->axes, (void *)ax_ptr, len*sizeof(my_type));

}

}
