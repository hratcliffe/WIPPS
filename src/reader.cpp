//
//  reader.cpp
//  
//
//  Created by Heather Ratcliffe on 02/10/2015.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include "support.h"
#include "reader.h"
#include "my_array.h"


/** This is going to take care of openeing and reading from SDF. We'll give it a pointer to a storage class instance, and it can then use that to obtain pointer to the bit of memory to fill. We'll also give it a numerical range, a file prefix, and a block name to grab. Perhaps we also have option to restrict on the ranges, so we can easily block in space.
*
*It will be sort of safe. It wont store the data array to write to, but will need to be fed it. So we can;t overwrite the memory.
*
*/

reader::reader(std::string file_prefix_in,  char * block_id_in){
/**We'll also have to work out how many zeros to use at some point. For nw we'll do it when we make the reader eh, as that's where we assign the prefix. We'll assume a 0th file exists, and default and minimum is 4. We'll also try 5-7.
*/
  strcpy(this->block_id, block_id_in);
  this->file_prefix = file_prefix_in;
  
  std::ifstream file;
  std::string name = file_prefix + "000.sdf";
  n_z=3;

  while(!file.is_open()){
   name.insert(file_prefix.size(), "0");
   file.open(name);
   n_z ++;
   if(n_z > 8) break;
  }
  file.close();

}

bool reader::read_data(data_array * my_data_in, int time_range[2], int space_range[2]){
/** This will open the files dictated by time range seuqentially, and populate them into the data_array. It'll stop when the end of rnage is reached, or it goes beyond the size available. Space range upper entry of -1 is taken as respective limit. @return 0 for success, 1 for error
*/

//First we check the space range is within range, and chnage the -1's to the repsective dimensions. If out of rnage we return error (1).
  int dim = my_data_in->get_dims(0);
  if((space_range[1] > dim) || ((space_range[1]> 0)  &&(space_range[0] > space_range[1]))) return 1;

  if(space_range[0]==-1) space_range[0] = 0;
  if(space_range[1]==-1) space_range[1] = dim;
  
  std::cout<<space_range[0]<<" "<<space_range[1]<<std::endl;
  //Now we start from time_range[0] and run through files to time_range[1], or until file not found. We construct with sprintf and the known n_z

  char fmt[5];
  char filename[40];
  //enough to hold any sane size of prefix...

  sprintf(fmt,"%s%d%c" , "%0", n_z, 'd');
  snprintf(filename, 40, fmt, 1);
  std::cout<<filename<<std::endl;



return 0;
}

